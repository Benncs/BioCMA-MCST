
#include "simulation/feed_descriptor.hpp"
#include <api/api.hpp>
#include <cli_parser.hpp>
#include <common/execinfo.hpp>
#include <core/case_data.hpp>
#include <core/simulation_parameters.hpp>
#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <rt_init.hpp>
#include <stream_io.hpp>
#include <string>
#include <string_view>
#include <utility>

#ifdef DECLARE_EXPORT_UDF
#  include <udf_includes.hpp>
std::shared_ptr<DynamicLibrary> wrap_non_scoped_udf(std::string_view path, bool load)
{
  if (load)
  {
    return UnsafeUDF::Loader::init_lib(path);
  }
  else
  {
    return nullptr;
  }
}
#  define DECLARE_LOADER(__path__) auto _ = wrap_non_scoped_udf(__path__, true);
#else

#  define DECLARE_LOADER(__path__) (void)__path__;
#endif

#ifndef NO_MPI
#  include <mpi_w/wrap_mpi.hpp>
#endif

#ifdef USE_PYTHON_MODULE
#  include <pymodule/import_py.hpp>
#  define INTERPRETER_INIT auto _interpreter_handle = init_python_interpreter();
#else
#  define INTERPRETER_INIT
#endif

#ifndef NO_MPI
#  define HANDLE_RC(__api_results__)                                                               \
    {                                                                                              \
      auto rc = (__api_results__);                                                                 \
      if (!rc)                                                                                     \
      {                                                                                            \
        std::cout << "ERROR " << #__api_results__ << " " << rc.get() << std::endl;                 \
        WrapMPI::critical_error();                                                                 \
        return -1;                                                                                 \
      }                                                                                            \
    }
#else
#  define HANDLE_RC(__api_results__)                                                               \
    {                                                                                              \
      auto rc = (__api_results__);                                                                 \
      if (!rc)                                                                                     \
      {                                                                                            \
        std::cout << "ERROR " << #__api_results__ << " " << rc.get() << std::endl;                 \
        return -1;                                                                                 \
      }                                                                                            \
    }
#endif

/**
 * @brief Check if result path exist or not and ask for overriding if yes
 * @return true if override results_path
 */
static bool override_result_path(const Core::UserControlParameters& params, const ExecInfo& exec);
static int parse_callback_ok(Core::UserControlParameters&& user_params,
                             std::optional<std::unique_ptr<Api::SimulationInstance>>& handle);

int main(int argc, char** argv)
{

  auto handle = Api::SimulationInstance::init(argc, argv);

  if (!handle)
  {
    std::cerr << "Error Handle init" << std::endl;
    return -1;
  }

  return parse_cli(argc, argv)
      .match(
          [&](auto&& user_params)
          { return parse_callback_ok(std::forward<decltype(user_params)>(user_params), handle); },
          [](auto&& val)
          {
            std::cout << "Err: " << val << std::endl;
            showHelp(std::cout);
            return 1;
          });
}

int parse_callback_ok(Core::UserControlParameters&& user_params,
                      std::optional<std::unique_ptr<Api::SimulationInstance>>& handle)
{
  DECLARE_LOADER("/home-local/casale/Documents/code/poc/builddir/host/apps/udf_model/"
                 "libudf_model.so");

  auto& h = *handle;
  if (!override_result_path(user_params, h->get_exec_info()))
  {
    return -1;
  }

  // auto sine_feed = Simulation::Feed::FeedFactory::pulse(20e-3 * 5 / 3600.,{5}, {0}, 
  // {0},0,0,2./3600.,true);

  const auto serde = user_params.serde;
  INTERPRETER_INIT
  REDIRECT_SCOPE({
    HANDLE_RC(h->register_parameters(std::forward<decltype(user_params)>(user_params)));
    h->set_feed_constant_from_rvalue(20e-3 * 0.5 / 3600., {300e-3}, {0}, {1}, true);
    // h->set_feed(sine_feed);
    // h->set_feed_constant_from_rvalue(0.031653119013143756, {0.}, {0}, {0},false);
    // h->set_feed_constant_from_rvalue(20e-3 * 0.8 / 3600., {5}, {0}, {0}, false);
    //   h->set_feed_constant_from_rvalue(20e-3 * 0.8 / 3600., {8e-3}, {0}, {1}, false);
    HANDLE_RC(h->apply(serde));
    HANDLE_RC(h->exec());
  })
  return 0;
}

bool override_result_path(const Core::UserControlParameters& params, const ExecInfo& exec)
{
  bool flag = true;
  if (exec.current_rank == 0)
  {
    if (std::filesystem::exists(params.results_file_name + std::string(".h5")) &&
        !params.force_override)
    {
      std::cout << "Override results ? (y/n)" << std::endl;
      std::string res;
      std::cin >> res;
      flag = res == "y";
    }
  }
#ifndef NO_MPI
  WrapMPI::barrier();
  WrapMPI::broadcast(flag, 0);
#endif
  return flag;
}

//   // h->set_feed_constant_from_rvalue(20e-3 * 0.5 / 3600., {300e-3}, {0}, {1}, true);

//   // Sanofi
//   //  h->set_feed_constant_from_rvalue(0.031653119013143756, {0.}, {0}, {0});
//   // h->set_feed_constant_from_rvalue(0.001 / 3600., {0.3}, {0}, {1}, true);
//   // h->set_feed_constant_from_rvalue(20e-3*0.9/3600., {10}, {0}, {0}, false);
//   // h->set_feed_constant_from_rvalue(100*0.01 / 3600., {0.3}, {0}, {1}, true);
