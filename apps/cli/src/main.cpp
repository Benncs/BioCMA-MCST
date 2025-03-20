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
#ifdef DECLARE_EXPORT_UDF
#  include <udf_includes.hpp>
#  define DECLARE_LOADER(__path__) auto _ = UnsafeUDF::Loader::init_lib(__path__)
#else

#  define DECLARE_LOADER(__path__) (void)__path__;
#endif
#include <utility>

#ifndef NO_MPI
#  include <mpi_w/wrap_mpi.hpp>
#endif

#ifdef USE_PYTHON_MODULE
#  include <pymodule/import_py.hpp>
#  define INTERPRETER_INIT auto _interpreter_handle = init_python_interpreter();
#else
#  define INTERPRETER_INIT
#endif

/**
 * @brief Check if result path exist or not and ask for overriding if yes
 * @return true if override results_path
 */
static bool override_result_path(const Core::UserControlParameters& params, const ExecInfo& exec);

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

int main(int argc, char** argv)
{
  // First manually retrieve argument from command line
  auto params_opt = parse_cli(argc, argv);
  if (!params_opt.has_value())
  {
    // If needed value are not given or invalid argument early return
    showHelp(std::cout);
    return -1;
  }

  // DECLARE_LOADER("/home-local/casale/Documents/thesis/code/BioCMA-MCST/builddir/test/apps/"
  //                "udf_model/libudf_model.so");

  auto user_params = params_opt.value(); // Deref value is safe  TODO: with
                                         // c++23 support use monadic

  /*Init environnement (MPI+Kokkos)
    Note that environnement should be the first action in the code to avoid
    conflict*/
  const ExecInfo exec_info = runtime_init(argc, argv, user_params);

  // Ask overring results file
  if (!override_result_path(user_params, exec_info))
  {
    return -1;
  }

  auto handle = Api::SimulationInstance::init(
      exec_info.n_rank, exec_info.current_rank, exec_info.run_id, exec_info.thread_per_process);

  if (!handle)
  {
    std::cerr << "Error Handle init" << std::endl;
    return -1;
  }

  auto& h = *handle;
  const auto serde = user_params.serde;
  INTERPRETER_INIT
  REDIRECT_SCOPE({
    HANDLE_RC(h->register_parameters(std::move(user_params)));
    // 1:20e-3*0.5/3600., {3}
    // 2:20e-3*0.9/3600., {10}
    h->set_feed_constant_from_rvalue(20e-3 * 0.5 / 3600., {300e-3}, {0}, {1}, true);

    // Sanofi
    //  h->set_feed_constant_from_rvalue(0.031653119013143756, {0.}, {0}, {0});
    // h->set_feed_constant_from_rvalue(0.001 / 3600., {0.3}, {0}, {1}, true);
    // h->set_feed_constant_from_rvalue(20e-3*0.4/3600., {10}, {0}, {0}, false);
    // h->set_feed_constant_from_rvalue(100*0.01 / 3600., {0.3}, {0}, {1}, true);

    // std::vector<double> c(100);
    // std::vector<std::size_t> p(100);
    // std::vector<std::size_t> ss(100);
    // for (int i = 0; i < 100; ++i)
    // {
    //   p[i] = i;
    //   c[i] = 0.3;
    //   ss[i] = i;
    // }

    // h->set_feed_constant(0.031653119013143756, c, p, ss, true);

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
