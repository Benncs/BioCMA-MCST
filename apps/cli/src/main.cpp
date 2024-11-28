#include "udf_includes.hpp"
#include <api/api.hpp>
#include <cli_parser.hpp>
#include <common/common.hpp>
#include <common/execinfo.hpp>
#include <core/case_data.hpp>
#include <core/simulation_parameters.hpp>
#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <rt_init.hpp>
#include <stream_io.hpp>

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
 * @brief Prepares the case data based on execution information and simulation
 * parameters.
 *
 * This function processes the provided execution information and simulation
 * parameters to prepare and return a `CaseData` object. It is designed to be
 * called before the main simulation execution to ensure all necessary data is
 * correctly initialized.
 *
 * @param exec_info The execution information containing details about the
 * current simulation run, such as environment settings and execution context.
 * @param params    The simulation parameters that configure various aspects of
 * the simulation, including time steps, boundary conditions, and other
 * relevant settings.
 * @return A `CaseData` object containing the prepared data for the simulation.
 */
static std::optional<Core::CaseData> prepare(const ExecInfo& exec_info,
                                             const Core::UserControlParameters&& param);

/**
 * @brief Check if result path exist or not and ask for overriding if yes
 * @return true if override results_path
 */
static bool override_result_path(const Core::UserControlParameters& params, const ExecInfo& exec);

#define HANDLE_RC(__api_results__)                                                                 \
  {                                                                                                \
    auto rc = (__api_results__);                                                                   \
    if (!rc)                                                                                       \
    {                                                                                              \
      std::cout << "ERROR " << #__api_results__ << " " << rc.get() << std::endl;                   \
      return -1;                                                                                   \
    }                                                                                              \
  }


#include "udf_includes.hpp"


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
  auto  _ = UnsafeUDF::Loader::init_lib("/home/benjamin/Documents/code/cpp/BioCMA-MCST/builddir/debug_python/apps/udf_model/libudf_model.so");

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

  auto handle = Handle::init(
      exec_info.n_rank, 
      exec_info.current_rank, 
      exec_info.run_id, 
      exec_info.thread_per_process);

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
