#include "core/global_initaliser.hpp"
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
 * the simulation, including time steps, boundary conditions, and other relevant
 * settings.
 * @return A `CaseData` object containing the prepared data for the simulation.
 */
static std::optional<Core::CaseData> prepare(const ExecInfo &exec_info, const Core::UserControlParameters &&param);

/**
 * @brief Wrapper to handle Excception raised in try/catch block
 */
template <typename ExceptionType> static int handle_catch(ExceptionType const &e) noexcept;

/**
 * @brief Check if result path exist or not and ask for overriding if yes
 * @return true if override results_path
 */
static bool override_result_path(const Core::UserControlParameters &params, const ExecInfo &exec);

int main(int argc, char **argv)
{
  // First manually retrieve argument from command line
  auto params_opt = parse_cli(argc, argv);
  if (!params_opt.has_value())
  {
    // If needed value are not given or invalid argument early return
    showHelp(std::cout);
    return -1;
  }

  auto user_params = params_opt.value(); // Deref value is safe  TODO: with c++23 support use monadic

  /*Init environnement (MPI+Kokkos)
    Note that environnement should be the first action in the code to avoid conflict*/
  const ExecInfo exec_info = runtime_init(argc, argv, user_params);

  // Ask overring results file
  if (!override_result_path(user_params, exec_info))
  {
    return -1;
  }
  std::optional<std::unique_ptr<Handle>> handle =
      Handle::init(exec_info.n_rank, exec_info.current_rank, exec_info.run_id, exec_info.thread_per_process);

  if (handle.has_value())
  {
    auto &h = *handle;
    h->register_parameters(std::move(user_params));
    h->apply(user_params.serde);
    try
    {
      INTERPRETER_INIT
      REDIRECT_SCOPE({ h->exec(); })
    }
    catch (std::exception const &e)
    {
      return handle_catch(e);
    }
    catch (...)
    {
      std::cerr << "Internal error" << std::endl;
      return -1;
    }
  }

  // const auto f_get_case_data = (user_params.serde) ? Core::load : prepare;

  // /*Main loop*/
  // try
  // {
  //   INTERPRETER_INIT

  //   REDIRECT_SCOPE({
  //     auto case_data = f_get_case_data(exec_info, std::move(user_params));
  //     if (case_data)
  //     {
  //       exec(std::move(*case_data));
  //     }
  //     else
  //     {
  //       std::cout << "ERROR" << std::endl;
  //       return -1;
  //     }
  //   })

  // }
  // catch (std::exception const &e)
  // {
  //   return handle_catch(e);
  // }
  // catch (...)
  // {
  //   std::cerr << "Internal error" << std::endl;
  //   return -1;
  // }

  return 0;
}

static std::optional<Core::CaseData> prepare(const ExecInfo &exec_info, const Core::UserControlParameters &&user_params)
{

  Core::GlobalInitialiser gi(exec_info, user_params);
  auto t = gi.init_transitionner();

  auto __simulation = gi.init_simulation();
  if ((!t.has_value() && !__simulation.has_value()) || !gi.check_init_terminate())
  {
    return std::nullopt;
  }

  Core::SimulationParameters params = gi.get_parameters();

  register_run(exec_info, user_params);

  return std::make_optional<Core::CaseData>(std::move(*__simulation), params, std::move(*t), exec_info);
}

template <typename ExceptionType> static int handle_catch(ExceptionType const &e) noexcept
{
#ifdef DEBUG
  std::cerr << "Caught Exception: " << e.what() << '\n';
  return -1;
#else
  std::cerr << "Internal error" << '\n';
  return -1;
#endif
}

bool override_result_path(const Core::UserControlParameters &params, const ExecInfo &exec)
{
  bool flag = true;
  if (exec.current_rank == 0)
  {
    if (std::filesystem::exists(params.results_file_name + std::string(".h5")) && !params.force_override)
    {
      std::cout << "Override results ? (y/n)" << std::endl;
      std::string res;
      std::cin >> res;
      flag = res == "y";
    }
  }
#ifndef NO_MPI
  MPI_W::barrier();
  MPI_W::broadcast(flag, 0);
#endif
  return flag;
}
