#include <cli_parser.hpp>
#include <common/common.hpp>
#include <common/execinfo.hpp>
#include <common/simulation_parameters.hpp>
#include <core/case_data.hpp>
#include <core/cp_flag.hpp>
#include <iostream>
#include <memory>
#include <rt_init.hpp>
#include <siminit.hpp>
#include <stream_io.hpp>


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
static Core::CaseData prepare(const ExecInfo &exec_info,
                              SimulationParameters params);

/**
 * @brief Wrapper to handle Excception raised in try/catch block
 */
template <typename ExceptionType>
static int handle_catch(ExceptionType const &e);

int main(int argc, char **argv)
{

  auto params_opt = parse_cli(argc, argv);
  if (!params_opt.has_value())
  {
    showHelp(std::cout);
    return -1;
  }

  const auto params = params_opt.value();
  const ExecInfo exec_info = runtime_init(argc, argv, params);
  try
  {

    INTERPRETER_INIT

    REDIRECT_BLOCK({
      auto case_data = prepare(exec_info, params);
      exec(std::move(case_data));
    })
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

  return 0;
}

static Core::CaseData prepare(const ExecInfo &exec_info,
                              SimulationParameters params)
{

  std::unique_ptr<Simulation::FlowMapTransitioner> transitioner = nullptr;
  auto simulation = init_simulation(exec_info, params, transitioner);
  return {std::move(simulation), params, std::move(transitioner), exec_info};
}

template <typename ExceptionType>
static int handle_catch(ExceptionType const &e)
{
#ifdef DEBUG
  std::cerr << e.what() << '\n';
  return -1;
#else
  std::cerr << "Internal error" << '\n';
  return -1;
#endif
}