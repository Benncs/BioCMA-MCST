#include <common/execinfo.hpp>
#include <common/simulation_parameters.hpp>
#include <cli_parser.hpp>
#include <common/common.hpp>
#include <simulation/simulation.hpp>
#include <simulation/update_flows.hpp>
#include <cma_read/flow_iterator.hpp>
#include <cma_read/reactorstate.hpp>
#include <host_specific.hpp>
#include <mpi_w/wrap_mpi.hpp>
#include <post_process.hpp>
#include <rt_init.hpp>
#include <siminit.hpp>
#include <sync.hpp>
#include <worker_specific.hpp>

#ifdef USE_PYTHON_MODULE
#  include <pymodule/import_py.hpp>
#endif

#include <iostream>
#include <memory>
#include <stream_io.hpp>

#ifdef USE_PYTHON_MODULE
#  define INTERPRETER_INIT auto _interpreter_handle = init_python_interpreter();
#else
#  define INTERPRETER_INIT
#endif

/**
 * @brief Holds the data required to execute a simulation case.
 *
 * The `CaseData` struct encapsulates all the necessary components needed to
 * perform a simulation. It includes the simulation unit, parameters,
 * transitioner, and execution information. This structure is typically used to
 * manage and transfer simulation data between different stages or components of
 * the simulation process.
 */
struct CaseData
{
  /**
   * @brief Unique pointer to the simulation unit.
   *
   * This is the core unit of the simulation, responsible for executing the
   * main simulation logic.
   */
  std::unique_ptr<Simulation::SimulationUnit> simulation;

  /**
   * @brief Parameters that configure the simulation.
   *
   * The simulation parameters control various aspects of the simulation,
   * such as time steps, boundary conditions, and solver configurations.
   */
  SimulationParameters params;

  /**
   * @brief Unique pointer to the flow map transitioner.
   *
   * The transitioner manages transitions between different flow states within
   * the simulation, ensuring consistency and accuracy in the simulation's
   * progression.
   */
  std::unique_ptr<Simulation::FlowMapTransitioner> transitioner;

  /**
   * @brief Information about the execution environment.
   *
   * The execution information contains details such as the execution context,
   * environment settings, and other relevant metadata necessary for the
   * simulation run.
   */
  ExecInfo exec_info;
};

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
static CaseData prepare(const ExecInfo &exec_info, SimulationParameters params);

/**
 * @brief Start simulation
 */
static void exec(CaseData &&case_data);

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

static CaseData prepare(const ExecInfo &exec_info, SimulationParameters params)
{

  std::unique_ptr<Simulation::FlowMapTransitioner> transitioner = nullptr;

  auto simulation = init_simulation(exec_info, params, transitioner);
  return {std::move(simulation), params, std::move(transitioner), exec_info};

}

static void exec(CaseData &&case_data)
{
  const auto f_run = (case_data.exec_info.current_rank == 0) ? &host_process
                                                             : &workers_process;

  auto *const sim = case_data.simulation.release();

  f_run(case_data.exec_info,
        std::move(*sim),
        case_data.params,
        std::move(case_data.transitioner));
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