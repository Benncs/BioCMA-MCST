#include "common/execinfo.hpp"
#include "common/simulation_parameters.hpp"

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

struct CaseData
{
  std::unique_ptr<Simulation::SimulationUnit> simulation;
  SimulationParameters params;
  std::unique_ptr<Simulation::FlowMapTransitioner> transitioner;
  ExecInfo exec_info;
};

constexpr bool verbose = true;   // TODO REMOVE
constexpr bool redirect = false; // TODO REMOVE

static CaseData prepare(const ExecInfo &exec_info, SimulationParameters params);

static void exec(CaseData &&case_data);

template <typename ExceptionType>
static int handle_catch(ExceptionType const &e);

int main(int argc, char **argv)
{

  init_environment();
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

    REDIRECT_BLOCK(
        {
          auto case_data = prepare(exec_info, params);
          exec(std::move(case_data));
        },
        verbose,
        redirect)
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