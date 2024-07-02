#include "common/execinfo.hpp"
#include "common/simulation_parameters.hpp"
#include "mpi.h"
#include "mpi_w/impl_op.hpp"
#include "mpi_w/message_t.hpp"
#include <cli_parser.hpp>

#include <common/common.hpp>
#include <models/models.hpp>
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

#include <model_list.hpp>


#ifdef USE_PYTHON_MODULE
#  include <pymodule/import_py.hpp>
#endif

#include <cstddef>
#include <exception>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>

#ifdef USE_PYTHON_MODULE
#define INTERPRETER_INIT auto _interpreter_handle = init_python_interpreter(); 
#else
#define INTERPRETER_INIT
#endif 

struct CaseData
{
  std::unique_ptr<Simulation::SimulationUnit> simulation;
  SimulationParameters params;
  std::unique_ptr<Simulation::FlowMapTransitioner> transitioner;
  ExecInfo exec_info;
};

/*The following should be remove or moved */
constexpr bool verbose = true;   // TODO REMOVE
constexpr bool redirect = false; // TODO REMOVE
static bool is_stdout_redirect = false;

static int redirect_stdout(std::streambuf *&original_buffer,
                           std::stringstream &variable_stream);
static void restore_stdout(int original_stdout,
                           std::streambuf *&original_buffer);

static void
workers_process(const ExecInfo &exec,
                Simulation::SimulationUnit &simulation,
                SimulationParameters &params,
                std::unique_ptr<Simulation::FlowMapTransitioner> transitioner);


static CaseData prepare(const ExecInfo &exec_info, SimulationParameters params);

static void exec(const ExecInfo &exec_info, CaseData &&cased);

int main(int argc, char **argv)
{
  init_environment();
  auto params_opt = parse_cli(argc, argv);
  if (!params_opt.has_value())
  {
    showHelp(std::cout);
    return -1;
  }
  auto params = params_opt.value();
  const ExecInfo exec_info = runtime_init(argc, argv, params);
  try
  {

    INTERPRETER_INIT

    // TODO REMOVE
    if constexpr (redirect)
    {
      std::stringstream output_variable;
      output_variable.str("");
      std::streambuf *original_buffer = nullptr;
      auto original_stdout_fd =
          redirect_stdout(original_buffer, output_variable);

      auto case_data = prepare(exec_info, std::move(params));
      exec(exec_info, std::move(case_data));

      restore_stdout(original_stdout_fd, original_buffer);
      if constexpr (verbose)
      {
        std::cout << output_variable.str() << std::endl;
      }
    }
    else
    {
      auto case_data = prepare(exec_info, std::move(params));
      exec(exec_info, std::move(case_data));
    }
  }
#ifdef DEBUG
  catch (std::exception const &e)

  {
    std::cerr << e.what() << '\n';
    return -1;
  }
#endif
  catch (...)
  {
    std::cerr << "Internal error" << '\n';
    return -1;
  }

  return 0;
}

static CaseData prepare(const ExecInfo &exec_info, SimulationParameters params)
{

  std::unique_ptr<Simulation::FlowMapTransitioner> transitioner = nullptr;
  const auto model = load_model_(params.user_params.model_name);
  MC::UniformLawINT law_param = {0, 0};

  auto simulation = init_simulation(
      exec_info, params, transitioner, model, std::move(law_param));
  return {std::move(simulation), params, std::move(transitioner), exec_info};
}

static void exec(const ExecInfo &exec_info, CaseData &&cased)
{
  if (cased.exec_info.current_rank == 0)
  {
    host_process(exec_info,
                 *cased.simulation,
                 cased.params,
                 std::move(cased.transitioner));
  }
  else
  {

    return workers_process(exec_info,
                           *cased.simulation,
                           cased.params,
                           std::move(cased.transitioner));
  }
}



static void
workers_process(const ExecInfo &exec,
                Simulation::SimulationUnit &simulation,
                SimulationParameters &params,
                std::unique_ptr<Simulation::FlowMapTransitioner> transitioner)
{

  double d_t = params.d_t;
  size_t n_compartments = simulation.mc_unit->domain.getNumberCompartments();
  MPI_Status status;

  MPI_W::IterationPayload payload(n_compartments * n_compartments,
                                  n_compartments);
  bool stop = false;
#pragma omp parallel
  while (!stop)
  {

    MPI_W::SIGNALS sign;
#pragma omp single
    {
      sign = MPI_W::try_recv<MPI_W::SIGNALS>(0, &status);
      if (sign == MPI_W::SIGNALS::STOP)
      {
        last_sync(exec, simulation);
        stop = true;
        
        // MPI_Abort(MPI_COMM_WORLD, 0);
      }
    }
    if(stop){break;}

#pragma omp single
    {
      payload.recv(0, &status);

      simulation.mc_unit->domain.setLiquidNeighbors(payload.neigbors);
      transitioner->update_flow(
          simulation, payload.liquid_flows, n_compartments);
      transitioner->advance(simulation);

      simulation.setVolumes(payload.gas_volumes, payload.liquid_volumes);
    }

    simulation.cycleProcess(d_t);

#pragma omp single
    {
      sync_step(exec, simulation);
      sync_prepare_next(exec, simulation);
    }
  }
}

static void restore_stdout(int original_stdout_fd,
                           std::streambuf *&original_buffer)
{
  if (is_stdout_redirect)
  {
    std::cout.rdbuf(original_buffer);
    fflush(stdout); // Flush the buffer to ensure all previous output is written
    dup2(original_stdout_fd,
         fileno(stdout)); // Restore the original file descriptor of stdout

    is_stdout_redirect = false;
  }
}

static int redirect_stdout(std::streambuf *&original_buffer,
                           std::stringstream &variable_stream)
{
  // Check if redirection is already active
  if (!is_stdout_redirect)
  {
    // Save the original buffer of std::cout
    original_buffer = std::cout.rdbuf();
    // Redirect std::cout to the stringstream
    std::cout.rdbuf(variable_stream.rdbuf());

    fflush(stdout); // Flush the buffer to ensure all previous output is written
    int original_stdout_fd =
        dup(fileno(stdout)); // Save the original file descriptor of stdout
    freopen("/dev/null", "w", stdout); // Redirect stdout to /dev/null

    // Set the flag to indicate that redirection is active
    is_stdout_redirect = true;
    return original_stdout_fd;
  }
}

