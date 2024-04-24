#include "messages/iteration_payload.hpp"
#include "messages/message_t.hpp"
#include <cli_parser.hpp>
#include <fstream>
#include <ios>
#include <simulation/update_flows.hpp>

#include <models/models.hpp>
#include <simulation/simulation.hpp>

#include <common/common.hpp>

#include <cma_read/flow_iterator.hpp>
#include <cma_read/reactorstate.hpp>

#include <host_specific.hpp>
#include <messages/wrap_mpi.hpp>
#include <post_process.hpp>
#include <rt_init.hpp>
#include <siminit.hpp>
#include <sync.hpp>

#ifdef USE_PYTHON_MODULE
#  include <pymodule/import_py.hpp>
#endif

#include <cstddef>
#include <exception>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>

#define SEND_MPI_SIG_STOP host_dispatch(exec, MPI_W::SIGNALS::STOP);

/*The following should be remove or moved */
constexpr bool verbose = true;                    // TODO REMOVE
constexpr bool redirect = false;                  // TODO REMOVE
static std::streambuf *original_buffer = nullptr; // TODO FIXME
static bool is_stdout_redirect = false;
static int original_stdout_fd = -1;
static void redirect_stdout(std::stringstream &variable_stream);
static void restore_stdout();
static void register_run(ExecInfo &exec, SimulationParameters &params);

static void workers_process(ExecInfo &exec,
                            Simulation::SimulationUnit &simulation,
                            SimulationParameters &params);

static void host_process(ExecInfo &exec,
                         Simulation::SimulationUnit &simulation,
                         SimulationParameters &params,
                         std::shared_ptr<FlowIterator> _flow_handle);

static KModel load_model();

static void exec(int argc, char **argv, SimulationParameters params);

int main(int argc, char **argv)
{
  init_environment();
  auto params_opt = parse_cli(argc, argv);
  if (!params_opt.has_value())
  {
    showHelp(std::cout);
    return -1;
  }

  try
  {

#ifdef USE_PYTHON_MODULE
    auto _interpreter_handle = init_python_interpreter();
#endif
    // TODO REMOVE
    if constexpr (redirect)
    {
      std::stringstream output_variable;
      output_variable.str("");
      redirect_stdout(output_variable);
      exec(argc, argv, std::move(params_opt.value()));
      restore_stdout();
      if constexpr (verbose)
      {
        std::cout << output_variable.str() << std::endl;
      }
    }
    else
    {
      exec(argc, argv, std::move(params_opt.value()));
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

static void exec(int argc, char **argv, SimulationParameters params)
{

  ExecInfo exec_info = runtime_init(argc, argv, params);

  std::shared_ptr<FlowIterator> _fd = nullptr;

  auto simulation = init_simulation(exec_info, params, _fd, load_model());

  if (exec_info.current_rank == 0)
  {

    if (_fd == nullptr)
    {
      throw std::runtime_error("Flow map are not loaded");
    }
    register_run(exec_info, params);
    host_process(exec_info, simulation, params, _fd);
  }
  else
  {

    workers_process(exec_info, simulation, params);
  }
}

static void show(Simulation::SimulationUnit &simulation)
{

  // auto distribution = simulation.mc_unit->domain.getDistribution();

  // for (auto &&i : distribution)
  // {
  //   std::cout << i << ", ";
  // }
  // std::cout << '\n';
  // try
  // {
  //   std::vector<double> totmas(simulation.unit->domain.n_compartments(), 0.);
  //   double cs = 0;
  //   for (auto &&p : simulation.container->to_process)
  //   {
  //     auto model = std::any_cast<std::shared_ptr<SimpleModel> &>(p.data);

  //     totmas[p.current_container] += p.weight * model->xi->mass;
  //     cs += p.weight * model->xi->mass;
  //   }
  //   std::cout << "mass: " << cs << std::endl;
  // }
  // catch (...)
  // {
  //   std::cout << std::endl;
  // }
  // std::cout << simulation.getCgas().row(1) << std::endl;
  // std::cout << "----Liquid---" << std::endl;
  // std::cout << simulation.getCliq().row(0) << std::endl;
}

static void host_process(ExecInfo &exec,
                         Simulation::SimulationUnit &simulation,
                         SimulationParameters &params,
                         std::shared_ptr<FlowIterator> _flow_handle)
{

  show(simulation);

  main_loop(params, exec, simulation, std::move(_flow_handle));

  SEND_MPI_SIG_STOP;
  last_sync(exec, simulation);

  show(simulation);

  post_process(params, simulation);
}

static void workers_process(ExecInfo &exec,
                            Simulation::SimulationUnit &simulation,
                            SimulationParameters &params)
{

  double d_t = params.d_t;
  size_t n_compartments = simulation.mc_unit->domain.n_compartments();
  MPI_Status status;

  size_t iteration_count = 0;
  const size_t n_loop = params.n_different_maps;
  auto liquid_flows = Simulation::BasicCacheMatflows(n_loop);
  std::vector<std::vector<size_t>> liquid_neighbors(n_compartments);

  MPI_W::IterationPayload payload(n_compartments * n_compartments,
                                  n_compartments);

  while (true)
  {

    auto sign = MPI_W::try_recv<MPI_W::SIGNALS>(0, &status);

    if (sign == MPI_W::SIGNALS::STOP)
    {
      last_sync(exec, simulation);
      return;
    }

    payload.recv(0, &status);

    // Neighbors could change during iteration so we have to allocate new
    // neighbors each time
    for (auto &&neighbors : liquid_neighbors)
    {
      neighbors = MPI_W::try_recv_v<size_t>(0, &status);
    }

    Simulation::update_flow(iteration_count,
                            n_loop,
                            simulation,
                            payload.liquid_flows,
                            n_compartments,
                            liquid_flows);

    simulation.mc_unit->domain.setLiquidNeighbors(liquid_neighbors);
    simulation.setVolumes(payload.gas_volumes, payload.liquid_volumes);

    simulation.cycleProcess(d_t);
    sync_step(exec, simulation);
    sync_prepare_next(exec, simulation);
  }
}

static KModel load_model()
{
#ifdef USE_PYTHON_MODULE
  std::cout << "LOADING modules.simple_model" << std::endl;
  return get_python_module("modules.simple_model_opt");
#else
  return get_simple_model();
#endif
}

static void register_run(ExecInfo &exec, SimulationParameters &params)
{
  // Open the file in append mode
  std::ofstream env(env_file_path(), std::ios_base::app);
  if (env.is_open())
  {
    append_date_time(env);
    env << exec;
    env << params;
    env << std::endl;
  }
  else
  {
    std::cerr << "Error: Unable to open file for writing\n";
  }
}

static void restore_stdout()
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

static void redirect_stdout(std::stringstream &variable_stream)
{
  // Check if redirection is already active
  if (!is_stdout_redirect)
  {
    // Save the original buffer of std::cout
    original_buffer = std::cout.rdbuf();
    // Redirect std::cout to the stringstream
    std::cout.rdbuf(variable_stream.rdbuf());

    fflush(stdout); // Flush the buffer to ensure all previous output is written
    original_stdout_fd =
        dup(fileno(stdout)); // Save the original file descriptor of stdout
    freopen("/dev/null", "w", stdout); // Redirect stdout to /dev/null

    // Set the flag to indicate that redirection is active
    is_stdout_redirect = true;
  }
}