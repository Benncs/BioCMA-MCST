#include "data_exporter.hpp"
#include "models/light_model.hpp"
#include "models/monod.hpp"
#include "models/simple_model.hpp"
#include "mpi_w/impl_op.hpp"
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
constexpr bool verbose = true;   // TODO REMOVE
constexpr bool redirect = false; // TODO REMOVE
static bool is_stdout_redirect = false;

static int redirect_stdout(std::streambuf *&original_buffer,
                           std::stringstream &variable_stream);
static void restore_stdout(int original_stdout,
                           std::streambuf *&original_buffer);

static void
workers_process(ExecInfo &exec,
                Simulation::SimulationUnit &simulation,
                SimulationParameters &params,
                std::unique_ptr<Simulation::FlowMapTransitioner> transitioner);

static KModel load_model_(size_t index_model);
static void exec(int argc, char **argv, SimulationParameters params);

// static ExportParameters export_i = {10,"result_.h5"};

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
      std::streambuf *original_buffer = nullptr;
      auto original_stdout_fd =
          redirect_stdout(original_buffer, output_variable);

      exec(argc, argv, std::move(params_opt.value()));

      restore_stdout(original_stdout_fd, original_buffer);
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

  std::unique_ptr<Simulation::FlowMapTransitioner> transitioner = nullptr;
  constexpr size_t i_model = 1;
  const auto model = load_model_(i_model);
  MC::UniformLawINT law_param = {0, 1};

  auto simulation = init_simulation(
      exec_info, params, transitioner, model, std::move(law_param));

  if (exec_info.current_rank == 0)
  {
    host_process(exec_info, simulation, params, std::move(transitioner));
  }
  else
  {

    workers_process(exec_info, simulation, params, std::move(transitioner));
  }
}

void host_process(
    ExecInfo &exec,
    Simulation::SimulationUnit &simulation,
    SimulationParameters &params,
    std::unique_ptr<Simulation::FlowMapTransitioner> &&transitioner)
{
  // TODO: clean
  std::string fn = sappend_date_time("result_") + std::string(".h5");
  ExportParameters export_i = {static_cast<size_t>(params.final_time) * 3, fn};
  std::string name = "./results/" + export_i.filename;

  auto d = simulation.mc_unit->domain.getDistribution();
  DataExporter de(exec, params, name, simulation.getDim(), export_i.n_save, d);

  show(simulation);

  main_loop(params, exec, simulation, std::move(transitioner), &de);

  show(simulation);

  SEND_MPI_SIG_STOP;
  last_sync(exec, simulation);

  post_process(exec, params, simulation, &de);
}

static void
workers_process(ExecInfo &exec,
                Simulation::SimulationUnit &simulation,
                SimulationParameters &params,
                std::unique_ptr<Simulation::FlowMapTransitioner> transitioner)
{

  double d_t = params.d_t;
  size_t n_compartments = simulation.mc_unit->domain.getNumberCompartments();
  MPI_Status status;

  MPI_W::IterationPayload payload(n_compartments * n_compartments,
                                  n_compartments);
  while (true)
  {

    auto sign = MPI_W::try_recv<MPI_W::SIGNALS>(0, &status);

    if (sign == MPI_W::SIGNALS::STOP)
    {
      last_sync(exec, simulation);
      break;
    }

    payload.recv(0, &status);

    simulation.mc_unit->domain.setLiquidNeighbors(payload.neigbors);
    transitioner->update_flow(simulation, payload.liquid_flows, n_compartments);
    transitioner->advance(simulation);

    simulation.setVolumes(payload.gas_volumes, payload.liquid_volumes);

    simulation.cycleProcess(d_t);
    sync_step(exec, simulation);
    sync_prepare_next(exec, simulation);
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

static KModel load_model_(size_t index_model)
{

#ifdef USE_PYTHON_MODULE
  constexpr bool use_python_module = USE_PYTHON_MODULE;
  if constexpr (use_python_module)
  {
    std::cout << "LOADING modules.simple_model" << std::endl;
    return get_python_module("modules.simple_model_opt");
  }
#endif

  if (index_model == 0)
  {
    return get_simple_model();
  }
  if (index_model == 1)
  {
    return get_light_model();
  }
  if (index_model == 2)
  {
    return get_mond_model();
  }
  throw std::invalid_argument("bad model");
}