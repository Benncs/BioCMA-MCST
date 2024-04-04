#include "cli_parser.hpp"
#include <simulation/models/models.hpp>
#include <simulation/simulation.hpp>
#include <simulation/transport.hpp>

#include <common/common.hpp>

#include <flow_iterator.hpp>
#include <reactorstate.hpp>

#include <host_specific.hpp>
#include <messages/wrap_mpi.hpp>
#include <post_process.hpp>
#include <rt_init.hpp>
#include <siminit.hpp>
#include <sync.hpp>
#include <cli_parser.hpp>
#ifdef BIO_DYNAMIC_MODULE
#  include <import_py.hpp>
#endif

#include <cstddef>
#include <exception>
#include <iostream>
#include <memory>
#include <stdexcept>


static void workers_process(ExecInfo &exec,
                            Simulation::SimulationUnit &simulation,
                            SimulationParameters &params);

static void host_process(ExecInfo &exec,
                         Simulation::SimulationUnit &simulation,
                         SimulationParameters &params,
                         std::shared_ptr<FlowIterator> _flow_handle);

static KModel load_model();

static void exec(int argc, char **argv);

int main(int argc, char **argv)
{
  try
  {

#ifdef BIO_DYNAMIC_MODULE
    auto _module_handle = init_dynamic_module();
#endif
    init_environement();
    exec(argc, argv);
  }
  catch(std::invalid_argument &e)
  {
    std::cerr<<e.what()<<std::endl;
    showHelp(std::cout);

    return -1;
  }
#ifdef DEBUG
  catch (std::exception &e)

  {
    std::cerr << e.what() << std::endl;
    return -1;
  }
#endif
  catch (...)
  {
    std::cerr << "Internal error" << std::endl;
    return -1;
  }

  return 0;
}

static void exec(int argc, char **argv)
{
  SimulationParameters params = parseCLI(argc, argv);
  ExecInfo exec_info = runtime_init(argc, argv, params);

  std::shared_ptr<FlowIterator> _fd = nullptr;
  auto simulation = sim_init(exec_info, params, _fd, load_model());

  if (exec_info.current_rank == 0)
  {
    if (_fd == nullptr)
    {
      throw std::runtime_error("Flow map are not loaded");
    }

    host_process(exec_info, simulation, params, _fd);
  }
  else
  {
    workers_process(exec_info, simulation, params);
  }
}

static void show(Simulation::SimulationUnit &simulation)
{
  try
  {
    std::vector<double> totmas(simulation.unit->domain.n_compartments(), 0.);
    double cs = 0;
    for (auto &&p : simulation.container->to_process)
    {
      auto model = std::any_cast<std::shared_ptr<SimpleModel> &>(p.data);

      totmas[p.current_container] += p.weight * model->xi->mass;
      cs += p.weight * model->xi->mass;
    }
    std::cout << "mass: " << cs << std::endl;
  }
  catch (...)
  {
    std::cout << std::endl;
  }
  std::cout << simulation.getCgas().row(1) << std::endl;
  std::cout<<"----Liquid---"<<std::endl;
  std::cout << simulation.getCliq().row(0) << std::endl;
}

static void host_process(ExecInfo &exec,
                         Simulation::SimulationUnit &simulation,
                         SimulationParameters &params,
                         std::shared_ptr<FlowIterator> _flow_handle)
{

  show(simulation);

  main_loop(params, exec, simulation, _flow_handle);

  show(simulation);

  host_dispatch(exec, MPI_W::SIGNALS::STOP);

  post_process(simulation);
}

static void workers_process(ExecInfo &exec,
                            Simulation::SimulationUnit &simulation,
                            SimulationParameters &params)
{

  double d_t = params.d_t;
  size_t n_compartments = simulation.unit->domain.n_compartments();
  MPI_Status status;
  while (true)
  {

    auto sign = MPI_W::try_recv<MPI_W::SIGNALS>(0, &status);

    if (sign == MPI_W::SIGNALS::STOP)
    {
      return;
    }

    auto flows = MPI_W::try_recv_v<double>(0);

    const auto mat_f = FlowmapToMat(flows, n_compartments);
    const auto mat_transition = Simulation::get_transition_matrix(mat_f);

    simulation.setLiquidFlow(
        Simulation::MatFlow(std::move(mat_f), std::move(mat_transition)));

    simulation.cycle_process(d_t);
    sync_step(exec, simulation);
    sync_prepare_next(exec, simulation);
  }
}


static KModel load_model()
{
#ifdef BIO_DYNAMIC_MODULE
  return get_python_module("modules.simple_model");
#else
  return simple_model;
#endif
}
