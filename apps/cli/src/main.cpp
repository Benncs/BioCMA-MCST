#include "messages/message_t.hpp"
#include "simulation/models/simple_model.hpp"
#include "simulation/models/types.hpp"
#include <any>
#include <common/common.hpp>

#include <host_specific.hpp>
#include <memory>
#include <siminit.hpp>

#include <flow_iterator.hpp>
#include <messages/wrap_mpi.hpp>
#include <mpi.h>
#include <reactorstate.hpp>
#include <simulation/simulation.hpp>
#include <simulation/transport.hpp>

#include <cstddef>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <sync.hpp>

#ifdef BIO_DYNAMIC_MODULE
#  include <import_py.hpp>
#endif

SimulationParameters moc_cli(int argc, char **argv);

void host_process(ExecInfo &exec,
                  Simulation::SimulationUnit &simulation,
                  SimulationParameters &params,
                  std::shared_ptr<FlowIterator> _flow_handle);

void workers_process(ExecInfo &exec,
                     Simulation::SimulationUnit &simulation,
                     SimulationParameters &params);

KModel load_model()
{
#ifdef BIO_DYNAMIC_MODULE
  return get_python_module("modules.test");
#else
  return simple_model;
#endif
}

int main(int argc, char **argv)
{

#ifdef BIO_DYNAMIC_MODULE
  auto _module_handle = init_dynamic_module();
#endif

  ExecInfo exec_info = MPI_W::init_mpi(argc, argv);

  SimulationParameters params = moc_cli(argc, argv);

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

  return 0;
}

void show(Simulation::SimulationUnit &simulation)
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

    std::cout << simulation.getCgas().row(1) << std::endl;
  }
  catch (...)
  {
  }
}

void host_process(ExecInfo &exec,
                  Simulation::SimulationUnit &simulation,
                  SimulationParameters &params,
                  std::shared_ptr<FlowIterator> _flow_handle)
{

  show(simulation);

  main_loop(params, exec, simulation, _flow_handle);

  show(simulation);

  host_dispatch(exec, MPI_W::SIGNALS::STOP);
}

void workers_process(ExecInfo &exec,
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

SimulationParameters moc_cli(int argc, char **argv)
{
  auto file = "/home/benjamin/Documenti/code/cpp/BIREM_generate/out/";
  return {10, 3, 1., {file}};
}