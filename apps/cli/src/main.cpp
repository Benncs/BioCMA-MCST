#include "simulation/models/simple_model.hpp"
#include <any>
#include <common/common.hpp>

#include <host_specific.hpp>
#include <memory>
#include <pinit.hpp>

#include <flow_iterator.hpp>
#include <messages/pmessage.hpp>
#include <mpi.h>
#include <reactorstate.hpp>
#include <simulation/simulation.hpp>
#include <simulation/transport.hpp>

#include <cstddef>
#include <iostream>
#include <messages/init.hpp>
#include <stdexcept>
#include <thread>

SimulationParameters moc_cli(int argc, char **argv);

void host_process(ExecInfo &exec,
                  Simulation::SimulationUnit &simulation,
                  SimulationParameters &params,
                  FlowIterator *_flow_handle);

void workers_process(ExecInfo &exec,
                     Simulation::SimulationUnit &simulation,
                     SimulationParameters &params);

void sync_step(ExecInfo &exec, Simulation::SimulationUnit &simulation);

void sync_step_2(ExecInfo &exec, Simulation::SimulationUnit &simulation);

int main(int argc, char **argv)
{

  ExecInfo exec_info = init_mpi(argc, argv);
  std::cout << "NUM thread per process " << exec_info.thread_per_process
            << std::endl;
  SimulationParameters params = moc_cli(argc, argv);

  FlowIterator *_fd = nullptr;
  auto simulation = pinit(exec_info, params, &_fd);

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

  if (exec_info.current_rank == 0)
  {
    delete _fd;
  }

  MPI_Finalize();

  return 0;
}

void host_process(ExecInfo &exec,
                  Simulation::SimulationUnit &simulation,
                  SimulationParameters &params,
                  FlowIterator *_flow_handle)
{
  ReactorState *f = _flow_handle->get();

  double d_t = params.d_t;
  // std::cout << simulation.getC() << std::endl;

  std::vector<double> totmas(simulation.unit->domain.n_compartments(),0.);

  double cs = 0;
  for(auto&& p : simulation.container->to_process){
    auto model = std::any_cast<std::shared_ptr<SimpleModel>&>(p.data);
    totmas[p.current_container]+= p.weight*model->xi->mass;
    cs+=p.weight*model->xi->mass;
  }
  std::cout<<"mass: "<<cs<<std::endl;
  for(auto&& i:totmas)
  {
    std::cout<<i<<"\t";
  }
  std::cout<<std::endl;

  while (f != nullptr)
  {
    
    auto &liq_flow = f->liquid_flow;
    auto &gas_flow = f->gas_flow;
    std::thread _h_liquid(computeLiquidFlow,
                          std::cref(exec),
                          std::ref(simulation),
                          std::ref(liq_flow));
    std::thread _h_gas(
        computeGasFlow, std::ref(simulation), std::ref(gas_flow));

    _h_liquid.join();
    _h_gas.join();



    simulation.cycle_process(d_t);

    sync_step(exec, simulation);
    simulation.step(d_t);
    sync_step_2(exec, simulation);

    _flow_handle->next(); // this could be done async ?
    f = _flow_handle->get();
  }
  // std::cout << simulation.getC() << std::endl;


  std::vector<double> totmas2(simulation.unit->domain.n_compartments(),0.);
  cs = 0;
  for(auto&& p : simulation.container->to_process){
    auto model = std::any_cast<std::shared_ptr<SimpleModel>&>(p.data);
    totmas2[p.current_container]+= p.weight*model->xi->mass;
    cs+=p.weight*model->xi->mass;;
  }

  std::cout<<"mass: "<<cs<<std::endl;
  for(auto&& i:totmas2)
  {
    std::cout<<i<<"\t";
  }
  std::cout<<std::endl;


  host_dispatch(exec, MPI_SIGNALS::STOP);
  // DO SOMETHING
}

void workers_process(ExecInfo &exec,
                     Simulation::SimulationUnit &simulation,
                     SimulationParameters &params)
{

  double d_t = params.d_t;
  size_t n_compartments = simulation.unit->domain.n_compartments();
  // DO SOMETHING
  MPI_Status status;
  while (true)
  {
    MPI_SIGNALS sign;
    MPI_Recv(
        &sign, sizeof(sign), MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    switch (sign)
    {
    case MPI_SIGNALS::STOP:
    {
      return;
    }
    case MPI_SIGNALS::RUN:
    {
      std::vector<double> flows;
      size_t sf;
      MPI_Recv(
          &sf, 1, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      flows.resize(sf);

      MPI_Recv(flows.data(),
               static_cast<int>(sf),
               MPI_DOUBLE,
               0,
               MPI_ANY_TAG,
               MPI_COMM_WORLD,
               &status);

      const auto mat_f = FlowmapToMat(flows, n_compartments);
      const auto mat_transition = Simulation::get_transition_matrix(mat_f);

      simulation.setLiquidFlow(
          Simulation::MatFlow(std::move(mat_f), std::move(mat_transition)));

      simulation.cycle_process(d_t);
      sync_step(exec, simulation);
      sync_step_2(exec, simulation);
      break;
    }
    }
  }
}

#include <Eigen/Dense>
void sync_step(ExecInfo &exec, Simulation::SimulationUnit &simulation)
{
  MPI_Barrier(MPI_COMM_WORLD);

  size_t nr = simulation.getC().rows();
  size_t nc = simulation.getC().cols();

  Eigen::MatrixXd &local_contribution = simulation.get_contribution();
  std::vector<double> total_contrib_data(local_contribution.size() *
                                         exec.n_rank);

  MPI_Gather(local_contribution.data(),
             static_cast<int>(local_contribution.size()),
             MPI_DOUBLE,
             total_contrib_data.data(),
             static_cast<int>(local_contribution.size()),
             MPI_DOUBLE,
             0,
             MPI_COMM_WORLD);

  if (exec.current_rank == 0)
  {
    Eigen::MatrixXd total_contrib = Eigen::MatrixXd(nr, nc);
    total_contrib.setZero();
    for (int i = 0; i < static_cast<int>(exec.n_rank); ++i)
    {
      total_contrib +=
          Eigen::Map<Eigen::MatrixXd>(&total_contrib_data[i * nr * nc],
                                      static_cast<int>(nr),
                                      static_cast<int>(nc));
    }

    simulation.get_contribution() = total_contrib;
  }
}
void sync_step_2(ExecInfo &exec, Simulation::SimulationUnit &simulation)
{
  MPI_Barrier(MPI_COMM_WORLD);
  simulation.clear_contribution();

  auto &data = simulation.getC();
  MPI_Bcast(data.data(),
            static_cast<int>(data.size()),
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD);
}

SimulationParameters moc_cli(int argc, char **argv)
{
  auto file = "/home/benjamin/Documenti/code/cpp/BIREM_generate/out/";
  return {1'000, 3, 1000., {file}};
}