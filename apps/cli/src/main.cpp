#include "common/execinfo.hpp"
#include <common/common.hpp>
#include <concepts>
#include <cstddef>
#include <flow_iterator.hpp>
#include <iostream>
#include <messages/pmessage.hpp>
#include <mpi.h>
#include <pinit.hpp>
#include <simulation/simulation.hpp>
#include <simulation/transport.hpp>
#include <stdexcept>

#include <messages/init.hpp>

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
  FlowInfo *f = _flow_handle->get();
  auto dis = simulation.unit->domain.getDistribution();
  for (auto &&i : dis)
  {
    std::cout << i << ", ";
  }
  std::cout << "----" << std::endl;

  // size_t n_t = _flow_handle->totalSteps();
  double d_t = params.d_t;
  std::cout << simulation.getC() << std::endl;

  while (f != nullptr)
  {

    const auto mat_f = FlowmapToMat(f->flows.data(), f->flows.getN());
    const auto mat_transition = Simulation::get_transition_matrix(mat_f);
    host_dispatch(exec, MPI_SIGNALS::RUN, f->flows.data());

    Simulation::MatFlow mf = {mat_f, mat_transition};
    simulation.setFLows(&mf);
    // TODO FIX UNIQUE PTR

    simulation.cycle_process(d_t);
    _flow_handle->next(); // this could be done async ?
    sync_step(exec, simulation);

    simulation.step(d_t);

    sync_step_2(exec, simulation);

    f = _flow_handle->get();
  }
  std::cout << simulation.getC() << std::endl;
  dis = simulation.unit->domain.getDistribution();
  size_t icc = 0;
  for (auto &&i : dis)
  {
    std::cout << i << ", ";
    icc += i;
  }

  std::cout << "---" << icc << std::endl;

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
               sf,
               MPI_DOUBLE,
               0,
               MPI_ANY_TAG,
               MPI_COMM_WORLD,
               &status);

      const auto mat_f = FlowmapToMat(flows, n_compartments);
      const auto mat_transition = Simulation::get_transition_matrix(mat_f);

      Simulation::MatFlow mf = {mat_f, mat_transition};
      simulation.setFLows(&mf);

      simulation.cycle_process(d_t);
      sync_step(exec, simulation);
      break;
    }
    }
  }
}

void sync_step(ExecInfo &exec, Simulation::SimulationUnit &simulation)
{
  MPI_Barrier(MPI_COMM_WORLD);

  // std::vector<int> full_rng;
  // if (exec.current_rank == 0)
  // {
  //   full_rng.resize(
  //       exec.n_rank); // Assuming exec.n_rank is the number of MPI processes
  // }
  // // Gather the 'rng' integer from all processes into 'full_rng' on process 0
  // MPI_Gather(&simulation.rng,
  //            1,
  //            MPI_INT,
  //            full_rng.data(),
  //            1,
  //            MPI_INT,
  //            0,
  //            MPI_COMM_WORLD);

  // // Broadcast 'full_rng' from rank 0 to all processes
  // MPI_Bcast(full_rng.data(),
  //           static_cast<int>(full_rng.size()),
  //           MPI_INT,
  //           0,
  //           MPI_COMM_WORLD);
}

void sync_step_2(ExecInfo &exec, Simulation::SimulationUnit &simulation)
{
}

SimulationParameters moc_cli(int argc, char **argv)
{
  auto file = "/home/benjamin/Documenti/code/cpp/BIREM_Project/Example/"
              "Sanofi/CMA_export/35rpm_flowL.raw";
  return {500'000, 4, 10., {file}};
}