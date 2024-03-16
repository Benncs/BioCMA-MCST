#include <messages/pmessage.hpp>
#include <simulation/simulation.hpp>
#include <common/common.hpp>
#include <flow_iterator.hpp>
#include <iostream>
#include <mpi.h>
#include <pinit.hpp>
#include <stdexcept>

ExecInfo init_mpi(int argc, char **argv);
SimulationParameters moc_cli(int argc, char **argv);

void master_dispatch(ExecInfo &info, MPI_SIGNALS &&sign)
{

  // Send message to all other processes
  for (size_t j = 1; j < info.n_rank; ++j)
  {
    MPI_Send(&sign, sizeof(sign), MPI_CHAR, j, 0, MPI_COMM_WORLD);
  }

  // // Broadcast each element in args to all processes
  // for (auto& arg : args) {
  //     MPI_Bcast(&arg, 1, MPI_DOUBLE, 0, _comm);
  // }
}

void master_process(ExecInfo &exec,
                    SimulationUnit &simulation,
                    SimulationParameters &params,
                    FlowIterator *_flow_handle);

void slave_process(ExecInfo &exec,
                   SimulationUnit &simulation,
                   SimulationParameters &params);

void sync_step(ExecInfo &exec, SimulationUnit &simulation);

int main(int argc, char **argv)
{

  auto exec_info = init_mpi(argc, argv);
  srand(time(NULL) + exec_info.current_rank);

  SimulationParameters params = moc_cli(argc, argv);

  FlowIterator *_fd = nullptr;
  auto simulation = pinit(exec_info, params, &_fd);

  test_common();

  if (exec_info.current_rank == 0)
  {
    if (_fd == nullptr)
    {
      throw std::runtime_error("Flow map are not loaded");
    }
    master_process(exec_info, simulation, params, _fd);
  }
  else
  {
    slave_process(exec_info, simulation, params);
  }

  if (exec_info.current_rank == 0)
  {
    delete _fd;
  }

  MPI_Finalize();

  return 0;
}

ExecInfo init_mpi(int argc, char **argv)
{
  ExecInfo info;

  int rank, size;

  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Get the rank of the current process
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Get the total number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  info.current_rank = static_cast<size_t>(rank);
  info.n_rank = static_cast<size_t>(size);

  return info;
}

void master_process(ExecInfo &exec,
                    SimulationUnit &simulation,
                    SimulationParameters &params,
                    FlowIterator *_flow_handle)
{
  FlowInfo *f = _flow_handle->get();

  while (f != nullptr)
  {
    // DO SOMETHING
    master_dispatch(exec, MPI_SIGNALS::RUN);
    simulation.cycle_process();
    _flow_handle->next(); // this could be done async ?
    sync_step(exec, simulation);
    f = _flow_handle->get();
  }
  master_dispatch(exec, MPI_SIGNALS::STOP);
  // DO SOMETHING
}

void slave_process(ExecInfo &exec,
                   SimulationUnit &simulation,
                   SimulationParameters &params)
{
  // DO SOMETHING
  MPI_Status status;
  while (true)
  {
    MPI_SIGNALS sign;
    MPI_Recv(
        &sign, sizeof(sign), MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    int sender_rank = status.MPI_SOURCE;

    switch (sign)
    {
    case MPI_SIGNALS::STOP:
    {
      return;
    }
    case MPI_SIGNALS::RUN:
    {
      simulation.cycle_process();
      sync_step(exec, simulation);
      break;
    }
    }
  }
}

void sync_step(ExecInfo &exec, SimulationUnit &simulation)
{
  MPI_Barrier(MPI_COMM_WORLD);

  std::vector<int> full_rng;
  if (exec.current_rank == 0)
  {
    full_rng.resize(
        exec.n_rank); // Assuming exec.n_rank is the number of MPI processes
  }
  // Gather the 'rng' integer from all processes into 'full_rng' on process 0
  MPI_Gather(&simulation.rng,
             1,
             MPI_INT,
             full_rng.data(),
             1,
             MPI_INT,
             0,
             MPI_COMM_WORLD);

  // Broadcast 'full_rng' from rank 0 to all processes
  MPI_Bcast(full_rng.data(), full_rng.size(), MPI_INT, 0, MPI_COMM_WORLD);

  // if (exec.current_rank == 0)
  // {
  //   for (auto i : full_rng)
  //   {
  //     std::cout << i << std::endl;
  //   }
  // }
}

SimulationParameters moc_cli(int argc, char **argv)
{
  auto file = "/home/benjamin/Documenti/code/cpp/BIREM_Project/Example/"
              "Sanofi/CMA_export/35rpm_flowL.raw";
  return {5'000'000, 10., {file}};
}