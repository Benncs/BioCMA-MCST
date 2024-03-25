#include "messages/impl_op.hpp"
#include "messages/message_t.hpp"
#include <Eigen/Dense>
#include <cstddef>
#include <messages/wrap_mpi.hpp>
#include <mpi.h>
#include <sync.hpp>
void sync_step(const ExecInfo &exec, Simulation::SimulationUnit &simulation)
{
  MPI_W::barrier();

  size_t nr = simulation.getCliq().rows();
  size_t nc = simulation.getCliq().cols();

  Eigen::MatrixXd &local_contribution = simulation.get_contribution();

  std::vector<double> total_contrib_data = MPI_W::gather<double>(
      std::span(local_contribution.data(), local_contribution.size()),
      exec.n_rank);

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
void sync_prepare_next(const ExecInfo &exec,
                       Simulation::SimulationUnit &simulation)
{
  MPI_W::barrier();
  simulation.clear_contribution();

  auto &data = simulation.getCliq();

  //We can use span here because we broadcast without changing size
  MPI_W::broadcast_span(std::span<double>(data.data(), data.size()), 0);
}