

#include <cstddef>
#include <messages/wrap_mpi.hpp>
#include <mpi.h>
#include <sync.hpp>

void sync_step(const ExecInfo &exec, Simulation::SimulationUnit &simulation)
{
  MPI_W::barrier();

  // Just use pointer to data wraped into span
  //  TODO: As we just gather we could use const data
  auto local_contribution = simulation.get_contributionData();

  // Move span is useless but to keep gather idea of "moving" data
  std::vector<double> total_contrib_data =
      MPI_W::gather<double>(std::move(local_contribution), exec.n_rank);

  if (exec.current_rank == 0)
  {
    simulation.reduceContribs(total_contrib_data, exec.n_rank);
  }
}
void sync_prepare_next(const ExecInfo & /*exec*/,
                       Simulation::SimulationUnit &simulation)
{
  MPI_W::barrier();
  simulation.clearContribution();

  auto data =
      simulation.getCliqData(); // Get concentration ptr wrapped into span

  // We can use span here because we broadcast without changing size
  MPI_W::broadcast_span(data, 0);
}