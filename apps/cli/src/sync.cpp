#include "mc/events.hpp"
#include <messages/wrap_mpi.hpp>
#include <mpi.h>
#include <sync.hpp>

void sync_step(const ExecInfo &exec, Simulation::SimulationUnit &simulation)
{
  MPI_W::barrier();

  // Just use pointer to data wraped into span
  //  TODO: As we just gather we could use const data
  auto local_contribution = simulation.getContributionData();

  // Move span is useless but keep im mind gather idea of "moving" data
  std::vector<double> total_contrib_data =
      MPI_W::gather<double>(local_contribution, exec.n_rank);

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

void last_sync(const ExecInfo &exec, Simulation::SimulationUnit &simulation)
{
  MPI_W::barrier();
  auto tot_events =
      MC::EventContainer::reduce_local(simulation.mc_unit->ts_events);

  std::vector<size_t> total_contrib_data =
      MPI_W::gather<size_t>(tot_events.events, exec.n_rank);
  if (exec.current_rank == 0)
  {
    tot_events = MC::EventContainer::reduce(total_contrib_data);
  }

  simulation.mc_unit->ts_events = {
      tot_events}; // FIX IT because we will reduce twice (here + post process)
}