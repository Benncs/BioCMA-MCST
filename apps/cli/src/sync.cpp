#include "mc/container_state.hpp"
#include "mc/domain.hpp"
#include "mc/events.hpp"
#include "mc/particles/mcparticles.hpp"
#include "rt_init.hpp"
#include <cstddef>
#include <mpi_w/wrap_mpi.hpp>
#include <sync.hpp>
#include <variant>

void sync_step(const ExecInfo &exec, Simulation::SimulationUnit &simulation)
{
  MPI_W::barrier();

  // Just use pointer to data wraped into span
  //  TODO: As we just gather we could use const data
  auto local_contribution = simulation.getContributionData();

  std::vector<double> total_contrib_data =
      MPI_W::gather<double>(local_contribution, exec.n_rank);

  if (exec.current_rank == 0)
  {
    simulation.reduceContribs(total_contrib_data, exec.n_rank);
  }
}
void sync_prepare_next(const ExecInfo &exec,
                       Simulation::SimulationUnit &simulation)
{
  if constexpr (RT::use_mpi)
  {

    MPI_W::barrier();
    simulation.clearContribution();

    auto data =
        simulation.getCliqData(); // Get concentration ptr wrapped into span

    // We can use span here because we broadcast without changing size
    MPI_W::broadcast_span(data, 0);
  }
  else
  {
    simulation.clearContribution();
  }
}

void last_sync(const ExecInfo &exec, Simulation::SimulationUnit &simulation)
{
  MPI_W::barrier();

  auto &tot_events = simulation.mc_unit->events;
  // MC::EventContainer::reduce_local(simulation.mc_unit->events);

  auto span_events = std::span<size_t>(tot_events._events.data(),MC::EventContainer::number_event_type);

  std::vector<size_t> total_contrib_data =
      MPI_W::gather<size_t>(span_events, exec.n_rank);

  auto local = simulation.mc_unit->domain.getDistribution();

  auto local_size = local.size();

  // auto local_particle = simulation.mc_unit->container.to_process.data();

  const auto visitor_sync = [&exec](auto &&container)
  {
    return MPI_W::gather_reduce<size_t>(container.process_size(), exec.n_rank);
  };

  auto total_particle = MPI_W::gather_reduce<size_t>(
      std::visit(visitor_sync, simulation.mc_unit->container), exec.n_rank);

  auto tot_distrib = MPI_W::gather<size_t>(local, exec.n_rank, 0);

  if (exec.current_rank == 0)
  {
    tot_events = MC::EventContainer::reduce(total_contrib_data);
    auto reduced =
        MC::ReactorDomain::reduce(tot_distrib, local_size, exec.n_rank);
    simulation.mc_unit->domain = std::move(reduced);
    std::cout << "nparticle " << total_particle << std::endl;
    // simulation.mc_unit->container.to_process.data() = total_particle;
  }
  simulation.mc_unit->events =
      tot_events; // FIX IT because we will reduce twice (here + post process)
}