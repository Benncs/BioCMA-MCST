#include <biocma_cst_config.hpp>
#include <cstddef>
#include <mc/domain.hpp>
#include <mc/events.hpp>
#include <sync.hpp>
#include <simulation/simulation.hpp>

#ifndef NO_MPI
#  include <mpi_w/wrap_mpi.hpp>
#else
// MOCK
// NOLINTBEGIN
namespace WrapMPI
{
  void barrier() {};
  template <typename T>
  std::vector<T> gather(auto s, auto n, int i = 0) {};

  void broadcast_span(auto data, auto n) {};

} // namespace WrapMPI
// NOLINTEND
#endif

void sync_step(const ExecInfo &exec, Simulation::SimulationUnit &simulation)
{
  PROFILE_SECTION("sync_step")
  // Sync is not needed in shared mode because only one unit is used in this
  // case, furthermore, barrier is already handled by kokkos at the end of
  // cycling process
  
  if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
  {
    // With multiple rank, we ensure than all finished kernel processing
    // afterwards, we rank 0 retrieves the local particle contribution in other
    // ranks

    WrapMPI::barrier();

    // Just use pointer to data wraped into span
    //  TODO: As we just gather we could use const data
    auto local_contribution = simulation.getContributionData();

    std::vector<double> total_contrib_data =
        WrapMPI::gather<double>(local_contribution, exec.n_rank);

    if (exec.current_rank == 0)
    {
      simulation.reduceContribs(total_contrib_data, exec.n_rank);
    }
  }
}

void sync_prepare_next(Simulation::SimulationUnit &simulation)
{
  PROFILE_SECTION("sync_prepare_next")
  simulation.clearContribution();
  if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
  {
    // In multiple rank context, we also need to broadcast the updated liquid
    // concentration computed by the host during the current step
    WrapMPI::barrier();

    auto data =
        simulation.getCliqData(); // Get concentration ptr wrapped into span

    // We can use span here because we broadcast without changing size
    WrapMPI::broadcast_span(data, 0);
  }
}

void last_sync(const ExecInfo &exec, Simulation::SimulationUnit &simulation)
{
  // For the last synchronisation, the aim for the host rank is to retrive all
  // local information of worker ranks
  if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
  {
    WrapMPI::barrier();

    // We first deal with event gathering, events can be broadcast and gather
    // easily as we deal with fixed sized array of integer

    auto &local_events = simulation.mc_unit->events;

    // We however have to declare events raw data as
    // vector to fit with WrapMPI API
    std::vector<size_t> total_events_data =
        WrapMPI::gather<size_t>(local_events.get_span(), exec.n_rank);

    auto local_distribution = simulation.mc_unit->domain.getRepartition();

    const std::size_t local_distribution_size = local_distribution.size();

    // Functor to get the total number of particle in all containers
    // const auto visitor_sync = [&exec](auto &&container) {
    //   return WrapMPI::gather_reduce<size_t>(container.process_size(),
    //                                       exec.n_rank);
    // };
    // auto total_particle_number = WrapMPI::gather_reduce<size_t>(
    //     std::visit(visitor_sync, simulation.mc_unit->container),
    //     exec.n_rank);

    // Distribution can be easily merged because all ranks have the same number
    // of compartment so 'local_distribution_size' is the same for all ranks
    const auto merged_distribution =
        WrapMPI::gather<size_t>(local_distribution, exec.n_rank, 0);

    if (exec.current_rank == 0)
    {
      // We could update 'local_events' for all ranks but it's useless as
      // simulation is finished and programm is going to exit
      local_events = MC::EventContainer::reduce(total_events_data);

      // simulation.mc_unit->domain.in_place_reduce(merged_distribution, local_distribution_size, exec.n_rank);
      
      // //WARN: reduce doesn't copy 
      // auto reduced_domain = MC::ReactorDomain::reduce(
      //     merged_distribution, local_distribution_size, exec.n_rank);

      // simulation.mc_unit->domain =
      //     std::move(reduced_domain); // Set new domain (uppdate particle
      //                                // location only on container side)

      // TODO: MERGE PARTICLELIST
    }
    simulation.mc_unit->events = local_events; // FIX IT because we will reduce
                                               // twice (here + post process)
  }
}