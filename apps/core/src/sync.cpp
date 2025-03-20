#include <biocma_cst_config.hpp>
#include <cstddef>
#include <mc/domain.hpp>
#include <mc/events.hpp>
#include <simulation/simulation.hpp>
#include <sync.hpp>

#ifndef NO_MPI
#  include <mpi_w/wrap_mpi.hpp>
#else
using MPI_Request = int;
// MOCK
// NOLINTBEGIN
namespace WrapMPI
{
  void barrier() {};
  template <typename T> std::vector<T> gather(auto, auto, int = 0) {};

  template <typename T> void gather_span(std::span<T>, std::span<const T>, size_t = 0){};

  void broadcast_span(auto data, auto n) {};
  template <typename T>
  int send_v(std::span<const T> data, size_t dest, size_t tag, bool send_size) noexcept {};

  namespace Async
  {

    template <typename T>
    int recv_span(MPI_Request& request, std::span<T> buf, size_t src, size_t tag) noexcept;
    void wait(MPI_Request& request);

  } // namespace Async
} // namespace WrapMPI
// NOLINTEND
#endif

void sync_step(const ExecInfo& exec, Simulation::SimulationUnit& simulation)
{
  PROFILE_SECTION("sync_step")
  // Sync is not needed in shared mode because only one unit is used in this
  // case, furthermore, barrier is already handled by kokkos at the end of
  // cycling process

  if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
  {
    // With multiple rank, we ensure that all finished kernel processing
    // then, rank 0 retrieves the local particle contribution in other
    // ranks
    // Reduce impl too long
    // #ifndef NO_MPI
    // if (exec.current_rank != 0)
    // {
    //   const auto local_contribution = simulation.getContributionData();
    //       std::vector<double> dummy_buffer(local_contribution.size(), 0.0); // Dummy buffer for
    //       non-root ranks
    //   MPI_Reduce(local_contribution.data(),
    //              dummy_buffer.data(),
    //              local_contribution.size(),
    //              MPI_DOUBLE,
    //              MPI_SUM,
    //              0,
    //              MPI_COMM_WORLD);
    // }
    // else
    // {
    //   auto host_contribution = simulation.getContributionData_mut();
    //   MPI_Reduce(MPI_IN_PLACE,host_contribution.data(),
    //              host_contribution.size(),
    //              MPI_DOUBLE,
    //              MPI_SUM,
    //              0,
    //              MPI_COMM_WORLD);
    // }
    // #endif

    // Old impl (GPU/CPU same node 432comps): - sync_step(REGION)   15.972034 1001
    // 0.015956 66.261864 36.022619 (6CPU ): - sync_step(REGION)   0.081290 1001 0.000081 0.918069
    // 0.712938

    const auto local_contribution = simulation.getContributionData();
    static std::vector<double> total_contrib_data;
    if (total_contrib_data.empty())
    {
      total_contrib_data.resize(local_contribution.size() * exec.n_rank, 0);
    }
    {
      PROFILE_SECTION("host:sync_step gathering")
      WrapMPI::gather_span<double>(total_contrib_data, local_contribution);
    }
    
    if (exec.current_rank == 0)
    {
      simulation.reduceContribs(total_contrib_data, exec.n_rank);
    }

    // Newimpl (GPU/CPU same node 432comps): - sync_step (REGION)   0.010833 1001 0.000011 0.079145
    // 0.025588
    //  (6CPU ): - sync_step (REGION)   0.109201 1001 0.000109 1.120788 0.845310
    

    // if (exec.current_rank != 0)
    // {
    //   const auto local_contribution = simulation.getContributionData();
    //   PROFILE_SECTION("worker:sync_step")
    //   WrapMPI::send_v(local_contribution, 0, 0, false);
    // }
    // else
    // {
    //   const auto local_contribution = simulation.getContributionData();
    //   static std::vector<MPI_Request> requests(exec.n_rank - 1);
    //   static std::vector<double> total_contrib_data;
    //   auto buffer_size = local_contribution.size();
    //   if (total_contrib_data.empty())
    //   {
    //     total_contrib_data.resize(buffer_size * (exec.n_rank), 0);
    //   }
    //   simulation.reduceContribs_per_rank(local_contribution);
    //   {
    //     PROFILE_SECTION("host:sync_step communication")
    //     for (int i = 1; i < exec.n_rank; ++i)
    //     {
    //       WrapMPI::Async::recv_span(
    //           requests[i - 1],
    //           std::span<double>(total_contrib_data.data() + i * buffer_size, buffer_size),
    //           i,
    //           0);
    //     }

    //     for (int i = 0; i < exec.n_rank - 1; ++i)
    //     {
    //       WrapMPI::Async::wait(requests[i]);
    //     }
    //   }
    //   simulation.reduceContribs(total_contrib_data, exec.n_rank);
    // }
  }
}

void sync_prepare_next(Simulation::SimulationUnit& simulation)
{
  PROFILE_SECTION("sync_prepare_next")
  simulation.clearContribution();
  if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
  {
    // In multiple rank context, we also need to broadcast the updated liquid
    // concentration computed by the host during the current step

    auto data = simulation.getCliqData(); // Get concentration ptr wrapped into span

    WrapMPI::barrier();
    // We can use span here because we broadcast without changing size
    WrapMPI::broadcast_span(data, 0);
  }
}

void last_sync(const ExecInfo& exec, Simulation::SimulationUnit& simulation)
{
  // For the last synchronisation, the aim for the host rank is to retrive all
  // local information of worker ranks
  if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
  {
    WrapMPI::barrier();

    // We first deal with event gathering, events can be broadcast and gather
    // easily as we deal with fixed sized array of integer

    auto& local_events = simulation.mc_unit->events;

    // We however have to declare events raw data as
    // vector to fit with WrapMPI API
    std::vector<size_t> total_events_data =
        WrapMPI::gather<size_t>(local_events.get_span(), exec.n_rank);

    if (exec.current_rank == 0)
    {
      // We could update 'local_events' for all ranks but it's useless as
      // simulation is finished and programm is going to exit
      local_events = MC::EventContainer::reduce(total_events_data);

      // TODO: MERGE PARTICLELIST
    }
    simulation.mc_unit->events = local_events; // FIX IT because we will reduce
                                               // twice (here + post process)
  }
}