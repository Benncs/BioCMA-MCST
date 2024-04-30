
#include <cma_read/reactorstate.hpp>
#include <common/common.hpp>
#include <cstddef>
#include <host_specific.hpp>
#include <mpi_w/wrap_mpi.hpp>
#include <simulation/simulation.hpp>
#include <simulation/update_flows.hpp>
#include <sync.hpp>

constexpr size_t n_particle_trigger_parralel = 1e6;

#ifdef DEBUG
#  define DEBUG_INSTRUCTION                                                    \
    // static size_t __debug_loop_cntr = 0; \
    // std::cout << __debug_loop_cntr << '\n'; \
    // __debug_loop_cntr++;
#else
#  define DEBUG_INSTRUCTION
#endif

#define FILL_PAYLOAD                                                           \
  payload.liquid_flows = reactor_state.liquid_flow.flows.data();               \
  payload.liquid_volumes = reactor_state.liquidVolume;                         \
  payload.gas_volumes = reactor_state.gasVolume;

#define MPI_DISPATCH_MAIN                                                      \
  for (size_t j = 1; j < exec.n_rank; ++j)                                     \
  {                                                                            \
    MPI_W::send(MPI_W::SIGNALS::RUN, j);                                       \
    payload.send(j);                                                           \
  }                                                                            \
  for (const auto &neighbor : reactor_state.liquid_flow.neigbors)              \
  { /*Send each neighbor vector to all processes*/                             \
    for (size_t j = 1; j < exec.n_rank; ++j)                                   \
    {                                                                          \
      MPI_W::send_v<size_t>(neighbor, j, 0);                                   \
    }                                                                          \
  }

void main_loop(const SimulationParameters &params,
               const ExecInfo &exec,
               Simulation::SimulationUnit &simulation,
               std::shared_ptr<FlowIterator> _flow_handle)
{

  size_t iteration_count = 0;
  const size_t n_loop = params.n_different_maps;
  const double d_t = params.d_t;

  auto liquid_flows = Simulation::BasicCacheMatflows(n_loop);
  auto gas_flows = Simulation::BasicCacheMatflows(n_loop);

  MPI_W::HostIterationPayload payload;

  std::cout << params.final_time << " " << d_t << '\n';

  auto iterator = _flow_handle->begin();
  const auto end = _flow_handle->end();
  auto &reactor_state = *iterator;
  bool end_flag = iterator != end;

#pragma omp parallel default(shared) num_threads(exec.thread_per_process) 
  {

    while (end_flag)
    {

#pragma omp single
      {

        reactor_state = *iterator;

        FILL_PAYLOAD

        DEBUG_INSTRUCTION

        MPI_DISPATCH_MAIN

        /*
        For the two following function calls, pass non-owning data types:
        - 'volume' as a std::span
        - 'flows' as a reference

        It's important to ensure that the lifetime of the volumes
        (reactor_state.gasVolume and reactor_state.liquidVolume)
        and 'liquid_flows', 'gas_flows' extends at least until the 'simulation'
        object finishes using them.

        - For 'liquid_flows' and 'gas_flows':
          These are used within the loop iteration (lifetime < 'loop_iteration')
          and stored as cache in owning types outside of the loop (lifetime >
        'loop_iteration'). This setup ensures that there are no lifetime issues.

        For 'volumes':
          'reactor_state' has at least the lifetime of 'loop_iteration',
          and the volumes are used during "cycleProcess" within the same
        'loop_iteration' (lifetime < 'loop_iteration'). This usage pattern
        ensures that the volumes are valid and accessible during the entire
        process of 'cycleProcess' within the loop iteration.
    */
        Simulation::update_flow(iteration_count,
                                n_loop,
                                simulation,
                                reactor_state,
                                liquid_flows,
                                gas_flows);

        simulation.setVolumes(reactor_state.gasVolume,
                              reactor_state.liquidVolume);
      }

      simulation.cycleProcess(d_t);

#pragma omp single
      {
        sync_step(exec, simulation);
        simulation.step(d_t, reactor_state);
        sync_prepare_next(exec, simulation);
        ++iterator;
        end_flag = iterator != end;
      }
    }
  }
}