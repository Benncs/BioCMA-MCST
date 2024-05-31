
#include <cma_read/reactorstate.hpp>
#include <common/common.hpp>
#include <cstddef>
#include <host_specific.hpp>
#include <mpi_w/wrap_mpi.hpp>
#include <simulation/simulation.hpp>
#include <simulation/update_flows.hpp>
#include <sync.hpp>

// constexpr size_t n_particle_trigger_parralel = 1e6;

#ifdef DEBUG
#  define DEBUG_INSTRUCTION                                                    \

#else
#  define DEBUG_INSTRUCTION
#endif

#define FILL_PAYLOAD                                                           \
  mpi_payload.liquid_flows = reactor_state->liquid_flow.getViewFlows().data(); \
  mpi_payload.liquid_volumes = reactor_state->liquidVolume;                    \
  mpi_payload.gas_volumes = reactor_state->gasVolume;                          \
  mpi_payload.neigbors =                                                       \
      reactor_state->liquid_flow.getViewNeighors().to_const();

#define MPI_DISPATCH_MAIN                                                      \
  for (size_t __macro_j = 1; __macro_j < exec.n_rank; ++__macro_j)             \
  {                                                                            \
    MPI_W::send(MPI_W::SIGNALS::RUN, __macro_j);                               \
    mpi_payload.send(__macro_j);                                               \
  }


void main_loop(const SimulationParameters &params,
               const ExecInfo &exec,
               Simulation::SimulationUnit &simulation,
               std::shared_ptr<FlowIterator> _flow_handle,
               DataExporter *exporter)
{

  size_t iteration_count = 0;
  const size_t n_loop = params.n_different_maps;
  const double d_t = params.d_t;
  const size_t n_per_flowmap = _flow_handle->n_p_element();

  auto cache_liquid_flows = Simulation::BasicCacheHydro(n_loop);
  auto cache_gas_flows = Simulation::BasicCacheHydro(n_loop);

  MPI_W::HostIterationPayload mpi_payload;

  auto flow_map_iterator = _flow_handle->begin();

  const size_t n_iter_simulation = _flow_handle->totalSteps();

  auto *reactor_state = &(*flow_map_iterator);

  size_t dump_counter = 0;

  const size_t dump_number = std::min(
      n_iter_simulation, static_cast<size_t>(exporter->n_iter))-1;

  const size_t dump_interval = (n_iter_simulation) / (dump_number)+1;

  double current_time = 0;

  Simulation::update_flow(iteration_count,
                          n_per_flowmap,
                          n_loop,
                          simulation,
                          *reactor_state,
                          cache_liquid_flows,
                          cache_gas_flows,
                          params.is_two_phase_flow);

  simulation.setVolumes(reactor_state->gasVolume, reactor_state->liquidVolume);

#pragma omp parallel default(shared) 
  {

    // while (end_flag)
    for (size_t i = 0; i < n_iter_simulation; ++i)
    {

      DEBUG_INSTRUCTION

#pragma omp single
      {

        reactor_state = &(*flow_map_iterator);

        FILL_PAYLOAD;

        MPI_DISPATCH_MAIN;

        Simulation::update_flow(iteration_count,
                                n_per_flowmap,
                                n_loop,
                                simulation,
                                *reactor_state,
                                cache_liquid_flows,
                                cache_gas_flows,
                                params.is_two_phase_flow);

        simulation.setVolumes(reactor_state->gasVolume,
                              reactor_state->liquidVolume);
      }

      simulation.cycleProcess(d_t);

#pragma omp master
      {
        dump_counter++;
        if (dump_counter == dump_interval )
        {
          exporter->append(current_time,simulation.getCliqData(),
                           simulation.mc_unit->domain.getDistribution());                   
          dump_counter = 0;
        }
      }

#pragma omp single
      {

        sync_step(exec, simulation);
        simulation.step(d_t, *reactor_state);
        sync_prepare_next(exec, simulation);
        ++flow_map_iterator;
        current_time+=d_t;
        // end_flag = iterator != iterator_end;
      }
    }
  }
}