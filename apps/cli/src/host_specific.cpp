
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
#  define DEBUG_INSTRUCTION

#else
#  define DEBUG_INSTRUCTION
#endif

#define FILL_PAYLOAD                                                           \
  mpi_payload.liquid_flows =                                                   \
      current_reactor_state->liquid_flow.getViewFlows().data();                \
  mpi_payload.liquid_volumes = current_reactor_state->liquidVolume;            \
  mpi_payload.gas_volumes = current_reactor_state->gasVolume;                  \
  mpi_payload.neigbors =                                                       \
      current_reactor_state->liquid_flow.getViewNeighors().to_const();

#define MPI_DISPATCH_MAIN                                                      \
  for (size_t __macro_j = 1; __macro_j < exec.n_rank; ++__macro_j)             \
  {                                                                            \
    MPI_W::send(MPI_W::SIGNALS::RUN, __macro_j);                               \
    mpi_payload.send(__macro_j);                                               \
  }


void main_loop(const SimulationParameters &params,
               const ExecInfo &exec,
               Simulation::SimulationUnit &simulation,
   std::unique_ptr<Simulation::FlowMapTransitioner> transitioner,
               DataExporter *exporter)
{

  const double d_t = params.d_t;

  MPI_W::HostIterationPayload mpi_payload;

  auto *current_reactor_state = &transitioner->get_unchecked_mut(0);

  ReactorState *next_reactor_state = nullptr;

  const size_t n_iter_simulation = transitioner->get_n_timestep();

  const size_t dump_number =
      std::min(n_iter_simulation, static_cast<size_t>(exporter->n_iter)) - 1;

  const size_t dump_interval = (n_iter_simulation) / (dump_number) + 1;

  size_t dump_counter = 0;
  double current_time = 0.;

  transitioner->update_flow(simulation, *current_reactor_state);

  simulation.setVolumes(current_reactor_state->gasVolume,
                        current_reactor_state->liquidVolume);

#pragma omp parallel default(shared) shared(current_time)
  {

    for (size_t __loop_counter = 0; __loop_counter < n_iter_simulation;
         ++__loop_counter)
    {

      DEBUG_INSTRUCTION

#pragma omp single
      {
        FILL_PAYLOAD;

        MPI_DISPATCH_MAIN;

        current_reactor_state =
            &transitioner->get_unchecked_mut();

        transitioner->update_flow(simulation, *current_reactor_state);
        simulation.setVolumes(current_reactor_state->gasVolume,
                              current_reactor_state->liquidVolume);
      }

      simulation.cycleProcess(d_t);

#pragma omp master
      {
        dump_counter++;
        if (dump_counter == dump_interval)
        {

          UPDATE_PROGRESS_BAR(n_iter_simulation, __loop_counter, true)
          exporter->append((double)__loop_counter * d_t,
                           simulation.getCliqData(),
                           simulation.mc_unit->domain.getDistribution());
          dump_counter = 0;
        }
      }

#pragma omp single
      {

        sync_step(exec, simulation);
        simulation.step(d_t, *current_reactor_state);
        sync_prepare_next(exec, simulation);
        current_time += d_t; // FIXME
      }
    }
  }
}