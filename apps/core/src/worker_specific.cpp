
#ifndef NO_MPI
#  include "biocma_cst_config.hpp"
#  include <csignal>
#  include <mpi_w/iteration_payload.hpp>
#  include <mpi_w/wrap_mpi.hpp>
#  include <impl_post_process.hpp>
#  include <simulation/simulation.hpp>
#  include <simulation/transitionner.hpp>
#  include <sync.hpp>
#  include <worker_specific.hpp>

void workers_process(const ExecInfo& exec,
                     Simulation::SimulationUnit& simulation,
                     const Core::SimulationParameters& params,
                     std::unique_ptr<Simulation::FlowMapTransitioner>&& transitioner,
                     Core::PartialExporter& partial_exporter)
{
  double d_t = params.d_t;
  size_t n_compartments = simulation.mc_unit->domain.getNumberCompartments();
  MPI_Status status;

  WrapMPI::IterationPayload payload(n_compartments * n_compartments, n_compartments);

  const auto loop_functor = [&](auto&& container)
  {
    // auto result = container.get_extra();
    // auto view_result = result.get_view();
    bool stop = false;
    WrapMPI::SIGNALS signal{};
    double current_time = 0;
    while (!stop)
    {
      signal = WrapMPI::try_recv<WrapMPI::SIGNALS>(0, &status);
      if (signal == WrapMPI::SIGNALS::STOP)
      {
        partial_exporter.write_number_particle(simulation.mc_unit->domain.getRepartition());
        if constexpr (AutoGenerated::FlagCompileTime::use_probe)
        {
          PostProcessing::save_probes(simulation, partial_exporter);
        }

        last_sync(exec, simulation);
        auto& list = container.get_compute();
        list.remove_dead(simulation.counter());
        PostProcessing::save_particle_state(simulation, partial_exporter);
        PostProcessing::reset_counter();
        stop = true;

        continue;
      }

      if (signal == WrapMPI::SIGNALS::DUMP)
      {
        PROFILE_SECTION("worker:dump")
        partial_exporter.write_number_particle(simulation.mc_unit->domain.getRepartition());

        // if constexpr (AutoGenerated::FlagCompileTime::use_probe)
        // {
        //   PostProcessing::save_probes(simulation, partial_exporter,true);
        // }

        if constexpr (AutoGenerated::FlagCompileTime::dump_particle_state)
        {
          PostProcessing::save_particle_state(simulation, partial_exporter);
        }
        continue;
      }

      payload.recv(0, &status);

      simulation.mc_unit->domain.setLiquidNeighbors(payload.neighbors);
      transitioner->update_flow(simulation, payload.liquid_flows, n_compartments);
      transitioner->advance(simulation);

      simulation.setVolumes(payload.gas_volumes, payload.liquid_volumes);

      simulation.cycleProcess(container, d_t);

      // result.clear(container.n_particle());
      // result.update_view(view_result);
      simulation.update_feed(current_time, d_t);
      current_time += d_t;

      sync_step(exec, simulation);

      sync_prepare_next(simulation);
    }
  };
  std::visit(loop_functor, simulation.mc_unit->container);
}
#endif // NO_MPI
