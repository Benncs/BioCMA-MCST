#include <csignal>
#include <mpi_w/iteration_payload.hpp>
#include <mpi_w/wrap_mpi.hpp>
#include <simulation/simulation.hpp>
#include <simulation/update_flows.hpp>
#include <sync.hpp>
#include <worker_specific.hpp>






void workers_process(
    const ExecInfo &exec,
    Simulation::SimulationUnit &&simulation,
    const SimulationParameters &params,
    std::unique_ptr<Simulation::FlowMapTransitioner> &&transitioner)
{
  double d_t = params.d_t;
  size_t n_compartments = simulation.mc_unit->domain.getNumberCompartments();
  MPI_Status status;

  MPI_W::IterationPayload payload(n_compartments * n_compartments,
                                  n_compartments);

  const auto loop_functor = [&](auto &&container)
  {
    auto result = container.get_extra();
    auto view_result = result.get_view();
    bool stop = false;
    MPI_W::SIGNALS signal{};
    while (!stop)
    {

      signal = MPI_W::try_recv<MPI_W::SIGNALS>(0, &status);
      if (signal == MPI_W::SIGNALS::STOP)
      {
        last_sync(exec, simulation);
        stop = true;
        continue;
      }

      payload.recv(0, &status);

      simulation.mc_unit->domain.setLiquidNeighbors(payload.neigbors);
      transitioner->update_flow(
          simulation, payload.liquid_flows, n_compartments);
      transitioner->advance(simulation);

      simulation.setVolumes(payload.gas_volumes, payload.liquid_volumes);

      simulation.cycleProcess(container, view_result, d_t);
      result.clear(container.n_particle());
      result.update_view(view_result);

      sync_step(exec, simulation);
      sync_prepare_next(exec, simulation);
    }
  };
  std::visit(loop_functor, simulation.mc_unit->container);
}