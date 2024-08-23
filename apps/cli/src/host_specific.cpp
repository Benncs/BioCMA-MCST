
#include "common/execinfo.hpp"
#include "mc/unit.hpp"
#include "rt_init.hpp"
#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include <cma_read/reactorstate.hpp>
#include <common/common.hpp>
#include <cstddef>
#include <dataexporter/factory.hpp>
#include <fstream>
#include <host_specific.hpp>
#include <memory>
#include <mpi_w/wrap_mpi.hpp>
#include <post_process.hpp>
#include <simulation/simulation.hpp>
#include <simulation/update_flows.hpp>
#include <sync.hpp>

// constexpr size_t n_particle_trigger_parralel = 1e6;

constexpr size_t PROGRESS_BAR_WIDTH = 100;

inline void
update_progress_bar(size_t total, size_t currentPosition, bool verbose)
{
  const auto default_precision{std::cout.precision()};
  std::ios::sync_with_stdio(false);
  if (verbose)
  {
    std::string progressBar(PROGRESS_BAR_WIDTH, ' ');
    size_t progress = (currentPosition * PROGRESS_BAR_WIDTH) / total;

    for (size_t i = 0; i < progress; ++i)
    {
      progressBar[i] = '*';
    }

    std::cout << "Progress: [" << progressBar << "] " << std::fixed
              << std::setprecision(2)
              << (static_cast<float>(currentPosition) * 100.0 /
                  static_cast<float>(total))
              << "%\r" << std::flush << std::setprecision(default_precision);
  }
}

#ifdef DEBUG
#  define DEBUG_INSTRUCTION // std::cout << "it: " <<__loop_counter<<
                            // std::endl;
#else
#  define DEBUG_INSTRUCTION
#endif

#define FILL_PAYLOAD                                                           \
  if constexpr (RT::use_mpi)                                                   \
  {                                                                            \
    mpi_payload.liquid_flows =                                                 \
        current_reactor_state->liquid_flow.getViewFlows().data();              \
    mpi_payload.liquid_volumes = current_reactor_state->liquidVolume;          \
    mpi_payload.gas_volumes = current_reactor_state->gasVolume;                \
    mpi_payload.neigbors =                                                     \
        current_reactor_state->liquid_flow.getViewNeighors().to_const();       \
  }
#define MPI_DISPATCH_MAIN                                                      \
  if constexpr (RT::use_mpi)                                                   \
  {                                                                            \
    for (size_t __macro_j = 1; __macro_j < exec.n_rank; ++__macro_j)           \
    {                                                                          \
      MPI_W::send(MPI_W::SIGNALS::RUN, __macro_j);                             \
      mpi_payload.send(__macro_j);                                             \
    }                                                                          \
  }

#define SEND_MPI_SIG_STOP                                                      \
  if constexpr (RT::use_mpi)                                                   \
  {                                                                            \
    host_dispatch(exec, MPI_W::SIGNALS::STOP);                                 \
  }

void host_process(
    const ExecInfo &exec,
    Simulation::SimulationUnit &&simulation,
    const SimulationParameters &params,
    std::unique_ptr<Simulation::FlowMapTransitioner> &&transitioner)
{
  std::unique_ptr<DataExporter> data_exporter;

  {
    auto initial_distribution = simulation.mc_unit->domain.getDistribution();

    data_exporter =
        data_exporter_factory(exec,
                              params,
                              params.results_file_name,
                              simulation.getDim(),
                              params.user_params.number_exported_result,
                              initial_distribution);

    PostProcessing::save_initial(simulation, data_exporter);
  }

  PostProcessing::show(simulation);

  main_loop(params, exec, simulation, std::move(transitioner), data_exporter);

  PostProcessing::show(simulation);

  SEND_MPI_SIG_STOP;
  if constexpr (RT::use_mpi)
  {
    last_sync(exec, simulation);
  }

  PostProcessing::post_process(
      exec, params, std::move(simulation), data_exporter);
}

// DEV

static bool save = false;

#include <csignal>
template <class Sim> void handle_sig(Sim &sim)
{

  if (save)
  {
    std::ofstream os("./out.cereal", std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
    archive(*sim.mc_unit);
    std::cout << "./out.cereal   " << sim.mc_unit->domain[0].n_cells
              << std::endl;
    int a;
  }
  save = false;
}

void handle_sig(int n)
{
  save = true;
}

// ENDDEV

void main_loop(const SimulationParameters &params,
               const ExecInfo &exec,
               Simulation::SimulationUnit &simulation,
               std::unique_ptr<Simulation::FlowMapTransitioner> transitioner,
               std::unique_ptr<DataExporter> &exporter)
{

  simulation.update_feed(0);

  // const size_t n_update_feed = 0; //TODO: move elsewhere

  const double d_t = params.d_t;

  const size_t n_iter_simulation = transitioner->get_n_timestep();

  const size_t dump_number =
      std::min(n_iter_simulation,
               static_cast<size_t>(exporter->expectecNiteration())) -
      1;

  const size_t dump_interval = (n_iter_simulation) / (dump_number) + 1;

  size_t dump_counter = 0;
  double current_time = 0.;
  // size_t update_feed_counter = 0;
  // const size_t update_feed_interval = (n_update_feed==0)? n_iter_simulation
  // : (n_iter_simulation) / (n_update_feed) + 1;

  MPI_W::HostIterationPayload mpi_payload;
  const auto *current_reactor_state = &transitioner->get_unchecked(0);

  transitioner->update_flow(simulation);
  signal(SIGUSR1, &handle_sig);

  exporter->append(current_time,
                   simulation.getCliqData(),
                   simulation.mc_unit->domain.getDistribution(),
                   current_reactor_state->liquidVolume,
                   current_reactor_state->gasVolume);

  auto loop_functor = [&](auto &&local_container)
  {
    for (size_t __loop_counter = 0; __loop_counter < n_iter_simulation;
         ++__loop_counter)
    {

      DEBUG_INSTRUCTION

      transitioner->update_flow(simulation);

      current_reactor_state = transitioner->getState();

      FILL_PAYLOAD;

      MPI_DISPATCH_MAIN;

      transitioner->advance(simulation);

      simulation.cycleProcess(container,d_t);

      dump_counter++;

      if (dump_counter == dump_interval)
      {
        update_progress_bar(n_iter_simulation, __loop_counter, true);
        exporter->append(current_time,
                         simulation.getCliqData(),
                         simulation.mc_unit->domain.getDistribution(),
                         current_reactor_state->liquidVolume,
                         current_reactor_state->gasVolume);
        dump_counter = 0;
      }

      if constexpr (RT::use_mpi)
      {
        sync_step(exec, simulation);
      }

      {

        simulation.update_feed(d_t);
        simulation.step(d_t, *current_reactor_state);
      }

      sync_prepare_next(exec, simulation);
      current_time += d_t;
      handle_sig(simulation);
    }
  };

  std::visit(loop_functor, simulation.mc_unit->container);

  exporter->append(current_time,
                   simulation.getCliqData(),
                   simulation.mc_unit->domain.getDistribution(),
                   current_reactor_state->liquidVolume,
                   current_reactor_state->gasVolume);

  transitioner.reset();
}