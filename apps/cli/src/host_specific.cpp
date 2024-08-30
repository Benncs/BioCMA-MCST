
#include <cma_read/reactorstate.hpp>
#include <common/common.hpp>
#include <common/execinfo.hpp>
#include <cstddef>
#include <dataexporter/factory.hpp>
#include <host_specific.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <mpi_w/wrap_mpi.hpp>
#include <post_process.hpp>
#include <rt_init.hpp>
#include <signal_handling.hpp>
#include <simulation/simulation.hpp>
#include <simulation/update_flows.hpp>
#include <stream_io.hpp>
#include <sync.hpp>
// #include <cereal/archives/binary.hpp>
// #include <cereal/types/memory.hpp>
// #include <cereal/types/vector.hpp>
// constexpr size_t n_particle_trigger_parralel = 1e6;

static void
main_loop(const SimulationParameters &params,
          const ExecInfo &exec,
          Simulation::SimulationUnit &simulation,
          std::unique_ptr<Simulation::FlowMapTransitioner> transitioner,
          std::unique_ptr<DataExporter> &exporter);

static constexpr size_t PROGRESS_BAR_WIDTH = 100;

inline static void update_progress_bar(size_t total, size_t current_position)
{
  const auto default_precision{std::cout.precision()};
  std::ios::sync_with_stdio(false);
  if constexpr (FlagCompileTIme::verbose)
  {
    std::string progressBar(PROGRESS_BAR_WIDTH, ' ');
    size_t progress = (current_position * PROGRESS_BAR_WIDTH) / total;

    for (size_t i = 0; i < progress; ++i)
    {
      progressBar[i] = '*';
    }

    std::cout << "Progress: [" << progressBar << "] " << std::fixed
              << std::setprecision(2)
              << (static_cast<float>(current_position) * 100.0 /
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

// In multirank context, fill the struct that will be broadcast to other workers
#define FILL_PAYLOAD                                                           \
  if constexpr (FlagCompileTIme::use_mpi)                                      \
  {                                                                            \
    mpi_payload.liquid_flows =                                                 \
        current_reactor_state->liquid_flow.getViewFlows().data();              \
    mpi_payload.liquid_volumes = current_reactor_state->liquidVolume;          \
    mpi_payload.gas_volumes = current_reactor_state->gasVolume;                \
    mpi_payload.neigbors =                                                     \
        current_reactor_state->liquid_flow.getViewNeighors().to_const();       \
  }

// In multirank context,  Send step payload to other workers
#define MPI_DISPATCH_MAIN                                                      \
  if constexpr (FlagCompileTIme::use_mpi)                                      \
  {                                                                            \
    for (size_t __macro_j = 1; __macro_j < exec.n_rank; ++__macro_j)           \
    {                                                                          \
      MPI_W::send(MPI_W::SIGNALS::RUN, __macro_j);                             \
      mpi_payload.send(__macro_j);                                             \
    }                                                                          \
  }

// In multirank context,  Send message to worker to ask them to return
#define SEND_MPI_SIG_STOP                                                      \
  if constexpr (FlagCompileTIme::use_mpi)                                      \
  {                                                                            \
    MPI_W::host_dispatch(exec, MPI_W::SIGNALS::STOP);                          \
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
                              initial_distribution,
                              simulation.mc_unit->init_weight);

    PostProcessing::save_initial_particle_state(simulation, data_exporter);
  }

  PostProcessing::show_sumup_state(simulation);

  main_loop(params, exec, simulation, std::move(transitioner), data_exporter);

  PostProcessing::show_sumup_state(simulation);

  SEND_MPI_SIG_STOP;
  if constexpr (FlagCompileTIme::use_mpi)
  {
    last_sync(exec, simulation);
  }

  PostProcessing::final_post_processing(
      exec, params, std::move(simulation), data_exporter);
}

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

  exporter->append(current_time,
                   simulation.getCliqData(),
                   simulation.mc_unit->domain.getDistribution(),
                   current_reactor_state->liquidVolume,
                   current_reactor_state->gasVolume);

  auto loop_functor = [&](auto &&local_container)
  {
    SignalHandler sig;
    auto result = local_container.get_extra();
    auto view_result = result.get_view();

    for (size_t __loop_counter = 0; __loop_counter < n_iter_simulation;
         ++__loop_counter)
    {

      DEBUG_INSTRUCTION

      transitioner->update_flow(simulation);

      current_reactor_state = transitioner->getState();

      FILL_PAYLOAD;

      MPI_DISPATCH_MAIN;

      transitioner->advance(simulation);

      simulation.cycleProcess(local_container, view_result, d_t);

      if (++dump_counter == dump_interval)
      {
        update_progress_bar(n_iter_simulation, __loop_counter);
        exporter->append(current_time,
                         simulation.getCliqData(),
                         simulation.mc_unit->domain.getDistribution(),
                         current_reactor_state->liquidVolume,
                         current_reactor_state->gasVolume);
        dump_counter = 0;
      }

      sync_step(exec, simulation);

      result.clear(local_container.n_particle());
      result.update_view(view_result);

      simulation.update_feed(d_t);
      simulation.step(d_t, *current_reactor_state);
      sync_prepare_next(exec, simulation);
      current_time += d_t;

      if (SignalHandler::is_usr1_raised())
      {
        PostProcessing::user_triggered_properties_export(simulation, exporter);
      }
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