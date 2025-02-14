
#include <algorithm> //for std::min
#include <biocma_cst_config.hpp>
#include <cma_read/reactorstate.hpp>
#include <common/common.hpp>
#include <common/execinfo.hpp>
#include <core/simulation_parameters.hpp>
#include <cstddef>
#include <dataexporter/main_exporter.hpp>
#include <dataexporter/partial_exporter.hpp>
#include <host_specific.hpp>
#include <impl_post_process.hpp>
#include <iomanip>
#include <ios>
#include <mc/unit.hpp>
#include <memory>
#include <optional>
#include <ostream>
#include <signal_handling.hpp>
#include <simulation/simulation.hpp>
#include <simulation/transitionner.hpp>
#include <string>
#include <sync.hpp>
#include <utility>
#include <variant>
#ifndef NO_MPI
#  include <mpi_w/wrap_mpi.hpp>
#endif

namespace
{

  // constexpr size_t n_particle_trigger_parralel = 1e6;
  void main_loop(const Core::SimulationParameters& params,
                 const ExecInfo& exec,
                 Simulation::SimulationUnit& simulation,
                 std::unique_ptr<Simulation::FlowMapTransitioner> transitioner,
                 std::unique_ptr<Core::MainExporter>& main_exporter,
                 Core::PartialExporter& partial_exporter);

  constexpr size_t PROGRESS_BAR_WIDTH = 100;

  void handle_export(const ExecInfo& exec,
                     size_t& dump_counter,
                     size_t dump_interval,
                     size_t n_iter_simulation,
                     size_t __loop_counter,
                     double current_time,
                     Simulation::SimulationUnit& simulation,
                     auto& current_reactor_state,
                     std::unique_ptr<Core::MainExporter>& main_exporter,
                     Core::PartialExporter& partial_exporter);

  inline void update_progress_bar(size_t total, size_t current_position)
  {

    if constexpr (AutoGenerated::FlagCompileTime::verbose)
    {
      const int default_precision = std::cout.precision(); // NOLINT Conversion long to int
      std::ios::sync_with_stdio(false);
      std::string progressBar(PROGRESS_BAR_WIDTH, ' ');
      size_t progress = (current_position * PROGRESS_BAR_WIDTH) / total;

      for (size_t i = 0; i < progress; ++i)
      {
        progressBar[i] = '*';
      }

      std::cout << "Progress: [" << progressBar << "] " << std::fixed << std::setprecision(2)
                << (static_cast<float>(current_position) * 100.0 / static_cast<float>(total))
                << "%\r" << std::flush << std::setprecision(default_precision);
    }
  }

  std::unique_ptr<Core::MainExporter> make_main_exporter(const ExecInfo& exec,
                                                         const Core::SimulationParameters& params)
  {
    std::unique_ptr<Core::MainExporter> main_exporter =
        std::make_unique<Core::MainExporter>(exec, params.results_file_name + std::string(".h5"));

    for (std::size_t i_rank = 0; i_rank < exec.n_rank; ++i_rank)
    {
      std::string group = "files/" + std::to_string(i_rank);

      auto filename = params.results_file_name + "_partial_" + std::to_string(i_rank) + ".h5";

      main_exporter->do_link(filename, group, "/");
    }

    return main_exporter;
  }

} // namespace

#ifndef NO_MPI
#  define BARRIER WrapMPI::barrier();
// In multirank context, fill the struct that will be broadcast to other workers
#  define FILL_PAYLOAD                                                                             \
    if constexpr (AutoGenerated::FlagCompileTime::use_mpi)                                         \
    {                                                                                              \
      mpi_payload.fill(current_reactor_state);                                                     \
    }
// In multirank context,  Send step payload to other workers
#  define MPI_DISPATCH_MAIN                                                                        \
    if constexpr (AutoGenerated::FlagCompileTime::use_mpi)                                         \
    {                                                                                              \
      for (size_t __macro_j = 1; __macro_j < exec.n_rank; ++__macro_j)                             \
      {                                                                                            \
        auto _ = WrapMPI::send(WrapMPI::SIGNALS::RUN, __macro_j);                                  \
        _ = mpi_payload.send(__macro_j);                                                           \
      }                                                                                            \
    }
#  define SEND_MPI_SIG_STOP                                                                        \
    if constexpr (AutoGenerated::FlagCompileTime::use_mpi)                                         \
    {                                                                                              \
      WrapMPI::host_dispatch(exec, WrapMPI::SIGNALS::STOP);                                        \
    }
#  define SEND_MPI_SIG_DUMP                                                                        \
    if constexpr (AutoGenerated::FlagCompileTime::use_mpi)                                         \
    {                                                                                              \
      WrapMPI::host_dispatch(exec, WrapMPI::SIGNALS::DUMP);                                        \
    }
#  define INIT_PAYLOAD WrapMPI::HostIterationPayload mpi_payload;
#else
// Without MPI do nothing
#  define BARRIER
#  define FILL_PAYLOAD
#  define MPI_DISPATCH_MAIN
#  define SEND_MPI_SIG_STOP
#  define SEND_MPI_SIG_DUMP
#  define INIT_PAYLOAD
#endif

#ifdef DEBUG
#  define DEBUG_INSTRUCTION // std::cout << "it: " <<__loop_counter<<std::endl;
#else
#  define DEBUG_INSTRUCTION
#endif

void host_process(const ExecInfo& exec,
                  Simulation::SimulationUnit& simulation,
                  const Core::SimulationParameters& params,
                  std::unique_ptr<Simulation::FlowMapTransitioner>&& transitioner,
                  Core::PartialExporter& partial_exporter)
{

  std::unique_ptr<Core::MainExporter> main_exporter = make_main_exporter(exec, params);

  const auto [n_species, n_compartment] = simulation.getDimensions();

  main_exporter->init_fields(
      params.number_exported_result, n_compartment, n_species, simulation.two_phase_flow());

  main_exporter->write_initial(
      simulation.mc_unit->init_weight, params, simulation.mc_unit->domain.getRepartition());

  PostProcessing::show_sumup_state(simulation);

  main_loop(params, exec, simulation, std::move(transitioner), main_exporter, partial_exporter);

  PostProcessing::show_sumup_state(simulation);

  SEND_MPI_SIG_STOP;

  auto clean_list = [&simulation](auto&& container)
  {
    auto& list = container.get_compute();
    if (simulation.counter() != 0)
    {
      list.remove_dead(simulation.counter());
    }
  };

  std::visit(clean_list, simulation.mc_unit->container);
  PostProcessing::save_particle_state(simulation, partial_exporter);
  last_sync(exec, simulation);
  PostProcessing::final_post_processing(exec, params, simulation, main_exporter);

  PostProcessing::reset_counter();
}

namespace
{

  void main_loop(const Core::SimulationParameters& params,
                 const ExecInfo& exec,
                 Simulation::SimulationUnit& simulation,
                 std::unique_ptr<Simulation::FlowMapTransitioner> transitioner,
                 std::unique_ptr<Core::MainExporter>& main_exporter,
                 Core::PartialExporter& partial_exporter)
  {

    simulation.update_feed(0, 0);

    // const size_t n_update_feed = 0; //TODO: move elsewhere

    const double d_t = params.d_t;

    const size_t n_iter_simulation = transitioner->get_n_timestep();

    const size_t dump_number =
        std::min(n_iter_simulation, static_cast<size_t>(params.number_exported_result)) - 1;
    // FIXME when number_exported_result==0 and number_exported_result==1
    const size_t dump_interval = (params.number_exported_result != 0)
                                     ? (n_iter_simulation) / (dump_number) + 1
                                     : n_iter_simulation + 1;

    size_t dump_counter = 0;
    double current_time = simulation.get_start_time_mut();
    // size_t update_feed_counter = 0;
    // const size_t update_feed_interval = (n_update_feed==0)? n_iter_simulation
    // : (n_iter_simulation) / (n_update_feed) + 1;

    INIT_PAYLOAD

    const auto* current_reactor_state = &transitioner->get_unchecked(0);

    transitioner->update_flow(simulation);

    // FIX CMTOOL
    const auto gas_concentration = simulation.getCgasData();
    auto gas_volume = (!gas_concentration.has_value())
                          ? std::nullopt
                          : std::make_optional(current_reactor_state->gasVolume);

    main_exporter->update_fields(current_time,
                                 simulation.getCliqData(),
                                 current_reactor_state->liquidVolume,
                                 gas_concentration,
                                 gas_volume,
                                 simulation.getMTRData());

    auto loop_functor = [&](auto&& local_container)
    {
      Core::SignalHandler sig;

      for (size_t __loop_counter = 0; __loop_counter < n_iter_simulation; ++__loop_counter)
      {

        DEBUG_INSTRUCTION

        {
          PROFILE_SECTION("host:update_flow")
          transitioner->update_flow(simulation);
        }

        current_reactor_state = transitioner->getState();

        FILL_PAYLOAD;

        MPI_DISPATCH_MAIN;

        {
          PROFILE_SECTION("host:update_flow::advance")
          transitioner->advance(simulation);
        }

        simulation.cycleProcess(local_container, d_t);

        handle_export(exec,
                      dump_counter,
                      dump_interval,
                      n_iter_simulation,
                      __loop_counter,
                      current_time,
                      simulation,
                      current_reactor_state,
                      main_exporter,
                      partial_exporter);

        sync_step(exec, simulation);
        {
          PROFILE_SECTION("host:sync_update")
          simulation.update_feed(current_time, d_t);
          simulation.step(d_t, *current_reactor_state);
          current_time += d_t;
        }
        sync_prepare_next(simulation);

        if (Core::SignalHandler::is_usr1_raised())
        {
          PostProcessing::save_particle_state(simulation, partial_exporter);
        }

        if (Core::SignalHandler::is_sigint_raised())
        {
          std::cout << "User interruption" << std::endl;
          break;
        }
      }
    };

    std::visit(loop_functor, simulation.mc_unit->container);

    gas_volume = (!simulation.getCgasData().has_value())
                     ? std::nullopt
                     : std::make_optional(current_reactor_state->gasVolume);

    main_exporter->update_fields(current_time,
                                 simulation.getCliqData(),
                                 current_reactor_state->liquidVolume,
                                 simulation.getCgasData(),
                                 gas_volume,
                                 simulation.getMTRData());
    partial_exporter.write_number_particle(simulation.mc_unit->domain.getRepartition());

    // if constexpr (AutoGenerated::FlagCompileTime::use_probe)
    // {
    //   PostProcessing::save_probes(simulation, partial_exporter,true);
    // }
    simulation.get_end_time_mut() = current_time;
    transitioner.reset();
  }

  // NOLINTBEGIN
  void handle_export(const ExecInfo& exec, /*/Exec is used for MPI **/
                     size_t& dump_counter,
                     const size_t dump_interval,
                     const size_t n_iter_simulation,
                     const size_t __loop_counter,
                     const double current_time,
                     Simulation::SimulationUnit& simulation,
                     auto& current_reactor_state,
                     std::unique_ptr<Core::MainExporter>& main_exporter,
                     Core::PartialExporter& partial_exporter)
  // NOLINTEND
  {

    PROFILE_SECTION("host:handle_export")
    if (++dump_counter == dump_interval)
    {
      SEND_MPI_SIG_DUMP
      auto vg = (simulation.getCgasData().has_value())
                    ? std::make_optional(current_reactor_state->gasVolume)
                    : std::nullopt;

      update_progress_bar(n_iter_simulation, __loop_counter);
      main_exporter->update_fields(current_time,
                                   simulation.getCliqData(),
                                   current_reactor_state->liquidVolume,
                                   simulation.getCgasData(),
                                   vg,
                                   simulation.getMTRData());

      partial_exporter.write_number_particle(simulation.mc_unit->domain.getRepartition());

      if constexpr (AutoGenerated::FlagCompileTime::use_probe)
      {
        PostProcessing::save_probes(simulation, partial_exporter);
      }
      if constexpr (AutoGenerated::FlagCompileTime::dump_particle_state)
      {
        PostProcessing::save_particle_state(simulation, partial_exporter);
      }

      dump_counter = 0;
    }
  }

} // namespace
