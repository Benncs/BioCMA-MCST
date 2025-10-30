#ifndef __CORE_HOST_EXPORT_HANDLER_HPP__
#define __CORE_HOST_EXPORT_HANDLER_HPP__

#include <dataexporter/main_exporter.hpp>
#include <dataexporter/partial_exporter.hpp>
#include <progress_bar.hpp>
#include <simulation/simulation.hpp>
#include <transitionner/transitionner.hpp>

class ExportHandler
{

public:
  ExportHandler() = default;

  ExportHandler(std::shared_ptr<Core::MainExporter> _main_exporter,
                ExecInfo _exec,
                std::size_t _dump_interval,
                std::size_t _n_iter_simulation);

  /**
   * @brief Handles periodic export of simulation data.
   *
   * This operator is called during the simulation loop to periodically export
   * simulation state.
   * It updates the main exporter, writes particle counts, and optionally saves
   * data based on compile-time flag.
   *
   * @param current_time The current simulation time.
   * @param loop_counter The current loop iteration counter.
   * @param simulation The simulation unit containing the state and data to
   * export.
   * @param partial_exporter The partial exporter for writing particle and probe
   * data.
   */
  void operator()(double current_time,
                  size_t loop_counter,
                  Simulation::SimulationUnit& simulation,
                  Core::PartialExporter& partial_exporter);

  /**
   * @brief Prepares and updates the exporter with the current simulation state
   * before/after main loop.
   *
   * This function gathers the necessary data from the simulation and reactor
   * state and updates the main exporter with this information.
   *
   * @param current_time The current simulation time.
   * @param simulation The simulation unit containing the state and data to
   * export.
   * @param transitioner The flow map transitioner providing access to the
   * reactor state.
   */
  void pre_post_export(
      double current_time,
      const Simulation::SimulationUnit& simulation,
      const std::unique_ptr<CmaUtils::FlowMapTransitionner>& transitioner);

private:
  size_t dump_counter{};
  size_t dump_interval{};
  size_t n_iter_simulation{};
  [[maybe_unused]] ExecInfo exec{};
  std::shared_ptr<Core::MainExporter> main_exporter;

  /**
   * @brief Prepares the event span if the event counter is enabled.
   * @param simulation The simulation unit.
   * @return Optional span of event counts.
   */
  static std::optional<std::span<std::size_t>>
  prepareEventSpan(Simulation::SimulationUnit& simulation);

  IO::ProgressBar progressbar;
};

#endif
