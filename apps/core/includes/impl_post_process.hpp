#ifndef __HOST_POST_PROCESS_HPP__
#define __HOST_POST_PROCESS_HPP__

#include "dataexporter/main_exporter.hpp"
#include "dataexporter/partial_exporter.hpp"
#include <common/execinfo.hpp>
#include <core/simulation_parameters.hpp>
#include <memory>

namespace Simulation
{
  class SimulationUnit;
} // namespace Simulation

namespace PostProcessing
{
  void save_results(const ExecInfo& exec,
                    const Core::SimulationParameters& params,
                    Simulation::SimulationUnit& simulation);

  void final_post_processing(const ExecInfo& exec,
                             const Core::SimulationParameters& params,
                             Simulation::SimulationUnit& simulation,
                             std::unique_ptr<Core::MainExporter>& mde);

  void show_sumup_state(const Simulation::SimulationUnit& simulation) noexcept;

  void save_particle_state(Simulation::SimulationUnit& simulation, Core::PartialExporter& pde);

  void save_probes(Simulation::SimulationUnit& simulation,
                   Core::PartialExporter& pde,
                   bool force = false);

  void reset_counter();
  // get_particle_properties(unit,

  //                         aggregated_values,
  //                         spatial,
  //                         distribution.size(),
  //                         model_properties,
  //                         true);

} // namespace PostProcessing

#endif //__HOST_POST_PROCESS _HPP__
