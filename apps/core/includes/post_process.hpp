#ifndef __HOST_POST_PROCESS_HPP__
#define __HOST_POST_PROCESS_HPP__

#include "dataexporter/main_exporter.hpp"
#include "dataexporter/partial_exporter.hpp"
#include <common/execinfo.hpp>
#include <common/simulation_parameters.hpp>
#include <memory>


namespace Simulation
{
  class SimulationUnit;
} // namespace Simulation

namespace PostProcessing
{
  void save_results(const ExecInfo &exec,
                    const SimulationParameters &params,
                    Simulation::SimulationUnit &simulation);

  void final_post_processing(const ExecInfo &exec,
                    const SimulationParameters &params,
                    Simulation::SimulationUnit &&simulation,std::unique_ptr<CORE_DE::MainExporter>& mde);

  void show_sumup_state(Simulation::SimulationUnit &simulation);

  void save_initial_particle_state(Simulation::SimulationUnit &simulation,CORE_DE::PartialExporter& pde);

  void save_final_particle_state(Simulation::SimulationUnit &simulation,CORE_DE::PartialExporter& pde);

  void save_probes(Simulation::SimulationUnit &simulation,CORE_DE::PartialExporter& pde);

  void user_triggered_properties_export(
      Simulation::SimulationUnit &sim,CORE_DE::PartialExporter& pde);

  // get_particle_properties(unit,

  //                         aggregated_values,
  //                         spatial,
  //                         distribution.size(),
  //                         model_properties,
  //                         true);

} // namespace PostProcessing

#endif //__HOST_POST_PROCESS _HPP__