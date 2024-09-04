#ifndef __HOST_POST_PROCESS_HPP__
#define __HOST_POST_PROCESS_HPP__

#include <common/execinfo.hpp>
#include <common/simulation_parameters.hpp>
#include <memory>


// Forward Declaration
class DataExporter;
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
                    Simulation::SimulationUnit &&simulation,
                    std::unique_ptr<DataExporter> &exporter);

  void show_sumup_state(Simulation::SimulationUnit &simulation);

  void save_initial_particle_state(Simulation::SimulationUnit &simulation,
                                   std::unique_ptr<DataExporter> &exporter);

  void user_triggered_properties_export(
      Simulation::SimulationUnit &sim,
      std::unique_ptr<DataExporter> &data_exporter);

  // get_particle_properties(unit,

  //                         aggregated_values,
  //                         spatial,
  //                         distribution.size(),
  //                         model_properties,
  //                         true);

} // namespace PostProcessing

#endif //__HOST_POST_PROCESS _HPP__