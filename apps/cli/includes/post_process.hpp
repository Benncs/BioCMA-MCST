#ifndef __HOST_POST_PROCESS_HPP__
#define __HOST_POST_PROCESS_HPP__

#include <dataexporter/data_exporter.hpp>
#include <simulation/simulation.hpp>

namespace PostProcessing
{
  void save_results(const ExecInfo &exec,
                    const SimulationParameters &params,
                    Simulation::SimulationUnit &simulation);

  void post_process(const ExecInfo &exec,
                    const SimulationParameters &params,
                    Simulation::SimulationUnit &&simulation,
                    std::unique_ptr<DataExporter> &exporter);

  void show(Simulation::SimulationUnit &simulation);

  void save_initial(Simulation::SimulationUnit &simulation,
                    std::unique_ptr<DataExporter> &exporter);

  void append_properties(int counter,
                         Simulation::SimulationUnit &simulation,
                         std::unique_ptr<DataExporter> &exporter);

  
  void user_triggered_properties_export(Simulation::SimulationUnit &sim, std::unique_ptr<DataExporter> &data_exporter);

  // get_particle_properties(unit,

  //                         aggregated_values,
  //                         spatial,
  //                         distribution.size(),
  //                         model_properties,
  //                         true);

} // namespace PostProcessing

#endif //__HOST_POST_PROCESS _HPP__