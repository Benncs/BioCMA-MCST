#include "mc/particles/mcparticles.hpp"
#include <common/execinfo.hpp>
#include <mc/unit.hpp>
#include <post_process.hpp>

#include <iostream>
#include <stdexcept>
#include <variant>

namespace PostProcessing
{

  static void get_particle_properties(
      Simulation::SimulationUnit &simulation,
      std::unordered_map<std::string, std::vector<model_properties_t>>
          &aggregated_values,
      std::unordered_map<std::string, std::vector<double>> &spatial,
      size_t size);

  void save_initial(Simulation::SimulationUnit &simulation,
                    std::unique_ptr<DataExporter> &exporter)
  {
    if (exporter != nullptr)
    {
      std::unordered_map<std::string, std::vector<model_properties_t>>
          aggregated_values;

      std::unordered_map<std::string, std::vector<double>> spatial;
      get_particle_properties(simulation,
                              aggregated_values,
                              spatial,
                              std::get<1>(simulation.getDim()));

      exporter->write_initial_particle_data(aggregated_values, spatial);
    }
  }

  void post_process(const ExecInfo &exec,
                    const SimulationParameters &params,
                    Simulation::SimulationUnit &simulation,
                    std::unique_ptr<DataExporter> &exporter)
  {
    auto distribution = simulation.mc_unit->domain.getDistribution();


    auto tot = std::accumulate(
        distribution.begin(), distribution.end(), static_cast<size_t>(0));
    auto removed =
        simulation.mc_unit->ts_events[0].get<MC::EventType::Death>() +
        simulation.mc_unit->ts_events[0].get<MC::EventType::Exit>();
    auto new_p =
        simulation.mc_unit->ts_events[0].get<MC::EventType::NewParticle>();


    if (exporter != nullptr)
    {
      std::unordered_map<std::string, std::vector<model_properties_t>>
          aggregated_values;

      std::unordered_map<std::string, std::vector<double>> spatial;

      get_particle_properties(
          simulation, aggregated_values, spatial, distribution.size());

      exporter->write_final_results(simulation, distribution);

      // Results depending on simulation are exported or copied into properties,
      // clear montecarlo state to save memory
      simulation.clear_mc();

      exporter->write_final_particle_data(aggregated_values, spatial);
    }

    if(tot == (new_p - removed + params.user_params.numper_particle))
    {
      throw std::runtime_error("Results are not coherent (Bad particle balance)");
    }

    
  }

  void show(Simulation::SimulationUnit &simulation)
  {

    std::vector<double> mass(
        simulation.mc_unit->domain.getNumberCompartments());

    auto d = simulation.mc_unit->domain.getDistribution();

    for (auto &&i : d)
    {
      std::cout << i << " ";
      // count += i;
    }
    std::cout << '\n';
  }

  void get_particle_properties(
      Simulation::SimulationUnit &simulation,
      std::unordered_map<std::string, std::vector<model_properties_t>>
          &aggregated_values,
      std::unordered_map<std::string, std::vector<double>> &spatial,
      size_t size)
  {
    std::cout << "POST PROCESSING" << std::endl;
    const auto &model_properties = simulation.getModel().get_properties;

    const auto &comp = simulation.mc_unit->domain.data();

    const auto &particles_data =
        simulation.mc_unit->container.to_process.data();

    const auto n_particle = particles_data.size();

    size_t i_p = 0;
    while (particles_data[i_p].status != MC::CellStatus::IDLE &&
           i_p != particles_data.size())
    {
      i_p++;
    }
    const auto prop_1 = model_properties(particles_data[i_p]);

    for (const auto &[key, _value] : prop_1)
    {
      aggregated_values[key].resize(n_particle);
      spatial[key].resize(size);
    }

#pragma omp parallel for shared(aggregated_values,                             \
                                    n_particle,                                \
                                    model_properties,                          \
                                    particles_data,                            \
                                    spatial,                                   \
                                    size,                                      \
                                    comp) default(none),                       \
    firstprivate(i_p)
    for (size_t i = i_p; i < n_particle; ++i)
    {
      if (particles_data[i].status == MC::CellStatus::IDLE)
      {
        auto prop = model_properties(particles_data[i]);
        const size_t i_container = particles_data[i].current_container;

        for (const auto &[key, value] : prop)
        {
          aggregated_values[key][i] = value;
#pragma omp critical
          {
            if (const double *val = std::get_if<double>(&value))
            {
              spatial[key][i_container] +=
                  *val / static_cast<double>(comp[i_container].n_cells);
            }
          }
        }
      }
    }
  }

} // namespace PostProcessing