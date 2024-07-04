#include <common/execinfo.hpp>
#include <mc/unit.hpp>
#include <post_process.hpp>

#include <iostream>
#include <variant>



static void get_particle_properties(
    Simulation::SimulationUnit &simulation,
    std::unordered_map<std::string, std::vector<model_properties_t>>
        &aggregated_values,
    std::unordered_map<std::string, std::vector<double>> &spatial,
    size_t size);


void post_process(const ExecInfo &exec,
                  const SimulationParameters &params ,
                  Simulation::SimulationUnit &simulation,
                  std::unique_ptr<DataExporter> &exporter)
{

  if (exporter != nullptr)
  {
    std::unordered_map<std::string, std::vector<model_properties_t>>
        aggregated_values;

    std::unordered_map<std::string, std::vector<double>> spatial;
    auto distribution = simulation.mc_unit->domain.getDistribution();

    get_particle_properties(simulation,aggregated_values,spatial,distribution.size());

    exporter->write_final_results(simulation, distribution);

    exporter->write_final_particle_data(aggregated_values, spatial);
  }
}




void show(Simulation::SimulationUnit &simulation)
{

  std::vector<double> mass(simulation.mc_unit->domain.getNumberCompartments());
  double totmass = 0.;

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

  const auto &data = simulation.mc_unit->container.to_process.data();

  const auto n_particle = data.size();

  const auto prop_1 = model_properties(data.front());

  for (const auto &[key, _value] : prop_1)
  {
    aggregated_values[key].resize(n_particle);
    spatial[key].resize(size);
  }

#pragma omp parallel for shared(aggregated_values,                             \
                                    n_particle,                                \
                                    model_properties,                          \
                                    data,                                      \
                                    spatial,                                   \
                                    size,                                      \
                                    comp) default(none)
  for (size_t i = 0; i < n_particle; ++i)
  {
    auto prop = model_properties(data[i]);
    const size_t i_container = data[i].current_container;

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