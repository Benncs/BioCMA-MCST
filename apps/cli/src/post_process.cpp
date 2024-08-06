#include "mc/particles/mcparticles.hpp"
#include "models/types.hpp"
#include <common/execinfo.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <post_process.hpp>

#include <iostream>
#include <stdexcept>
#include <variant>

namespace PostProcessing
{
  // static void get_particle_properties2(
  //     std::unique_ptr<DataExporter> &exporter,
  //     std::unique_ptr<MC::MonteCarloUnit> &mc_unit,
  //     size_t size,
  //     std::unordered_map<std::string, std::vector<double>> &spatial,
  //     ModelGetProperties model_properties,
  //     bool clean);

  static void get_particle_properties(
      std::unique_ptr<MC::MonteCarloUnit> &mc_unit,
      std::unordered_map<std::string, std::vector<model_properties_t>>
          &aggregated_values,
      std::unordered_map<std::string, std::vector<double>> &spatial,
      size_t size,
      ModelGetProperties model_properties,
      bool clean = false);

  void save_initial(Simulation::SimulationUnit &simulation,
                    std::unique_ptr<DataExporter> &exporter)
  {
    if (exporter != nullptr)
    {
      std::unordered_map<std::string, std::vector<model_properties_t>>
          aggregated_values;

      std::unordered_map<std::string, std::vector<double>> spatial;
      get_particle_properties(simulation.mc_unit,
                              aggregated_values,
                              spatial,
                              std::get<1>(simulation.getDim()),
                              simulation.getModel().get_properties);

      exporter->write_initial_particle_data(aggregated_values, spatial);
    }
  }

  void post_process(const ExecInfo &exec,
                    const SimulationParameters &params,
                    Simulation::SimulationUnit &&simulation,
                    std::unique_ptr<DataExporter> &exporter)
  {
    auto distribution = simulation.mc_unit->domain.getDistribution();

    auto tot = std::accumulate(
        distribution.begin(), distribution.end(), static_cast<size_t>(0));
    auto removed = simulation.mc_unit->events.get<MC::EventType::Death>() +
                   simulation.mc_unit->events.get<MC::EventType::Exit>();
    auto new_p = simulation.mc_unit->events.get<MC::EventType::NewParticle>();

    if (exporter != nullptr)
    {

      std::unordered_map<std::string, std::vector<model_properties_t>>
          aggregated_values;

      std::unordered_map<std::string, std::vector<double>> spatial;

      exporter->write_final_results(simulation, distribution);

      // Results depending on simulation are exported or copied into properties,
      // clear montecarlo state to save memory
      auto unit = std::move(simulation.mc_unit);

      const auto model_properties = simulation.getModel().get_properties;

      get_particle_properties(unit,

                              aggregated_values,
                              spatial,
                              distribution.size(),
                              model_properties,
                              true);
      unit.reset();
      exporter->write_final_particle_data(aggregated_values, spatial);
    }

    if (tot != (new_p - removed + params.user_params.numper_particle))
    {
      std::cerr << ("Results are not coherent (Bad particle balance): ");
      std::cerr << tot << "=" << new_p << "-" << removed << "+"
                << params.user_params.numper_particle << std::endl;
      ;
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
      std::unique_ptr<MC::MonteCarloUnit> &mc_unit,
      std::unordered_map<std::string, std::vector<model_properties_t>>
          &aggregated_values,
      std::unordered_map<std::string, std::vector<double>> &spatial,
      size_t size,
      const ModelGetProperties model_properties,
      bool clean)
  {
    auto find_first_idle_particle = [](const auto &particles_data)
    {
      size_t i_p = 0;
      while (particles_data[i_p].status != MC::CellStatus::IDLE &&
             i_p != particles_data.size())
      {
        i_p++;
      }
      return i_p;
    };

    std::cout << "POST PROCESSING" << std::endl;

    const auto compartments = mc_unit->domain.data();

    auto particles_data = mc_unit->container.to_process.data_span();

    const auto n_particle = particles_data.size();

    const size_t i_p = find_first_idle_particle(particles_data);

    const auto first_property = model_properties(particles_data[i_p]);

    std::for_each(first_property.begin(),
                  first_property.end(),
                  [&](auto &&_tuple)
                  {
                    const auto &[key, _value] = _tuple;
                    aggregated_values[key].resize(n_particle);
                    spatial[key].resize(size);
                  });

#pragma omp parallel for shared(aggregated_values,                             \
                                    n_particle,                                \
                                    model_properties,                          \
                                    size,                                      \
                                    spatial) default(none),                    \
    firstprivate(i_p, particles_data, compartments, clean)

    for (size_t i = i_p; i < n_particle; ++i)
    {
      auto &particle = particles_data[i];
      if (particle.status == MC::CellStatus::IDLE)
      {
        auto prop = model_properties(particle);
        const size_t i_container = particle.current_container;

        for (const auto &[key, value] : prop)
        {
            aggregated_values[key][i] = value;

          if (const double *val = std::get_if<double>(&value))
          {
            // double weighted_val = *val*particle.weight;
            
#pragma omp atomic
            spatial[key][i_container] +=
                *val / static_cast<double>(compartments[i_container].n_cells);
          }
        }
        if (clean)
        {
          particle.clearState();
        }
      }
    }
  }

  //   void get_particle_properties2(
  //       std::unique_ptr<DataExporter> &exporter,
  //       std::unique_ptr<MC::MonteCarloUnit> &mc_unit,
  //       size_t size,
  //       std::unordered_map<std::string, std::vector<double>> &spatial,
  //       const ModelGetProperties model_properties,
  //       bool clean)
  //   {
  //     auto find_first_idle_particle = [](const auto &particles_data)
  //     {
  //       size_t i_p = 0;
  //       while (particles_data[i_p].status != MC::CellStatus::IDLE &&
  //              i_p != particles_data.size())
  //       {
  //         i_p++;
  //       }
  //       return i_p;
  //     };

  //     std::cout << "POST PROCESSING" << std::endl;

  //     const auto compartments = mc_unit->domain.data();

  //     auto particles_data = mc_unit->container.to_process.data_span();

  //     const auto n_particle = particles_data.size();

  //     const size_t i_p = find_first_idle_particle(particles_data);

  //     const auto first_property = model_properties(particles_data[i_p]);

  //     auto *handle = exporter->start_model_dataset();

  //     std::for_each(first_property.begin(),
  //                   first_property.end(),
  //                   [&](auto &&_tuple)
  //                   {
  //                     const auto &[key, _value] = _tuple;
  //                     exporter->init_fill_model_dataset(
  //                         handle, "biological_model/final2/" + key,
  //                         n_particle);
  //                   });

  //     // using model_properties_t = std::variant<double, int, std::string>;

  //     // using model_properties_detail_t =
  //     //     std::unordered_map<std::string, model_properties_t>;

  //     // for (size_t idx = i_p; idx < n_particle; ++idx)
  //     // {
  //     //   auto &particle = particles_data[idx];
  //     //   if (particle.status == MC::CellStatus::IDLE)
  //     //   {
  //     //     model_properties_detail_t properties =
  //     model_properties(particle);

  //     //   }
  //     // }

  // #pragma omp parallel for
  //     for (size_t idx = i_p; idx < n_particle; ++idx)
  //     {
  //       auto &particle = particles_data[idx];
  //       if (particle.status == MC::CellStatus::IDLE)
  //       {
  //         auto properties = model_properties(particle);
  //         const size_t container_index = particle.current_container;

  //         for (auto &&property : properties)
  //         {
  //           auto [key, value_variant] = property;
  //           std::visit(
  //               [&](auto &&sample_val)
  //               {
  //                 using T = std::decay_t<decltype(sample_val)>;
  //                 if constexpr (std::is_same_v<T, double>)
  //                 {
  //                   double value = sample_val;
  // #pragma omp critical
  //                   {
  //                     exporter->fill_model_dataset(
  //                         idx, handle, "biological_model/final2/" + key,
  //                         value);
  //                   }
  //                 }
  //               },
  //               value_variant);
  //         }

  //         if (clean)
  //         {
  //           particle.clearState();
  //         }
  //       }
  //     }

  // #pragma omp barrier
  //     exporter->stop_fill_model_dataset(handle);
  //   }

} // namespace PostProcessing