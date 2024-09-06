#include <Kokkos_Printf.hpp>
#include <common/execinfo.hpp>
#include <dataexporter/data_exporter.hpp>
#include <iostream>
#include <mc/unit.hpp>
#include <memory>
#include <post_process.hpp>
#include <simulation/simulation.hpp>
#include <variant>


namespace PostProcessing
{
  static void append_properties(int counter,
                                Simulation::SimulationUnit &simulation,
                                std::unique_ptr<DataExporter> &exporter);

  static void get_particle_properties(
      std::unique_ptr<MC::MonteCarloUnit> &mc_unit,
      std::unordered_map<std::string, std::vector<model_properties_t>>
          &aggregated_values,
      std::unordered_map<std::string, std::vector<double>> &spatial,
      size_t size,
      bool clean = true);

  void save_initial_particle_state(Simulation::SimulationUnit &simulation,
                                   std::unique_ptr<DataExporter> &exporter)
  {
    if (exporter != nullptr)
    {
      std::unordered_map<std::string, std::vector<model_properties_t>>
          aggregated_values;

      std::unordered_map<std::string, std::vector<double>> spatial;
      auto distribution = simulation.mc_unit->domain.getRepartition();
      std::cout << "EXPORTING PARTICLE DATA" << std::endl;
      get_particle_properties(simulation.mc_unit,
                              aggregated_values,
                              spatial,
                              distribution.size(),
                              false);

      exporter->write_initial_particle_data(aggregated_values, spatial);
    }
  }

  void final_post_processing(const ExecInfo &exec,
                             const SimulationParameters &params,
                             Simulation::SimulationUnit &&simulation,
                             std::unique_ptr<DataExporter> &exporter)
  {
        std::cout << "POST PROCESSING" << std::endl;
    auto distribution = simulation.mc_unit->domain.getRepartition();

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
      std::cout << "EXPORTING PARTICLE DATA" << std::endl;
      get_particle_properties(
          unit, aggregated_values, spatial, distribution.size(), true);
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

  void show_sumup_state(Simulation::SimulationUnit &simulation)
  {
    // Assuming domain data is in sharedSpace
    for (auto &&c : simulation.mc_unit->domain)
    {
      Kokkos::printf("%d ", c.n_cells);
    }
    Kokkos::printf("\r\n");
  }

  void append_properties(int counter,
                         Simulation::SimulationUnit &simulation,
                         std::unique_ptr<DataExporter> &exporter)
  {
    if (exporter != nullptr)
    {
      std::unordered_map<std::string, std::vector<model_properties_t>>
          aggregated_values;

      std::unordered_map<std::string, std::vector<double>> spatial;
      auto distribution = simulation.mc_unit->domain.getRepartition();
      get_particle_properties(simulation.mc_unit,
                              aggregated_values,
                              spatial,
                              distribution.size(),
                              false);
      exporter->append_particle_properties(counter, aggregated_values, spatial);
    }
  }

  void
  user_triggered_properties_export(Simulation::SimulationUnit &sim,
                                   std::unique_ptr<DataExporter> &data_exporter)
  {
    static int counter = 0;
    counter++;
    PostProcessing::append_properties(counter, sim, data_exporter);
    // std::ofstream os("./out.cereal", std::ios::binary);
    // cereal::BinaryOutputArchive archive(os);
    // archive(*sim.mc_unit);
    // std::cout << "./out.cereal   " << sim.mc_unit->domain[0].n_cells
    //           << std::endl;
  }

  void get_particle_properties(
      std::unique_ptr<MC::MonteCarloUnit> &mc_unit,
      std::unordered_map<std::string, std::vector<model_properties_t>>
          &aggregated_values,
      std::unordered_map<std::string, std::vector<double>> &spatial,
      size_t size,
      bool clean)
  {
    auto find_first_idle_particle = [](const auto &particles_data)
    {
      size_t i_p = 0;
      while (particles_data[i_p].properties.status != MC::CellStatus::IDLE &&
             i_p != particles_data.size())
      {
        i_p++;
      }
      return i_p;
    };



    const auto compartments = mc_unit->domain.data();

    auto functor = [&](auto &&container)
    {
      auto particles_data = container.get_host();

      const auto n_particle = particles_data.size();

      const size_t i_p = find_first_idle_particle(particles_data);

      const auto first_property = particles_data[i_p].data.get_properties();

      std::for_each(first_property.begin(),
                    first_property.end(),
                    [&](auto &&_tuple)
                    {
                      const auto &[key, _value] = _tuple;
                      aggregated_values[key].resize(n_particle);
                      spatial[key].resize(size);
                    });

      for (size_t i = i_p; i < n_particle; ++i)
      {
        auto &particle = particles_data[i];
        if (particle.properties.status == MC::CellStatus::IDLE)
        {
          auto prop = particle.data.get_properties();
          const size_t i_container = particle.properties.current_container;

          for (const auto &[key, value] : prop)
          {
            aggregated_values[key][i] = value;

            if (const double *val = std::get_if<double>(&value))
            {
              // double weighted_val = *val*particle.weight;
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
    };

    std::visit(functor, mc_unit->container);
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