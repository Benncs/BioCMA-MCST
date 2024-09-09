#include "common/kokkos_vector.hpp"
#include "mc/particles/particle_model.hpp"
#include "mc/particles/particles_container.hpp"
#include <Kokkos_Atomic.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_DynamicView.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Printf.hpp>
#include <common/execinfo.hpp>
#include <dataexporter/data_exporter.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>
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

  static void get_particle_properties_opti(
      Kokkos::View<std::string *, HostSpace> &host_names,
      Kokkos::View<double **, HostSpace> &host_particle_values,
      Kokkos::View<double **, HostSpace> &host_spatial_values,
      std::unique_ptr<MC::MonteCarloUnit> &mc_unit,
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

      Kokkos::View<std::string *, HostSpace> names;
      Kokkos::View<double **, HostSpace> particle_values;
      Kokkos::View<double **, HostSpace> spatial_values;

      get_particle_properties_opti(
          names, particle_values, spatial_values, simulation.mc_unit, false);
      std::string ds_name = "biological_model/initial/";
      exporter->write_particle_data(
          names, particle_values, spatial_values, ds_name);
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

      exporter->write_final_results(simulation, distribution);

      // Results depending on simulation are exported or copied into properties,
      // clear montecarlo state to save memory
      auto unit = std::move(simulation.mc_unit);
      std::cout << "EXPORTING PARTICLE DATA" << std::endl;

      Kokkos::View<std::string *, HostSpace> names;
      Kokkos::View<double **, HostSpace> particle_values;
      Kokkos::View<double **, HostSpace> spatial_values;

      get_particle_properties_opti(
          names, particle_values, spatial_values, unit, true);
      std::string ds_name = "biological_model/final/";
      exporter->write_particle_data(
          names, particle_values, spatial_values, ds_name);
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
      Kokkos::View<std::string *, HostSpace> names;
      Kokkos::View<double **, HostSpace> particle_values;
      Kokkos::View<double **, HostSpace> spatial_values;

      get_particle_properties_opti(
          names, particle_values, spatial_values, simulation.mc_unit, false);
      std::string ds_name = "biological_model/" + std::to_string(counter) + "/";
      exporter->write_particle_data(
          names, particle_values, spatial_values, ds_name);
    }
  }

  void
  user_triggered_properties_export(Simulation::SimulationUnit &sim,
                                   std::unique_ptr<DataExporter> &data_exporter)
  {
    static int counter = 1;
    counter++;
    std::string ds_name = "biological_model/" + std::to_string(counter) + "/";
    PostProcessing::append_properties(counter, sim, data_exporter);
    // std::ofstream os("./out.cereal", std::ios::binary);
    // cereal::BinaryOutputArchive archive(os);
    // archive(*sim.mc_unit);
    // std::cout << "./out.cereal   " << sim.mc_unit->domain[0].n_cells
    //           << std::endl;
  }

  void get_particle_properties_opti(
      Kokkos::View<std::string *, HostSpace> &host_names,
      Kokkos::View<double **, HostSpace> &host_particle_values,
      Kokkos::View<double **, HostSpace> &host_spatial_values,
      std::unique_ptr<MC::MonteCarloUnit> &mc_unit,
      bool clean)
  {

    auto get_names = [&host_names](auto &host_particles_data)
    {
      const auto first_property = host_particles_data[0].data.get_properties();
      const size_t n_keys = first_property.size();
      host_names = Kokkos::View<std::string *, HostSpace>(
          "host_particle_property_name", n_keys);

      size_t i = 0;
      for (auto &[k, v] : first_property)
      {
        host_names(i) = std::move(k);
        i++;
      }
    };

    const size_t n_compartment = mc_unit->domain.getNumberCompartments();
    const auto compartments = mc_unit->domain.data();
    auto visitor = [clean,
                    n_compartment,
                    compartments,
                    get_names,
                    &host_names,
                    &host_particle_values,
                    &host_spatial_values](auto &&container)
    {
      auto &compute_particles_data = container.get_compute();
      auto &host_particles_data = container.get_host();

      get_names(host_particles_data);

      auto property_names =
          Kokkos::create_mirror_view_and_copy(ComputeSpace(), host_names);

      size_t n_p = compute_particles_data.size();
      Kokkos::View<double **> particle_values(
          "host_property_values", property_names.size(), n_p);

      Kokkos::View<double **> spatial_values(
          "host_spatial_values", property_names.size(), n_compartment);

      Kokkos::parallel_for(
          "get_particle_properties", n_p, KOKKOS_LAMBDA(const int i) {
            auto &particle = compute_particles_data._owned_data(i);
            if (particle.properties.status == MC::CellStatus::IDLE)
            {
              auto prop = particle.data.get_properties();
              const size_t i_container = particle.properties.current_container;

              size_t i_key = 0;
              for (const auto &[key, value] : prop)
              {

                particle_values(i_key, i) = value;

                Kokkos::atomic_add(
                    &spatial_values(i_key, i_container),
                    value /
                        static_cast<double>(compartments[i_container].n_cells));
                i_key++;
              }
            }
            if (clean)
            {
              particle.clearState();
            }
          });

      host_particle_values =
          Kokkos::create_mirror_view_and_copy(HostSpace(), particle_values);
      host_spatial_values =
          Kokkos::create_mirror_view_and_copy(HostSpace(), spatial_values);
    };

    std::visit(visitor, mc_unit->container);
  }

} // namespace PostProcessing