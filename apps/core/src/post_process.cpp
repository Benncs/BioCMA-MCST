#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>
#include <algorithm>
#include <common/execinfo.hpp>
#include <common/kokkos_vector.hpp>
#include <dataexporter/data_exporter.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>
#include <iostream>
#include <mc/unit.hpp>
#include <memory>
#include <numeric>
#include <post_process.hpp>
#include <simulation/simulation.hpp>
#include <variant>

namespace
{
  void append_properties(int counter,
                         Simulation::SimulationUnit& simulation,
                         Core::PartialExporter& partial_exporter);

  void get_particle_properties_opti(Kokkos::View<std::string*, HostSpace>& host_names,
                                    Kokkos::View<double**, HostSpace>& host_particle_values,
                                    Kokkos::View<double**, HostSpace>& host_spatial_values,
                                    std::unique_ptr<MC::MonteCarloUnit>& mc_unit,
                                    bool clean = true);
} // namespace

namespace PostProcessing
{

  void save_probes(Simulation::SimulationUnit& simulation, Core::PartialExporter& pde)
  {
    auto& probes = simulation.get_probes();

    //TODO: Find out if comment is necessary or not 
    // if (probes.need_export())
    {
      pde.write_probe(probes.get());
      probes.clear();
    }
  }

  // FIXME
  void save_final_particle_state(Simulation::SimulationUnit& simulation, Core::PartialExporter& pde)
  {
    // save_particle_sate(simulation, pde, "final", false);
    user_triggered_properties_export(simulation, pde);
  }

  void save_initial_particle_state(Simulation::SimulationUnit& simulation,
                                   Core::PartialExporter& pde)
  {
    user_triggered_properties_export(simulation, pde);
    // save_particle_sate(simulation, pde, "initial", false);
  }

  void final_post_processing(const ExecInfo& exec,
                             const Core::SimulationParameters& params,
                             Simulation::SimulationUnit& simulation,
                             std::unique_ptr<Core::MainExporter>& mde)
  {
    std::cout << "POST PROCESSING" << std::endl;
    auto distribution = simulation.mc_unit->domain.getRepartition();

    auto tot = std::accumulate(distribution.begin(), distribution.end(), static_cast<size_t>(0));
    auto removed = simulation.mc_unit->events.get<MC::EventType::Death>() +
                   simulation.mc_unit->events.get<MC::EventType::Exit>();
    auto new_p = simulation.mc_unit->events.get<MC::EventType::NewParticle>();

    // TODO Core::MainExporter& mde

    mde->write_final(simulation, distribution);
    // exporter->write_final_results(simulation, distribution);

    if (tot != (new_p - removed + params.number_particle))
    {
      std::cerr << ("Results are not coherent (Bad particle balance): ");
      std::cerr << tot << "=" << new_p << "-" << removed << "+" << params.number_particle
                << std::endl;
    }
  }

  void show_sumup_state(const Simulation::SimulationUnit& simulation) noexcept
  {
    // Assuming domain data is in sharedSpace
    for (auto&& c : simulation.mc_unit->domain)
    {
      Kokkos::printf("%d ", c.n_cells);
    }
    Kokkos::printf("\r\n");
  }

  void user_triggered_properties_export(Simulation::SimulationUnit& sim, Core::PartialExporter& pde)
  {
    static int counter = 0;
    std::string ds_name = "biological_model/" + std::to_string(counter) + "/";
    ::append_properties(counter, sim, pde);
    counter++;
  }

} // namespace PostProcessing

namespace
{
  void append_properties(int counter,
                         Simulation::SimulationUnit& simulation,
                         Core::PartialExporter& partial_exporter)
  {

    {
      Kokkos::View<std::string*, HostSpace> names;
      Kokkos::View<double**, HostSpace> particle_values;
      Kokkos::View<double**, HostSpace> spatial_values;

      get_particle_properties_opti(
          names, particle_values, spatial_values, simulation.mc_unit, false);
      std::string ds_name = "biological_model/" + std::to_string(counter) + "/";

      partial_exporter.write_particle_data(
          {names.data(), names.extent(0)}, particle_values, spatial_values, ds_name);
    }
  }

  void get_particle_properties_opti(Kokkos::View<std::string*, HostSpace>& host_names,
                                    Kokkos::View<double**, HostSpace>& host_particle_values,
                                    Kokkos::View<double**, HostSpace>& host_spatial_values,
                                    std::unique_ptr<MC::MonteCarloUnit>& mc_unit,
                                    bool clean)
  {
    PROFILE_SECTION("get_particle_properties_opti")
    auto get_names = [&host_names](auto& host_particles_data)
    {
      const auto first_property = host_particles_data[0].data.get_properties();
      const size_t n_keys = first_property.size() + 2; //+2 because of hydraulic and divison time
                                                       // in properties common for all model
      host_names = Kokkos::View<std::string*, HostSpace>("host_particle_property_name", n_keys);

      size_t i = 0;
      for (auto& [k, v] : first_property)
      {
        host_names(i) = std::move(k);
        i++;
      }
      host_names(i++) = "hydraulic_time";
      host_names(i) = "age";
    };

    const size_t n_compartment = mc_unit->domain.getNumberCompartments();
    const auto compartments = mc_unit->domain.data();
    auto visitor = [clean,
                    n_compartment,
                    compartments,
                    get_names,
                    &host_names,
                    &host_particle_values,
                    &host_spatial_values](auto&& container)
    {
      auto& compute_particles_data = container.get_compute();
      if (compute_particles_data.size() == 0)
      {
        return;
      }
      auto& host_particles_data = container.get_host();

      get_names(host_particles_data);

      auto property_names = Kokkos::create_mirror_view_and_copy(ComputeSpace(), host_names);

      size_t n_p = compute_particles_data.size();
      Kokkos::View<double**, ComputeSpace> particle_values(
          "host_property_values", property_names.size(), n_p);

      Kokkos::View<double**> _spatial_values(
          "host_spatial_values", property_names.size(), n_compartment);

      Kokkos::Experimental::ScatterView<double**> scatter_spatial_values(_spatial_values);
      Kokkos::parallel_for(
          "get_particle_properties", n_p, KOKKOS_LAMBDA(const int i) {
            auto& particle = compute_particles_data._owned_data(i);
            if (particle.properties.status == MC::CellStatus::IDLE)
            {
              auto spatial_values = scatter_spatial_values.access();
              auto prop = particle.data.get_properties();
              const size_t i_container = particle.properties.current_container;
              // const auto cast_n_cell =
              //     static_cast<double>(compartments[i_container].n_cells);
              size_t i_key = 0;
              for (const auto& [key, value] : prop)
              {

                particle_values(i_key, i) = value;

                // Kokkos::atomic_add(&spatial_values(i_key, i_container),
                //                    value / cast_n_cell);
                spatial_values(i_key, i_container) += value;
                ++i_key;
              }

              particle_values(i_key, i) = particle.properties.hydraulic_time;
              spatial_values(i_key, i_container) += particle.properties.hydraulic_time;
              // Kokkos::atomic_add(&spatial_values(i_key, i_container),
              //                    particle.properties.hydraulic_time /
              //                        cast_n_cell);
              ++i_key;
              particle_values(i_key, i) = particle.properties.interdivision_time;
              // Kokkos::atomic_add(&spatial_values(i_key, i_container),
              //                    particle.properties.interdivision_time /
              //                        cast_n_cell);
              spatial_values(i_key, i_container) += particle.properties.interdivision_time;
            }
            if (clean)
            {
              particle.clearState();
            }
          });
      Kokkos::fence();
      Kokkos::Experimental::contribute(_spatial_values, scatter_spatial_values);
      host_particle_values = Kokkos::create_mirror_view_and_copy(HostSpace(), particle_values);
      host_spatial_values = Kokkos::create_mirror_view_and_copy(HostSpace(), _spatial_values);
    };

    std::visit(visitor, mc_unit->container);
  }
} // namespace