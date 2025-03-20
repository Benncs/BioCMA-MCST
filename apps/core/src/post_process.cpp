#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>
#include <common/execinfo.hpp>
#include <common/kokkos_vector.hpp>
#include <core/post_process.hpp>
#include <dataexporter/data_exporter.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>
#include <impl_post_process.hpp>
#include <iostream>
#include <mc/unit.hpp>
#include <memory>
#include <numeric>

#include <simulation/simulation.hpp>
#include <variant>

namespace
{
  void append_properties(int counter,
                         Simulation::SimulationUnit& simulation,
                         Core::PartialExporter& partial_exporter);
} // namespace

namespace PostProcessing
{

  void save_probes(Simulation::SimulationUnit& simulation, Core::PartialExporter& pde, bool force)
  {
    auto& probes = simulation.get_probes();

    // TODO: Find out if comment is necessary or not
    if (probes.need_export() || force)
    {
      pde.write_probe(probes.get()); // probe.get only returns the used chunk of memory id:
                                     // buffer_size if need export else internal counter
      probes.clear();
    }
  }
  static int counter = 0; // TODO Remove static and reset to 0 when new simulation. If handle is
                          // reused for two simulation as itś static counter is not reset

  // FIXME
  void save_particle_state(Simulation::SimulationUnit& simulation, Core::PartialExporter& pde)
  {
    ::append_properties(counter, simulation, pde);
    counter++;
  }
  void final_post_processing(const ExecInfo& exec,
                             const Core::SimulationParameters& params,
                             Simulation::SimulationUnit& simulation,
                             std::unique_ptr<Core::MainExporter>& mde)
  {
    std::cout << "POST PROCESSING" << std::endl;
    
    auto removed = simulation.mc_unit->events.get<MC::EventType::Death>() +
                   simulation.mc_unit->events.get<MC::EventType::Exit>();
    auto new_p = simulation.mc_unit->events.get<MC::EventType::NewParticle>();

    mde->write_final(simulation, simulation.mc_unit->n_particle());

    auto distribution = simulation.mc_unit->getRepartition();

    auto tot = std::accumulate(distribution.begin(), distribution.end(), static_cast<size_t>(0));
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
    for (auto&& c : simulation.mc_unit->getRepartition())
    {
      Kokkos::printf("%d ", c);
    }
    Kokkos::printf("\r\n");
  }

  void reset_counter()
  {
    counter = 0;
  }

} // namespace PostProcessing

namespace
{

  PostProcessing::BonceBuffer
  get_particle_properties_device(const std::unique_ptr<MC::MonteCarloUnit>& mc_unit);

  void append_properties(int counter,
                         Simulation::SimulationUnit& simulation,
                         Core::PartialExporter& partial_exporter)
  {

    {
      auto dump = get_particle_properties_device(simulation.mc_unit);

      std::string ds_name = "biological_model/" + std::to_string(counter) + "/";
      partial_exporter.write_particle_data(
          dump.vnames, dump.particle_values, dump.spatial_values, ds_name);
    }
  }

  PostProcessing::BonceBuffer
  get_particle_properties_device(const std::unique_ptr<MC::MonteCarloUnit>& mc_unit)
  {

    // BonceBuffer properties;
    const size_t n_compartment = mc_unit->domain.getNumberCompartments();

    // auto visitor = [&properties, n_compartment](auto& container) mutable
    // {
    //   auto list = container.get_compute();
    //   const std::size_t n_p = list.size(); // USE list size not Kokkos View size. ParticleList
    //                                        // allocates more particles than needed

    //   auto view = list._owned_data;
    //   using ListType = decltype(view);
    //   using ViewParticleType = ListType::value_type;
    //   using Model = ViewParticleType::Model;
    //   properties.vnames = Model::names();
    //   properties.vnames.emplace_back("hydraulic_time");
    //   properties.vnames.emplace_back("age");
    //   constexpr std::size_t n_properties = Model::get_number() + 2; //+2 for hydraulic time and
    //   age
    //   // TODO: Find out if we cna use compile time size or not
    //   // For the moment double** is used for simplicyt with particle model but algortihm remains
    //   // exactly the same if we know the size of view.
    //   // Kokkos::View<double* [n_properties], ComputeSpace> particle_values(
    //   //     "device_property_values", names.size(), n_p);

    //   ParticlePropertyViewType<ComputeSpace> particle_values("property_values", n_properties,
    //   n_p); ParticlePropertyViewType<ComputeSpace> spatial_values(
    //       "property_spatial", n_properties, n_compartment);

    //   inner<ListType, Kokkos::DefaultExecutionSpace>(n_p, view, particle_values, spatial_values);

    //   properties.particle_values =
    //       Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particle_values);
    //   properties.spatial_values =
    //       Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), spatial_values);
    // };

    // auto properties=std::visit(func, mc_unit->container);

    return std::visit(
        [n_compartment](auto& container)
        {
          auto list = container.get_compute();
          const std::size_t n_p = list.size(); // USE list size not Kokkos View size. ParticleList
                                               // allocates more particles than needed

          return PostProcessing::get_properties(n_p, list._owned_data, n_compartment);
        },
        mc_unit->container);

    // return std::visit([n_compartment](auto& container)
    //                   { PostProcessing::get_properties(container.get_compute(), n_compartment) },
    //                   mc_unit->container);
  }

} // namespace
