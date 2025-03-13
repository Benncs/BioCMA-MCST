#ifndef __SIMULATIONS_UNIT_HPP__
#define __SIMULATIONS_UNIT_HPP__

#include "mc/domain.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>
#include <biocma_cst_config.hpp>
#include <cassert>
#include <cma_utils/iteration_state.hpp>
#include <common/common.hpp>
#include <common/kokkos_vector.hpp>
#include <cstddef>
#include <cstdint>
#include <mc/events.hpp>
#include <mc/prng/prng.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <optional>
#include <simulation/alias.hpp>
#include <simulation/feed_descriptor.hpp>
#include <simulation/mass_transfer.hpp>
#include <simulation/move_kernel.hpp>
#include <simulation/probe.hpp>
#include <simulation/scalar_initializer.hpp>
#include <simulation/simulation_kernel.hpp>
// TODO Clean
static constexpr size_t trigger_const_particle_number = 1e6;

namespace CmaUtils
{
  class PreCalculatedHydroState;
} // namespace CmaUtils

/**
 * @namespace Simulation
 * @brief Namespace that contains classes and structures related to simulation
 * handling
 */
namespace Simulation
{

  struct Dimensions
  {
    std::size_t n_species{};
    std::size_t n_compartment{};

    template <class Archive> void serialize(Archive& archive)
    {
      archive(n_species, n_compartment);
    }
  };

  class ScalarSimulation;

  class SimulationUnit
  {
  public:
    SimulationUnit(std::unique_ptr<MC::MonteCarloUnit>&& _unit,
                   const ScalarInitializer& scalar_init,
                   std::optional<Feed::SimulationFeed> _feed = std::nullopt);

    ~SimulationUnit();

    SimulationUnit(SimulationUnit&& other) noexcept;
    SimulationUnit(const SimulationUnit& other) = delete;
    SimulationUnit& operator=(SimulationUnit&& rhs) = delete;
    SimulationUnit& operator=(const SimulationUnit& rhs) = delete;

    std::unique_ptr<MC::MonteCarloUnit> mc_unit;

    [[nodiscard]] Kokkos::View<const double**,
                               Kokkos::LayoutLeft,
                               ComputeSpace,
                               Kokkos::MemoryTraits<Kokkos::RandomAccess>>
    getkernel_concentration() const;
    [[nodiscard]] std::span<double> getCliqData() const;
    [[nodiscard]] std::optional<std::span<const double>> getCgasData() const;
    [[nodiscard]] std::optional<std::span<const double>> getMTRData() const;
    [[nodiscard]] Dimensions getDimensions() const noexcept;
    [[nodiscard]] std::span<const double> getContributionData() const;
    std::span<double> getContributionData_mut();
    [[nodiscard]] const Simulation::Feed::SimulationFeed& get_feed() const;
    ProbeAutogeneratedBuffer& get_probes();

    void set_probes(ProbeAutogeneratedBuffer&& _probes);

    // Simulation methods
    void cycleProcess(auto&& container, double d_t);

    void step(double d_t) const;

    void reduceContribs(std::span<const double> data, size_t n_rank) const;

    [[deprecated("perf:not useful")]] void
    reduceContribs_per_rank(std::span<const double> data) const;

    void clearContribution() const noexcept;

    void update_feed(double t, double d_t, bool update_scalar = true) noexcept;

    [[nodiscard]] bool two_phase_flow() const;

    [[nodiscard]] std::size_t counter() const;

    [[nodiscard]] double& get_start_time_mut();

    [[nodiscard]] double& get_end_time_mut();

    void update(CmaUtils::IterationState&& new_state);

    const CmaUtils::IterationState& get_state() const;

    // template <class ListType>
    // inline void post_kernel_process(ListType& list,
    //                                 ListType& process_buffer,
    //                                 // MC::ReactorDomain::n_cells_view_type& n_cells,
    //                                 size_t _waiting_allocation_particle);
    // Memory management
    void clear_mc();

    void reset();
    void post_init_compartments();

  private:
    void setVolumes() const;
    Kokkos::TeamPolicy<ComputeSpace> _policy;
    bool set_policy = false;
    void setLiquidFlow(CmaUtils::PreCalculatedHydroState* _flows_l);
    void setGasFlow(CmaUtils::PreCalculatedHydroState* _flows_g);
    // Attributes
    std::size_t internal_counter_dead;
    ProbeAutogeneratedBuffer probes;
    Simulation::Feed::SimulationFeed feed;

    CmaUtils::IterationState state;
    bool const_number_simulation = true;
    const bool is_two_phase_flow;
    double starting_time = 0.; // Not used within calculation, only for export purposes
    double end_time{};         // Not used within calculation, only for export purposes

    // Bounce methods to pimpl
    [[nodiscard]] DiagonalView<ComputeSpace> get_kernel_diagonal() const;
    CumulativeProbabilityView<ComputeSpace> get_kernel_cumulative_proba();
    [[nodiscard]] kernelContribution get_kernel_contribution() const;

    void set_kernel_contribs_to_host();

    void post_init_concentration(const ScalarInitializer& scalar_init);
    void post_init_concentration_functor(const ScalarInitializer& scalar_init);
    void post_init_concentration_file(const ScalarInitializer& scalar_init);

    std::shared_ptr<ScalarSimulation> liquid_scalar;
    std::shared_ptr<ScalarSimulation> gas_scalar;
    MassTransfer::MassTransferModel mt_model;
    KernelInline::MoveInfo<ComputeSpace> move_info;
  };

  inline ProbeAutogeneratedBuffer& SimulationUnit::get_probes()
  {
    return probes;
  }

  inline const CmaUtils::IterationState& SimulationUnit::get_state() const
  {
    return state;
  }

  void SimulationUnit::cycleProcess(auto&& container, double d_t)
  {

    PROFILE_SECTION("cycleProcess")
    using CurrentModel = typename std::remove_reference<decltype(container)>::type::UsedModel;
    const size_t n_particle = container.n_particles();

    auto& local_rng = mc_unit->rng;
    auto events = mc_unit->events;
    auto contribs = get_kernel_contribution();

    MC::ContributionView contribs_scatter(contribs);

    this->move_info.cumulative_probability = get_kernel_cumulative_proba();
    this->move_info.diag_transition = get_kernel_diagonal();
    this->move_info.neighbors = mc_unit->domain.getNeighbors();

    int out_total = 0;
    int dead_total = 0;
    int _waiting_allocation_particle = 0;

    auto f = Simulation::KernelInline::CycleFunctor<CurrentModel>(
        d_t, container, local_rng.random_pool, getkernel_concentration(),contribs_scatter);


    if (!set_policy)
    {
      _policy = MC::get_policty(f, n_particle,true);
        set_policy = true;
    }

    // if (mc_unit->domain.getNumberCompartments() != 1)
    {
      Kokkos::parallel_reduce("cycle_move",
                              _policy,
                              Simulation::KernelInline::MoveFunctor(d_t,
                                                                    container.position,
                                                                    container.status,
                                                                    n_particle,
                                                                    move_info,
                                                                    local_rng.random_pool),
                              out_total);
    }

    Kokkos::parallel_reduce("cycle_model", _policy, f, _waiting_allocation_particle, dead_total);
    Kokkos::fence("fence_mc_cycle_process");

    internal_counter_dead += out_total;
    internal_counter_dead += dead_total;

    Kokkos::Experimental::contribute(contribs, contribs_scatter);

    static constexpr uint64_t minimum_dead_particle_removal = 100;
    const auto threshold =
        std::max(minimum_dead_particle_removal,
                 static_cast<uint64_t>(static_cast<double>(n_particle) *
                                       AutoGenerated::dead_particle_ratio_threshold));

    if (internal_counter_dead > threshold)
    {
      container.clean_dead(internal_counter_dead);
      internal_counter_dead = 0;
    }

    container.merge_buffer();

    if (_waiting_allocation_particle != 0)
    {
      // TODO add overflow
    }

    set_kernel_contribs_to_host();
  }

  //   template <class ListType>
  //   inline void SimulationUnit::post_kernel_process(ListType& list,
  //                                                   ListType& process_buffer,
  //                                                   // MC::ReactorDomain::n_cells_view_type&
  //                                                   n_cells, size_t _waiting_allocation_particle)
  //   {

  //     // Kokkos::parallel_for(
  //     //     "update_compartment_number", process_buffer.size(), KOKKOS_LAMBDA(const int i) {
  //     //       Kokkos::atomic_increment(
  //     //           &n_cells(process_buffer._owned_data(i).properties.current_container));
  //     //     });

  //     static constexpr uint64_t minimum_dead_particle_removal = 100;
  //     const auto threshold =
  //         std::max(minimum_dead_particle_removal,
  //                  static_cast<uint64_t>(static_cast<double>(list.size()) *
  //                                        AutoGenerated::dead_particle_ratio_threshold));

  //     // Use commented line to remove dead when there is no particle.
  //     // This lead to stop exporting state for next time steps. (id: n_step=50 but only 40
  //     biological
  //     // export)
  //     //  if (_internal_counter > threshold || _internal_counter==list.size())
  //     if (internal_counter_dead > threshold)
  //     {
  // #ifndef NDEBUG
  //       const auto old_size = list.size();
  // #endif
  //       list.remove_dead(internal_counter_dead);
  // #ifndef NDEBUG
  //       KOKKOS_ASSERT(list.size() == old_size - internal_counter_dead);
  // #endif
  //       internal_counter_dead = 0;
  //     }

  //     list.insert(process_buffer);
  //     process_buffer.clear();
  //     const size_t n_new_alloc = _waiting_allocation_particle;
  //     if (n_new_alloc != 0)
  //     {
  //       const double new_weight = this->mc_unit->init_weight;
  //       list._spawn_alloc(n_new_alloc, new_weight);

  //       // Kokkos::parallel_for(
  //       //     "add_new_alloc",
  //       //     Kokkos::RangePolicy<ComputeSpace>(0, n_new_alloc),
  //       //     KOKKOS_LAMBDA(const int) { Kokkos::atomic_increment(&n_cells(0)); });

  //       Kokkos::fence("Fence cycle process ");
  //     }

  //     set_kernel_contribs_to_host();
  //   }

} // namespace Simulation

#endif //__SIMULATIONS_UNIT_HPP__
