#ifndef __SIMULATIONS_UNIT_HPP__
#define __SIMULATIONS_UNIT_HPP__

#include <simulation/feed_descriptor.hpp>
#include <simulation/probe.hpp>
#include <Kokkos_Assert.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <biocma_cst_config.hpp>
#include <cassert>
#include <cma_read/reactorstate.hpp>
#include <common/common.hpp>
#include <common/kokkos_vector.hpp>
#include <cstdint>
#include <mc/events.hpp>
#include <mc/particles/extra_process.hpp>
#include <mc/prng/prng.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <optional>
#include <simulation/alias.hpp>
#include <simulation/scalar_initializer.hpp>
#include <simulation/simulation_kernel.hpp>

#include <Kokkos_ScatterView.hpp>

// TODO Clean
static constexpr size_t trigger_const_particle_number = 1e6;

/**
 * @namespace Simulation
 * @brief Namespace that contains classes and structures related to simulation
 * handling
 */
namespace Simulation
{

  class PreCalculatedHydroState;
  class ScalarSimulation;

  class SimulationUnit
  {
  public:
    SimulationUnit(std::unique_ptr<MC::MonteCarloUnit> &&_unit,
                   const ScalarInitializer &scalar_init,std::optional<Feed::SimulationFeed> _feed=std::nullopt);

    ~SimulationUnit() = default;

    SimulationUnit(SimulationUnit &&other) noexcept;
    SimulationUnit(const SimulationUnit &other) = delete;
    SimulationUnit &operator=(SimulationUnit &&rhs) = delete;
    SimulationUnit &operator=(const SimulationUnit &rhs) = delete;

    std::unique_ptr<MC::MonteCarloUnit> mc_unit;

    [[nodiscard]] std::span<double> getCliqData() const;
    [[nodiscard]] std::tuple<size_t, size_t> getDim() const noexcept;

    [[nodiscard]] std::optional<std::span<double>> getCgasData() const;
    [[nodiscard]] std::span<double> getContributionData() const;

    void setVolumes(std::span<const double> volumesgas,
                    std::span<const double> volumesliq) const;

    void step(double d_t, const CmaRead::ReactorState &state) const;

    void cycleProcess(auto &&container, auto &&rview, double d_t);

    void setLiquidFlow(PreCalculatedHydroState *_flows_l);

    void setGasFlow(PreCalculatedHydroState *_flows_g);

    void reduceContribs(std::span<double> data, size_t n_rank) const;

    void clearContribution() const noexcept;

    void update_feed(double t, double d_t,bool update_scalar=true)noexcept;

    

    void clear_mc();

    void reset();

    [[nodiscard]] bool two_phase_flow()const{return is_two_phase_flow;}

    // FIXME
    Probes probes;

    void set_probes(Probes &&_probes)
    {
      probes = std::move(_probes);
    }

    [[nodiscard]] auto counter() const
    {
      return _internal_counter;
    }

  private:
    Simulation::Feed::SimulationFeed feed;

    bool const_number_simulation = true;
    uint64_t _internal_counter = 0;
    LeavingFlowIndexType index_leaving_flow;
    LeavingFlowType leaving_flow;
    PreCalculatedHydroState *flow_liquid; // TODO OPTI
    PreCalculatedHydroState *flow_gas;    // TODO OPTI
    bool is_two_phase_flow;

    DiagonalViewCompute get_kernel_diagonal();

    CumulativeProbabilityViewCompute get_kernel_cumulative_proba();

    cv get_kernel_contribution();

    void set_kernel_contribs_to_host(cv c);

    [[nodiscard]] NeighborsViewCompute get_kernel_neighbors() const;

    void post_init_compartments();
    void post_init_concentration(const ScalarInitializer &scalar_init);

    void post_init_concentration_functor(const ScalarInitializer &scalar_init);
    void post_init_concentration_file(const ScalarInitializer &scalar_init);

    template <class ListType,
              class ResultViewType,
              class CompartmentListType,
              typename ContribType>
    inline void post_kernel_process(ListType &list,
                                    ResultViewType &rview,
                                    CompartmentListType &local_compartments,
                                    ContribType &contribs);

    struct pimpl_deleter
    {
      void operator()(ScalarSimulation *) const;
    };
    using pimp_ptr_t = std::unique_ptr<ScalarSimulation, pimpl_deleter>;
    pimp_ptr_t liquid_scalar;
    pimp_ptr_t gas_scalar;
  };

  inline void SimulationUnit::reset()
  {
    liquid_scalar.reset();
    gas_scalar.reset();
    flow_liquid = nullptr;
    flow_gas = nullptr;
  }

  inline void SimulationUnit::setLiquidFlow(PreCalculatedHydroState *_flows_l)
  {
    flow_liquid = _flows_l;
  }

  inline void SimulationUnit::setGasFlow(PreCalculatedHydroState *_flows_g)
  {
    flow_gas = _flows_g;
  }

  inline void SimulationUnit::clear_mc()
  {
    mc_unit.reset();
  }

  void SimulationUnit::cycleProcess(auto &&container, auto &&rview, double d_t)
  {
    PROFILE_SECTION("cycleProcess")
    auto &list = container.get_compute();
    const size_t n_particle = list.size();

    const auto diag_transition = get_kernel_diagonal();
    const auto &local_leaving_flow = leaving_flow;
    const auto &local_index_leaving_flow = index_leaving_flow;
    const auto neighbors = get_kernel_neighbors();
    const auto cumulative_probability = get_kernel_cumulative_proba();

    auto &local_compartments = mc_unit->domain.data();
    auto &local_rng = mc_unit->rng;
    auto events = mc_unit->events;
    auto contribs = get_kernel_contribution();

    ContributionView contribs_scatter(contribs);

    const_number_simulation = (n_particle > trigger_const_particle_number);
    Kokkos::View<size_t, Kokkos::SharedSpace> internal_counter_dead(
        "internal_counter");

    Kokkos::deep_copy(internal_counter_dead, _internal_counter);
    auto k = KernelInline::Kernel(d_t,
                                  list,
                                  rview,
                                  internal_counter_dead,
                                  local_compartments,
                                  neighbors,
                                  local_rng,
                                  diag_transition,
                                  cumulative_probability,
                                  events,
                                  contribs_scatter,
                                  local_leaving_flow,
                                  local_index_leaving_flow,
                                  probes);
    Kokkos::parallel_for(
        "mc_cycle_process", Kokkos::RangePolicy<>(0, n_particle), k);
    Kokkos::fence("fence_mc_cycle_process");

    Kokkos::Experimental::contribute(contribs, contribs_scatter);

    _internal_counter = internal_counter_dead();
    post_kernel_process(list, rview, local_compartments, contribs);
  }

  template <class ListType,
            class ResultViewType,
            class CompartmentListType,
            typename ContribType>
  inline void
  SimulationUnit::post_kernel_process(ListType &list,
                                      ResultViewType &rview,
                                      CompartmentListType &local_compartments,
                                      ContribType &contribs)
  {
    Kokkos::parallel_for(
        "update_compartment_number",
        rview().extra_process.size(),
        KOKKOS_LAMBDA(const int i) {
          Kokkos::atomic_increment(
              &local_compartments(rview()
                                      .extra_process._owned_data(i)
                                      .properties.current_container)
                   .n_cells);
        });

    static constexpr uint64_t minimum_dead_particle_removal = 100;
    const auto threshold = std::max(
        minimum_dead_particle_removal,
        static_cast<uint64_t>(static_cast<double>(list.size()) *
                              AutoGenerated::dead_particle_ratio_threshold));

    if (_internal_counter > threshold)
    {
#ifndef NDEBUG
      const auto old_size = list.size();
#endif
      list.remove_dead(_internal_counter);
#ifndef NDEBUG
      KOKKOS_ASSERT(list.size() == old_size - _internal_counter);
#endif

      _internal_counter = 0;
    }

    list.insert(rview().extra_process);
    const auto n_new_alloc = rview().waiting_allocation_particle;
    const double new_weight =
        list._owned_data(0).properties.weight; // Weight is random, try to
                                               // find other initialisation
    list._spawn_alloc(n_new_alloc, new_weight);

    Kokkos::parallel_for(
        "add_new_alloc", n_new_alloc, KOKKOS_LAMBDA(const int) {
          Kokkos::atomic_increment(&local_compartments(0).n_cells);
        });
    Kokkos::fence("Fence cycle process ");

    set_kernel_contribs_to_host(contribs);
  }

} // namespace Simulation

#endif //__SIMULATIONS_UNIT_HPP__
