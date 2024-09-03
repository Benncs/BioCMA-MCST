#ifndef __SIMULATIONS_UNIT_HPP__
#define __SIMULATIONS_UNIT_HPP__

#include <Kokkos_Core.hpp>
#include <cassert>
#include <cma_read/reactorstate.hpp>
#include <common/common.hpp>
#include <common/kokkos_vector.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>
#include <mc/events.hpp>
#include <mc/particles/extra_process.hpp>
#include <mc/prng/prng.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <simulation/scalar_initializer.hpp>
#include <simulation/simulation_kernel.hpp>

static constexpr size_t trigger_const_particle_number = 1e6;
static constexpr bool const_number_simulation = false;
namespace Simulation
{

  class PreCalculatedHydroState;
  class ScalarSimulation;

  class SimulationUnit
  {
  public:
    SimulationUnit(
                   std::unique_ptr<MC::MonteCarloUnit> &&_unit,
                   ScalarInitializer scalar_init);

    ~SimulationUnit() = default;

    SimulationUnit(SimulationUnit &&other) noexcept;
    SimulationUnit(const SimulationUnit &other) = delete;
    SimulationUnit &operator=(SimulationUnit &&rhs) = delete;
    SimulationUnit &operator=(const SimulationUnit &rhs) = delete;

    std::unique_ptr<MC::MonteCarloUnit> mc_unit;

    [[nodiscard]] std::span<double> getCliqData() const;
    [[nodiscard]] std::tuple<size_t, size_t> getDim() const noexcept;
    [[nodiscard]] std::span<double> getCgasData() const;
    [[nodiscard]] std::span<double> getContributionData() const;

    void setVolumes(std::span<const double> volumesgas,
                    std::span<const double> volumesliq) const;

    void step(double d_t, const CmaRead::ReactorState &state) const;

    void cycleProcess(auto &&container, auto &&rview, double d_t);

    void setLiquidFlow(PreCalculatedHydroState *_flows_l);

    void setGasFlow(PreCalculatedHydroState *_flows_g);

    void reduceContribs(std::span<double> data, size_t n_rank) const;

    void clearContribution() const noexcept;

    void update_feed(double d_t) const;

    void clear_mc();

    void reset();

  private:
    bool const_number_simulation = true;

    /*Following definition are related to kokkos specifc view, those types are
     * only used during cycleprocess kernel, associated functions (get/set)
     * handle data transfer if necessary */
    using DiagonalViewCompute =
        Kokkos::View<double *,
                     Kokkos::LayoutLeft,
                     ComputeSpace,
                     Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

    using CumulativeProbabilityViewCompute =
        Kokkos::View<double **,
                     Kokkos::LayoutLeft,
                     ComputeSpace,
                     Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

    using ContributionViewCompute =
        Kokkos::View<double **, Kokkos::LayoutLeft, ComputeSpace>;

    using NeighborsViewCompute =
        Kokkos::View<const size_t **, Kokkos::LayoutStride, ComputeSpace>;

    using LeavingFlowIndexType = Kokkos::View<size_t *, ComputeSpace>;
    using LeavingFlowType = Kokkos::View<double *, ComputeSpace>;

    DiagonalViewCompute get_kernel_diagonal();

    CumulativeProbabilityViewCompute get_kernel_cumulative_proba();

    ContributionViewCompute get_kernel_contribution();

    void set_kernel_contribs_to_host(ContributionViewCompute c);

    [[nodiscard]] NeighborsViewCompute get_kernel_neighbors() const;

    LeavingFlowIndexType index_leaving_flow;
    LeavingFlowType leaving_flow;

    struct pimpl_deleter
    {
      void operator()(ScalarSimulation *) const;
    };

    using pimp_ptr_t = std::unique_ptr<ScalarSimulation, pimpl_deleter>;

    void post_init_compartments();

    void execute_process_knrl(const auto &kernel);

    bool is_two_phase_flow;

    PreCalculatedHydroState *flow_liquid; // TODO OPTI
    PreCalculatedHydroState *flow_gas;    // TODO OPTI

    pimp_ptr_t liquid_scalar;
    pimp_ptr_t gas_scalar;
    void post_init_concentration(init_scalar_f_t liquid, init_scalar_f_t gas);
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

  static inline void impl_cycle_process(double d_t,
                                        auto &&list,
                                        auto &&rview,
                                        auto &&local_compartments,
                                        auto &&neighbors,
                                        auto &&local_rng,
                                        auto &&_diag_transition,
                                        auto &&cumulative_probability,
                                        auto &&events,
                                        auto &&local_index_leaving_flow,
                                        auto &&local_leaving_flow,
                                        auto &&contribs) noexcept
  {

    Kokkos::parallel_for(
        "process", list.size(), KOKKOS_LAMBDA(const int i_particle) {
          auto &particle = list._owned_data(i_particle);
          auto &status = particle.properties.status;
          if (status == MC::CellStatus::DEAD)
          {
            return;
          }

          const size_t i_compartment = particle.properties.current_container;

          KernelInline::handle_move(i_compartment,
                                    particle,
                                    local_compartments,
                                    neighbors,
                                    local_rng,
                                    _diag_transition,
                                    cumulative_probability,
                                    events,
                                    d_t);

          const double random_number = local_rng.double_uniform();

          for (size_t i = 0; i < local_index_leaving_flow.size(); ++i)
          {
            const auto &index = local_index_leaving_flow(i);
            const auto &flow = local_leaving_flow(i);
            KernelInline::handle_exit(
                random_number, d_t, particle, local_compartments, index, flow);
          }

          if (status == MC::CellStatus::OUT)
          {
            events.template incr<MC::EventType::Exit>();
            status = MC::CellStatus::DEAD;
            return;
          }

          particle.update(d_t,
                          local_compartments(i_compartment).concentrations,
                          list.rng_instance);
          particle.contribution(contribs);

          if (status == MC::CellStatus::CYTOKINESIS)
          {
            events.template incr<MC::EventType::NewParticle>();

            const auto new_particle = particle.division();
            if constexpr (!const_number_simulation)
            {
              auto *np = rview().extra_process.spawn();
              if (np != nullptr)
              {
                *np = new_particle;
              }
              else
              {
                Kokkos::printf("SPAWNING OVERFLOW\r\n");
                Kokkos::atomic_increment(&rview().waiting_allocation_particle);
              }
            }
            else
            {
              particle.properties.weight *= 2;
              Kokkos::atomic_increment(
                  &local_compartments(particle.properties.current_container)
                       .n_cells);
            }

            assert(particle.properties.status!=MC::CellStatus::CYTOKINESIS);
            assert(new_particle.properties.status!=MC::CellStatus::CYTOKINESIS);
          }
          // assert(particle.data)
        });
    Kokkos::fence();
  }

  void SimulationUnit::cycleProcess(auto &&container, auto &&rview, double d_t)
  {
    auto &list = container.get_compute();

    const_number_simulation = (list.size() > trigger_const_particle_number);
    const auto diag_transition = get_kernel_diagonal();
    auto &local_compartments = mc_unit->domain.data();
    const LeavingFlowType &local_leaving_flow = leaving_flow;
    const LeavingFlowIndexType &local_index_leaving_flow = index_leaving_flow;

    const auto neighbors = get_kernel_neighbors();
    const auto cumulative_probability = get_kernel_cumulative_proba();

    auto &local_rng = mc_unit->rng;
    auto events = mc_unit->events;
    auto contribs = get_kernel_contribution();
    impl_cycle_process(d_t,
                       list,
                       rview,
                       local_compartments,
                       neighbors,
                       local_rng,
                       diag_transition,
                       cumulative_probability,
                       events,
                       local_index_leaving_flow,
                       local_leaving_flow,
                       contribs);

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
    list.insert(rview().extra_process);
    const auto n_new_alloc = rview().waiting_allocation_particle;
    const double new_weight =
        list._owned_data(0).properties.weight; // Weight is random, try to
                                               // find other initialisation
    list._spawn_alloc(n_new_alloc, new_weight);

    Kokkos::parallel_for(
        "add_new_alloc", n_new_alloc, KOKKOS_LAMBDA(const int ) {
          Kokkos::atomic_increment(&local_compartments(0).n_cells);
        });
    Kokkos::fence();

    set_kernel_contribs_to_host(contribs);
  }

} // namespace Simulation

#endif //__SIMULATIONS_UNIT_HPP__
