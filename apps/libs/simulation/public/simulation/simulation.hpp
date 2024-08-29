#ifndef __SIMULATIONS_UNIT_HPP__
#define __SIMULATIONS_UNIT_HPP__

#include "common/kokkos_vector.hpp"
#include "mc/events.hpp"
#include "mc/particles/data_holder.hpp"
#include "mc/particles/extra_process.hpp"
#include "mc/prng/prng.hpp"
#include <Kokkos_Atomic.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Printf.hpp>
#include <cma_read/reactorstate.hpp>
#include <common/common.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>
#include <mc/unit.hpp>
#include <memory>

#include <simulation/simulation_kernel.hpp>

constexpr size_t trigger_const_particle_number = 1e6;
namespace Simulation
{

  using init_scalar_f_t = void (*)(size_t, CmaRead::L2DView<double> &);
  static constexpr init_scalar_f_t default_gas_init =
      [](size_t, CmaRead::L2DView<double> &) {};

  struct ScalarInitializer
  {
    std::span<double> volumesgas;
    std::span<double> volumesliq;
    init_scalar_f_t liquid_f_init;
    init_scalar_f_t gaz_f_init = default_gas_init;
  };

  class PreCalculatedHydroState;
  class ScalarSimulation;

  class SimulationUnit
  {
  public:
    explicit SimulationUnit(const ExecInfo &info,
                            std::unique_ptr<MC::MonteCarloUnit> &&_unit,
                            std::span<double> volumesgas,
                            std::span<double> volumesliq,
                            size_t n_species,
                            bool _gas_flow = false);

    ~SimulationUnit() = default;

    SimulationUnit(SimulationUnit &&other) noexcept;
    SimulationUnit(const SimulationUnit &other) = delete;
    SimulationUnit &operator=(SimulationUnit &&rhs) = delete;
    SimulationUnit &operator=(const SimulationUnit &rhs) = delete;

    std::unique_ptr<MC::MonteCarloUnit> mc_unit;

    [[nodiscard]] std::span<double> getCliqData() const;
    [[nodiscard]] std::tuple<size_t, size_t> getDim() const;
    [[nodiscard]] std::span<double> getCgasData() const;
    [[nodiscard]] std::span<double> getContributionData() const;

    void setVolumes(std::span<const double> volumesgas,
                    std::span<const double> volumesliq) const;

    void step(double d_t, const CmaRead::ReactorState &state) const;

    void _cycleProces(auto &&container, auto &&rview, double d_t)
    {
      auto &list = container.get_compute();

      const_number_simulation = (list.size() > trigger_const_particle_number);
      const auto diag_transition = get_dg();
      auto &local_compartments = mc_unit->domain.data();
      const auto &local_leaving_flow = leaving_flow;
      const auto &local_index_leaving_flow = index_leaving_flow;

      const auto neighbors = get_view_neighbor();
      const auto cumulative_probability = get_cp();

      auto &local_rng = mc_unit->rng;
      auto events = mc_unit->events;
      auto contribs = get_contribs();
      auto& particle_rng = mc_unit->particle_rng;
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
                         contribs,particle_rng);

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
          "add_new_alloc", n_new_alloc, KOKKOS_LAMBDA(const int i) {
            Kokkos::atomic_increment(&local_compartments(0).n_cells);
          });
      Kokkos::fence();

      set_contribs(contribs);
    }

    void setLiquidFlow(PreCalculatedHydroState *_flows_l);

    void setGasFlow(PreCalculatedHydroState *_flows_g);

    void reduceContribs(std::span<double> data, size_t n_rank) const;

    void clearContribution() const;

    void update_feed(double d_t) const;

    void clear_mc();

    void reset()
    {
      liquid_scalar.reset();
      gas_scalar.reset();
      flow_liquid = nullptr;
      flow_gas = nullptr;
    }

  private:
    bool const_number_simulation = true;
    inline void impl_cycle_process(double d_t,
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
                                   auto &&contribs,auto&& particle_rng)
    {

      // Kokkos::Random_XorShift64_Pool<> p_rng(2024);

      auto local_const_number_simulation = const_number_simulation;
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

            const double random_number = local_rng.double_unfiform();

            for (size_t i = 0; i < local_index_leaving_flow.size(); ++i)
            {
              const auto &index = local_index_leaving_flow(i);
              const auto &flow = local_leaving_flow(i);
              KernelInline::handle_exit(random_number,
                                        d_t,
                                        particle,
                                        local_compartments,
                                        index,
                                        flow);
            }

            if (status == MC::CellStatus::OUT)
            {
              events.template incr<MC::EventType::Exit>();
              status = MC::CellStatus::DEAD;
              return;
            }

            particle.update(
                d_t, local_compartments(i_compartment).concentrations, particle_rng);
            particle.contribution(contribs);

            if (status == MC::CellStatus::CYTOKINESIS)
            {
              events.template incr<MC::EventType::NewParticle>();

              const auto new_particle = particle.division();
              // if (!local_const_number_simulation)
              // {
                auto *np = rview().extra_process.spawn();
                if (np != nullptr)
                {
                  *np = new_particle;
                }
                else
                {
                  Kokkos::printf("SPAWNING OVERFLOW\r\n");
                  Kokkos::atomic_increment(
                      &rview().waiting_allocation_particle);
                }
              // }
              // else
              // {
              //   particle.properties.weight *= 2;
              //   Kokkos::atomic_increment(
              //       &local_compartments(particle.properties.current_container)
              //            .n_cells);
              // }
            }
          });
      Kokkos::fence();
    }

    Kokkos::View<double *,
                 Kokkos::LayoutLeft,
                 ComputeSpace,
                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>
    get_dg();

    Kokkos::View<double **,
                 Kokkos::LayoutLeft,
                 ComputeSpace,
                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>
    get_cp();

    Kokkos::View<double **, Kokkos::LayoutLeft, ComputeSpace> get_contribs();

    void
    set_contribs(Kokkos::View<double **, Kokkos::LayoutLeft, ComputeSpace> c);

    [[nodiscard]] Kokkos::
        View<const size_t **, Kokkos::LayoutStride, ComputeSpace>
        get_view_neighbor() const;

    Kokkos::View<double *,
                 Kokkos::LayoutLeft,
                 ComputeSpace,
                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>
        diag_transition;

    Kokkos::View<size_t *, ComputeSpace> index_leaving_flow;
    Kokkos::View<double *, ComputeSpace> leaving_flow;

    struct pimpl_deleter
    {
      void operator()(ScalarSimulation *) const;
    };

    using pimp_ptr_t = std::unique_ptr<ScalarSimulation, pimpl_deleter>;

    void post_init_compartments();

    void execute_process_knrl(const auto &kernel);

    bool is_two_phase_flow;
    size_t n_thread;

    PreCalculatedHydroState *flow_liquid; // TODO OPTI
    PreCalculatedHydroState *flow_gas;    // TODO OPTI

    pimp_ptr_t liquid_scalar;
    pimp_ptr_t gas_scalar;
    void post_init_concentration();
  };

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

} // namespace Simulation

#endif //__SIMULATIONS_UNIT_HPP__
