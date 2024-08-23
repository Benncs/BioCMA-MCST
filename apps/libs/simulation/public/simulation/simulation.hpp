#ifndef __SIMULATIONS_UNIT_HPP__
#define __SIMULATIONS_UNIT_HPP__

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <cma_read/reactorstate.hpp>
#include <common/common.hpp>
#include <mc/unit.hpp>
#include <memory>

struct Kernel;

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

    void cycleProces(auto&& container,double d_t)
    {
      auto &list = container.get_compute();

      Kokkos::parallel_for(
          "process", list.size(), KOKKOS_LAMBDA(auto&& i){
            //Do something 
          });
      Kokkos::fence();

    }

    void cycleProcess(double d_t);

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
    Kokkos::View<MC::ContainerState *, Kokkos::LayoutRight> domain_view;
    // MC::Results kernel_results;
    struct pimpl_deleter
    {
      void operator()(ScalarSimulation *) const;
    };

    struct pimpl_deleter_
    {
      void operator()(Kernel *) const;
    };

    using pimp_ptr_t = std::unique_ptr<ScalarSimulation, pimpl_deleter>;

    void post_init_compartments();

    void post_process_reducing();
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
