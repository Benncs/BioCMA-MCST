#ifndef __CMA_UTILS_TRANSITIONNER_HPP__
#define __CMA_UTILS_TRANSITIONNER_HPP__

#include "cma_read/reactorstate.hpp"
#include "cma_utils/cache_hydro_state.hpp"
#include "cma_utils/iteration_state.hpp"
#include <cma_read/flow_iterator.hpp>
#include <cma_read/neighbors.hpp>
#include <cstddef>
#include <memory>
namespace CmaUtils
{
  enum class FlowmapTransitionMethod : char
  {
    Discontinuous,
    InterpolationFO,
  };

  class FlowMapTransitionner
  {
  public:
    static NeighborsView<HostSpace>
    get_neighbors_view(const CmaRead::Neighbors::Neighbors_const_view_t& liquid_neighbors);

    FlowMapTransitionner() = default;

    FlowMapTransitionner(std::size_t _n_flowmap,
                         std::size_t _n_per_flowmap,
                         std::size_t number_time_step,
                         std::unique_ptr<CmaRead::FlowIterator>&& _iterator,
                         bool is_two_phase_flow);

    virtual ~FlowMapTransitionner() = default;

    IterationState advance();

    IterationState advance_worker(std::span<double> flows,
                                  std::span<double> volumeLiq,
                                  std::span<double> volumeGas,
                                  const CmaRead::Neighbors::Neighbors_const_view_t& neighbors);

    [[nodiscard]] bool need_liquid_state() const noexcept;

    [[nodiscard]] bool is_two_phase_flow() const noexcept;

    [[nodiscard]] std::size_t get_n_timestep() const noexcept;

    [[nodiscard]] virtual const CmaRead::ReactorState&
    get_current_reactor_state() const noexcept = 0;

  protected:
    IterationState common_advance(NeighborsView<HostSpace> host_view);
    virtual void update_flow() = 0;
    void update_flow_worker(std::span<double> flows,
                            std::span<double> volumeLiq,
                            std::span<double> volumeGas,
                            const CmaRead::Neighbors::Neighbors_const_view_t& neighbors);
    virtual CmaUtils::PreCalculatedHydroState& current_liq_hydro_state() = 0;
    virtual CmaUtils::PreCalculatedHydroState& current_gas_hydro_state() = 0;

    void calculate_full_state(const CmaRead::ReactorState& reactor_state,
                              CmaUtils::PreCalculatedHydroState& liq_hydro_state,
                              CmaUtils::PreCalculatedHydroState& gas_hydro_state) const;

    [[nodiscard]] size_t getFlowIndex() const noexcept;

    [[nodiscard]] size_t size() const noexcept;

    [[nodiscard]] const CmaRead::ReactorState& get_unchecked(size_t index) const noexcept;

    std::vector<CmaUtils::PreCalculatedHydroState> liquid_pc;
    std::vector<CmaUtils::PreCalculatedHydroState> gas_pc;

  private:
    void update_counters();
    bool two_phase_flow;
    std::size_t n_per_flowmap;
    std::size_t n_flowmap;
    std::size_t n_timestep;
    std::size_t repetition_count;
    std::size_t current_flowmap_count;
    std::unique_ptr<CmaRead::FlowIterator> iterator;
  };

  [[nodiscard]] inline size_t FlowMapTransitionner::getFlowIndex() const noexcept
  {
    return this->repetition_count % this->n_flowmap;
  }

  [[nodiscard]] inline size_t FlowMapTransitionner::size() const noexcept
  {
    return iterator->size();
  }

  [[nodiscard]] inline bool FlowMapTransitionner::is_two_phase_flow() const noexcept
  {
    return two_phase_flow;
  }

  [[nodiscard]] inline bool FlowMapTransitionner::need_liquid_state() const noexcept
  {
    return repetition_count < n_flowmap && this->current_flowmap_count == 0;
  }

  [[nodiscard]] inline const CmaRead::ReactorState&
  FlowMapTransitionner::get_unchecked(const std::size_t index) const noexcept
  {
    return iterator->get_unchecked(index);
  };

}; // namespace CmaUtils

#endif