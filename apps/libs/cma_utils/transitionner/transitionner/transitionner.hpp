#ifndef __CMA_UTILS_TRANSITIONNER_HPP__
#define __CMA_UTILS_TRANSITIONNER_HPP__

#include <cma_read/reactorstate.hpp>
#include <cma_utils/cache_hydro_state.hpp>
#include <cma_utils/iteration_state.hpp>
#include <cma_read/flow_iterator.hpp>
#include <cma_read/neighbors.hpp>
#include <cstddef>
#include <memory>
#include <transitionner/proxy_cache.hpp>

namespace CmaUtils
{
  /**
   * @enum FlowmapTransitionMethod
   * @brief Defines the transition methods used between flowmaps.
   *   */
  enum class FlowmapTransitionMethod : char
  {
    Discontinuous,  ///< Represents a sharp, abrupt transition.
    InterpolationFO ///< Represents a smooth transition using first-order interpolation.
  };

  // Forward declaration
  class ProxyPreCalculatedHydroState;

  /**
   * @class FlowMapTransitionner
   * @brief Manages the reading, caching, and transitioning of flowmaps for simulation timesteps.
   *
   * This class reads flowmaps from a file and converts them into simulation states used at each
   * timestep. It handles transitions between consecutive flowmaps (from instant t to t + dt) and
   * supports periodic looping over the same flowmaps. To optimize performance, caching ensures that
   * flowmap-to-state conversions are not repeated.
   */
  class FlowMapTransitionner
  {
  public:
    FlowMapTransitionner() = default;
    FlowMapTransitionner& operator=(FlowMapTransitionner&&) = default;
    FlowMapTransitionner(FlowMapTransitionner&&) = default;

    FlowMapTransitionner& operator=(const FlowMapTransitionner&) = delete;
    FlowMapTransitionner(const FlowMapTransitionner&) = delete;

    FlowMapTransitionner(std::size_t _n_flowmap,
                         std::size_t _n_per_flowmap,
                         std::size_t number_time_step,
                         std::unique_ptr<CmaRead::FlowIterator>&& _iterator,
                         bool is_two_phase_flow);

    virtual ~FlowMapTransitionner() = default;

    /**
     * @brief Advances the simulation by one timestep.
     *
     * This method updates the current flowmap state, transitions to the next flowmap if needed,
     * and handles periodic looping.
     *
     * @return IterationState representing the current iteration's state.
     */
    IterationState advance();

    /**
     * @brief Advances the simulation by one timestep from a MPI Wroker.
     * @todo Use template to have conditional build
     * This method updates the current flowmap state, transitions to the next flowmap if needed,
     * and handles periodic looping.
     *
     * @return IterationState representing the current iteration's state.
     */
    IterationState advance_worker(std::span<double> flows,
                                  std::span<double> volumeLiq,
                                  std::span<double> volumeGas,
                                  const CmaRead::Neighbors::Neighbors_const_view_t& neighbors);

    /**
     * @brief Checks whether a gas flowmap is provided.
     *
     * This method returns true if a gas flowmap is available, indicating
     * a two-phase flow simulation.
     */
    [[nodiscard]] bool is_two_phase_flow() const noexcept;

    /**
     * @brief Returns the total number of expected simulation timesteps.
     */
    [[nodiscard]] std::size_t get_n_timestep() const noexcept;

    /**
     * @brief Provides a const reference to the current reactor state.
     *
     * This method returns the reactor state currently read from the flowmap.
     * @todo: Consider making this method protected or moving it to MPI_WRAP.
     */
    [[nodiscard]] virtual const CmaRead::ReactorState&
    get_current_reactor_state() const noexcept = 0;

    /**
     * @brief Converts the provided neighbors' data into a Kokkos view.
     *
     * This static method transforms neighbor information into a Kokkos view,
     * enabling efficient access and usage during the simulation.
     *
     * @param liquid_neighbors Constant view of the liquid neighbors.
     * @return Kokkos view of the neighbors.
     */
    static NeighborsView<HostSpace>
    get_neighbors_view(const CmaRead::Neighbors::Neighbors_const_view_t& liquid_neighbors);

  protected:
    /**
     * @brief Check if liquid state has to calculated or can be use via caching.
     */
    [[nodiscard]] bool need_liquid_state_calculation() const noexcept;

    /**
     * @brief Operation required by all methods in order to perform advance function .
     */
    IterationState common_advance(NeighborsView<HostSpace> host_view,
                                  std::unordered_map<std::string, std::span<const double>>&& info);

    /**
     * @brief Calculate the current state based on the selected method, may be caching or direct
     * calculation
     */
    virtual void update_flow() = 0;

    /**
     * @brief ONLY for MPI_WORKER Calculate the current state based on the selected method, may be
     * caching or direct calculation
     */
    void update_flow_worker(std::span<double> flows,
                            std::span<double> volumeLiq,
                            std::span<double> volumeGas,
                            const CmaRead::Neighbors::Neighbors_const_view_t& neighbors);

    /**
     * @brief Get current liquid state read from CmaRead
     */
    virtual CmaUtils::ProxyPreCalculatedHydroState& current_liq_hydro_state() = 0;
    /**
     * @brief Get current gas state read from CmaRead
     */
    virtual CmaUtils::ProxyPreCalculatedHydroState& current_gas_hydro_state() = 0;

    /**
     * @brief Calculate the whole state (liquid+gas+neighbor) for the current time step
     */
    void calculate_full_state(const CmaRead::ReactorState& reactor_state,
                              CmaUtils::ProxyPreCalculatedHydroState& liq_hydro_state,
                              CmaUtils::ProxyPreCalculatedHydroState& gas_hydro_state) const;

    /**
     * @brief Return the index of the current ProxyPreCalculatedHydroState
     */
    [[nodiscard]] size_t getFlowIndex() const noexcept;

    /**
     * @brief Get the number of load flowmap
     */
    [[nodiscard]] size_t size() const noexcept;

    /**
     * @brief Direct read state from iterator
     */
    [[nodiscard]] const CmaRead::ReactorState& get_unchecked(size_t index) const noexcept;

    // Buffer for caching, vector have the same size which is the number of flowmap read from case
    std::vector<CmaUtils::ProxyPreCalculatedHydroState> liquid_pc; ///< Buffer for liquid state
    std::vector<CmaUtils::ProxyPreCalculatedHydroState> gas_pc;    ///< Buffer for gas state

  private:
    /**
     * @brief Update flow and iteration counters
     */
    void update_counters();

    bool two_phase_flow{};          ///< Is two_phase_flow
    std::size_t n_per_flowmap{};    ///< Number of iteration per flowmap
    std::size_t n_flowmap{};        ///< Number of flowmap
    std::size_t n_timestep{};       ///< Total number of timestep
    std::size_t repetition_count{}; ///< How many flowmaps (different or not) have been used 
    std::size_t current_flowmap_count{}; ///< How many iteration have been performed since last
                                         ///< flowmap change (<n_per_flowmap)
    std::unique_ptr<CmaRead::FlowIterator> iterator; ///< Iterator that reads flowmap
  };

  [[nodiscard]] inline size_t FlowMapTransitionner::getFlowIndex() const noexcept
  {
    //Modulo is used to select correct index into buffers.
    //During the first loop we have repetition_count<n_flowmap but after looping we have repetition_count>n_flowmap 
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

  [[nodiscard]] inline bool FlowMapTransitionner::need_liquid_state_calculation() const noexcept
  {
    // We need to calculated new state if:
    //- it is the first time we have this flowmap (this->current_flowmap_count)
    //  We have not already looped through all the flowmap
    return repetition_count < n_flowmap && this->current_flowmap_count == 0;
  }

  [[nodiscard]] inline const CmaRead::ReactorState&
  FlowMapTransitionner::get_unchecked(const std::size_t index) const noexcept
  {
    return iterator->get_unchecked(index);
  };

}; // namespace CmaUtils

#endif