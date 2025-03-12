#ifndef __CORE_GLOBAL_INITIALLISER_HPP__
#define __CORE_GLOBAL_INITIALLISER_HPP__

#include "transitionner/transitionner.hpp"
#include <core/scalar_factory.hpp>
#include <core/simulation_parameters.hpp>
#include <array>
#include <cma_read/flow_iterator.hpp>
#include <cma_read/neighbors.hpp>
#include <cmt_common/cma_case.hpp>
#include <common/execinfo.hpp>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <mc/unit.hpp>
#include <memory>
#include <optional>
#include <simulation/feed_descriptor.hpp>
#include <simulation/scalar_initializer.hpp>
#include <simulation/simulation.hpp>
#include <string>
#include <vector>


class ILoadBalancer;

namespace Core
{

  /**
   * @class GlobalInitialiser
   * @brief A class responsible for initializing various components of a simulation framework.
   *
   */
  class GlobalInitialiser
  {
  public:
    /**
     * @brief Type alias for an optional unique pointer.
     *
     * This template alias provides a way to represent an optional ownership of
     * a unique pointer of type T.
     */
    template <typename T> using OptionalPtr = std::optional<std::unique_ptr<T>>;

    /**
     * @brief Constructs a GlobalInitialiser instance.
     *
     * @param _info Execution information containing context for the simulation.
     * @param _params Simulation parameters that configure the simulation behavior.
     */
    GlobalInitialiser(const ExecInfo& _info, UserControlParameters _user_params);

    /**
     * @brief Initializes a flow iterator.
     *
     * @return An optional unique pointer to the initialized flow iterator.
     */
    OptionalPtr<CmaRead::FlowIterator> init_flow_iterator();

    /**
     * @brief Initializes the simulation state.
     *
     * @param flow_handle A unique pointer to the flow iterator used for state initialization.
     * @return An optional boolean indicating success or failure of the initialization.
     */
    std::optional<bool> init_state(std::unique_ptr<CmaRead::FlowIterator>& flow_handle);

    /**
     * @brief Initializes a transitioner with the provided flow iterator.
     *
     * @param flow_handle A unique pointer to the flow iterator.
     * @return An optional unique pointer to the initialized flow map transitioner.
     */
    OptionalPtr<CmaUtils::FlowMapTransitionner>
    init_transitionner(std::unique_ptr<CmaRead::FlowIterator>&& flow_handle);

    /**
     * @brief Initializes a transitioner without parameters.
     *
     * @return An optional unique pointer to the initialized flow map transitioner.
     */
    OptionalPtr<CmaUtils::FlowMapTransitionner> init_transitionner();

    /**
     * @brief Initializes a simulation unit.
     *
     * @return An optional unique pointer to the initialized simulation unit.
     */
    OptionalPtr<Simulation::SimulationUnit> init_simulation(std::optional<Core::ScalarFactory::ScalarVariant> variant=std::nullopt);


    /**
     * @brief Initializes a simulation unit with specified parameters.
     *
     * @param _unit A unique pointer to the Monte Carlo unit initialized.
     * @param scalar_init A scalar initializer for the simulation.
     * @param _feed An optional simulation feed for initialization.
     * @return An optional unique pointer to the initialized simulation unit.
     */
    OptionalPtr<Simulation::SimulationUnit>
    init_simulation(std::unique_ptr<MC::MonteCarloUnit> _unit,
                    const Simulation::ScalarInitializer& scalar_init);

    /**
     * @brief Initializes a scalar component of the simulation.
     *
     * @return An optional scalar initializer.
     */
    std::optional<Simulation::ScalarInitializer> init_scalar();

    std::optional<Simulation::ScalarInitializer> init_scalar(Core::ScalarFactory::ScalarVariant&& variant);

    /**
     * @brief Initializes a simulation feed.
     *
     * @return An optional simulation feed.
     */
    void init_feed(std::optional<Simulation::Feed::SimulationFeed> feed = std::nullopt);

    /**
     * @brief Initializes a Monte Carlo unit.
     *
     * @return An optional unique pointer to the initialized Monte Carlo unit.
     */
    OptionalPtr<MC::MonteCarloUnit> init_monte_carlo();

    /**
     * @brief Checks if all initialization steps have been validated.
     *
     * @return A boolean indicating whether initialization can be terminated.
     */
    [[nodiscard]] bool check_init_terminate() const;

    [[nodiscard]] SimulationParameters get_parameters() const;

    void set_initial_number_particle(uint64_t np) noexcept;

  private:
    /**
     * @brief Enum to define initialization steps.
     *
     * Each member represents a distinct step in the initialization process,
     * except for the last member, `Count`, which signifies the total number
     * of initialization steps.
     *
     * @note This enum class leverages a C++ feature where if the first member
     *       is initialized to 0, the subsequent members are automatically
     *       assigned incrementing integer values. Thus, the `Count` member
     *       effectively represents the total number of steps (n), minus one
     *       for zero-based indexing.
     */
    enum class InitStep : std::size_t
    {
      FlowIterator = 0, ///< Step for initializing the flow iterator.
      Transitioner,     ///< Step for initializing the transitioner.
      InitState,        ///< Step for initializing the simulation state.
      Scalar,           ///< Step for initializing scalar values.
      Feed,             ///< Step for initializing the simulation feed.
      MC,               ///< Step for initializing the Monte Carlo unit.
      SimulationUnit,   ///< Step for initializing the simulation unit.
      Count             ///< Total number of steps in the initialization sequence.
    };

    /**
     * @brief Array to track validated initialization steps.
     */
    std::array<bool, static_cast<size_t>(InitStep::Count)> validated_steps{};

    /**
     * @brief Checks the validity of specified initialization steps.
     *
     * @param step The current initialization step to check.
     * @param args Additional steps to check.
     * @return A boolean indicating whether all specified steps are validated.
     */
    template <typename... Args> [[nodiscard]] bool check_steps(InitStep step, Args... args) const
    {
      if (!validated_steps[static_cast<size_t>(step)]) // NOLINT
      {
        return false;
      }

      return (... && validated_steps[static_cast<size_t>(args)]);
    }

    /**
     * @brief Validates the specified initialization step and any additional steps.
     *
     * @param step The current initialization step to validate.
     * @param args Additional steps to validate.
     */
    template <typename... Args> void validate_step(InitStep step, Args... args)
    {
      validated_steps[static_cast<size_t>(step)] = true; // NOLINT
      (void)std::initializer_list<int>{(validated_steps[static_cast<size_t>(args)] = true, 0)...};
    }

    /**
     * @brief Retrieves path files based on the provided case path.
     *
     * @param cma_case_path The path for the CMA case.
     * @return A structure containing exported paths.
     */
    CmtCommons::cma_exported_paths_t get_path_files(const std::string& cma_case_path);

    /**
     * @brief Host speciffic state init
     */
    std::optional<bool> host_init_state(std::unique_ptr<CmaRead::FlowIterator>& flow_handle);

    /**
     * @brief Performs MPI broadcast for synchronization.
     * @note: this function compiles and doesn't throw error even if no_mpi is not defined
     */
    void mpi_broadcast(); ///< Method for handling MPI broadcast communication.



    ExecInfo info;
    SimulationParameters params;

    UserControlParameters user_params;

    /////INNER STRUCT
    uint64_t particle_per_process;
    std::vector<double> liquid_volume;
    std::vector<double> gas_volume;
    CmaRead::Neighbors::Neighbors_const_view_t liquid_neighbors;
    double t_per_flowmap{};
    std::vector<size_t> worker_neighbor_data;
    bool f_init_gas_flow;
    std::optional<Simulation::Feed::SimulationFeed> feed;
    /////

    bool is_host; ///< Flag indicating if this instance is the host.
  };

} // namespace Core

#endif //__CORE_GLOBAL_INITIALLISER_HPP__