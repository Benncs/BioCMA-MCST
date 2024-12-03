#ifndef __BIOMC_API_HPP__
#define __BIOMC_API_HPP__

#include <core/case_data.hpp>
#include <core/simulation_parameters.hpp>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <optional>
#include <simulation/feed_descriptor.hpp>
#include <span>
#include <string_view>

#include <api/results.hpp>
/**
 * @namespace Api
 * This namespace contains classes and functions related to the simulation API.
 */
namespace Api
{
  /**
   * @class SimulationInstance
   * @brief Represents an instance of a simulation with customizable parameters and behaviors.
   *
   * The SimulationInstance class is designed to encapsulate the configuration, execution,
   * and management of a simulation process. It supports initialization for single-threaded
   * or multi-threaded (optionally MPI-enabled) scenarios, parameter registration,
   * and execution control.
   *
   * The class is non-copyable but movable to ensure unique ownership and efficient resource
   * management.
   */
  class SimulationInstance
  {
  public:
    /**
     * @brief Deleted copy constructor to prevent copying.
     */
    SimulationInstance(const SimulationInstance&) = delete;
    /**
     * @brief Defaulted move constructor for efficient resource transfer.
     */
    SimulationInstance(SimulationInstance&&) = default;
    /**
     * @brief Deleted copy assignment operator to prevent copying.
     */
    SimulationInstance& operator=(const SimulationInstance&) = delete;
    /**
     * @brief Defaulted move assignment operator for efficient resource transfer.
     */
    SimulationInstance& operator=(SimulationInstance&&) = default;

    // TODO Enable if def USE_MPI
    /**
     * @brief Initialize a simulation instance with MPI support.
     *
     * @param n_rank The total number of ranks in the MPI group.
     * @param current_rank The rank ID for this instance.
     * @param id A unique identifier for the simulation instance.
     * @param thread_per_process Number of threads allocated per process.
     * @return An optional containing a unique pointer to the instance if successful,
     *         or std::nullopt if initialization failed.
     */
    static std::optional<std::unique_ptr<SimulationInstance>>
    init(uint32_t n_rank, uint32_t current_rank, uint64_t id, uint32_t thread_per_process) noexcept;

    /**
     * @brief Initialize a simulation instance without MPI support.
     *
     * @param id A unique identifier for the simulation instance.
     * @param thread_per_process Number of threads allocated per process.
     * @return An optional containing a unique pointer to the instance if successful,
     *         or std::nullopt if initialization failed.
     */
    static std::optional<std::unique_ptr<SimulationInstance>>
    init(uint64_t id, uint32_t thread_per_process) noexcept;

    /**
     * @brief Default constructor.
     */
    SimulationInstance() = default;
    /**
     * @brief Default destructor.
     */
    ~SimulationInstance() = default;

    /**
     * @brief Apply the simulation configuration and prepare for execution.
     *
     * @param to_load Flag indicating whether the configuration should be loaded.
     * @return An ApiResult indicating the success or failure of the operation.
     */
    ApiResult apply(bool to_load) noexcept;

    /**
     * @brief Apply the simulation configuration and prepare for execution.
     *
     * @return An ApiResult indicating the success or failure of the operation.
     */
    ApiResult apply() noexcept;

    /**
     * @brief Load ad apply the simulation configuration and prepare for execution based on file.
     *
     * @return An ApiResult indicating the success or failure of the operation.
     */
    ApiResult apply_load() noexcept;

    /**
     * @brief Register user control parameters for the simulation.
     *
     * @param params A rvalue reference to a UserControlParameters object.
     * @return An ApiResult indicating the success or failure of the operation.
     */
    ApiResult register_parameters(Core::UserControlParameters&& params) noexcept;

    /**
     * @brief Register a result output path.
     *
     * @param path The file system path to store simulation results.
     * @return True if the path was successfully registered; false otherwise.
     */
    bool register_result_path(std::string_view path);

    /**
     * @brief Register a path for CMA data.
     *
     * @param path The file system path for CMA data.
     * @param recursive Flag indicating if recursive registration is required.
     * @return True if the path was successfully registered; false otherwise.
     */
    bool register_cma_path(std::string_view path, bool recursive = false);

    /**
     * @brief Register a serialization/deserialization (serde) path.
     *
     * @param path The file system path for serde operations.
     * @return True if the path was successfully registered; false otherwise.
     */
    bool register_serde(std::string_view path);

    /**
     * @brief Register the model name for the simulation.
     *
     * @param path The model name or path.
     * @return An ApiResult indicating the success or failure of the operation.
     */
    ApiResult register_model_name(std::string_view path);

    /**
     * @brief Configure feed constants for the simulation.
     *
     * @param _f A constant feed value.
     * @param _target A span of target values.
     * @param _position A span of position indices.
     * @param _species A span of species indices.
     * @param gas Flag indicating whether the feed is gas-phase.
     * @return True if the feed was successfully configured; false otherwise.
     */
    bool set_feed_constant(double _f,
                           std::span<double> _target,
                           std::span<std::size_t> _position,
                           std::span<std::size_t> _species,
                           bool gas = false);



    bool set_feed_constant_from_rvalue(double _f,
                           std::vector<double>&& _target,
                           std::vector<std::size_t>&& _position,
                           std::vector<std::size_t>&& _species,
                           bool gas = false);

    /**
     * @brief Retrieve the simulation instance's unique identifier.
     *
     * @return The ID of the simulation instance.
     */
    [[nodiscard]] int get_id() const;

    /**
     * @brief Execute the simulation.
     *
     * @return An ApiResult indicating the success or failure of the operation.
     */
    ApiResult exec() noexcept;

  private:
    int id{}; ///< The unique identifier to connect with c api.

    /**
     * @brief Private constructor for internal initialization.
     *
     * @param n_rank The total number of ranks in the MPI group.
     * @param current_rank The rank ID for this instance.
     * @param id A unique identifier for the simulation instance.
     * @param thread_per_process Number of threads allocated per process.
     */
    SimulationInstance(uint32_t n_rank,
                       uint32_t current_rank,
                       uint64_t id,
                       uint32_t thread_per_process);

    Core::CaseData _data;               ///< Case data for the simulation.
    Core::UserControlParameters params; ///< User-defined control parameters.
    bool loaded = false;                ///< Flag indicating if the instance is loaded.
    bool applied = false;               ///< Flag indicating if the configuration is applied.
    bool registered = false;            ///< Flag indicating if resources are registered.
    std::optional<Simulation::Feed::SimulationFeed> feed =
        std::nullopt; ///< Optional feed configuration.
  };

} // namespace Api
#endif
