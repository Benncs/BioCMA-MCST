#ifndef __BIOMC_API_HPP__
#define __BIOMC_API_HPP__

#include "common/logger.hpp"
#include <api/results.hpp>
#include <common/execinfo.hpp>
#include <core/case_data.hpp>
#include <core/scalar_factory.hpp>
#include <core/simulation_parameters.hpp>
#include <cstdlib>
#include <memory>
#include <optional>
#include <simulation/feed_descriptor.hpp>
#include <simulation/mass_transfer.hpp>
#include <span>
#include <string_view>
#include <utility>

/**
 * @brief Api
 * This namespace contains classes and functions related to the simulation API.
 */
namespace Api
{

  void finalise();

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
   * Note: Define this struct as struct and not class to be coherent with the C opaque type defined
   * as typedef struct This change absolutely nothing to the code behaviour it's just a name
   */
  struct SimulationInstance
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

    /**
     * @brief Initialize a simulation instance
     *
     * @param argc Number of runtime argument.
     * @param argv Tuntime arguments.
     * @param id A unique identifier for the simulation instance (optional).
     * @return An optional containing a unique pointer to the instance if successful,
     *         or std::nullopt if initialization failed.
     */
    static std::optional<std::unique_ptr<SimulationInstance>>
    init(int argc, char** argv, std::optional<std::size_t> run_id = std::nullopt) noexcept;

    /**
     * @brief Default constructor.
     */
    SimulationInstance() = delete;
    /**
     * @brief Default destructor.
     */
    ~SimulationInstance();

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

    ApiResult register_initial_condition(Core::ScalarFactory::ScalarVariant&& type);
    ApiResult register_initialiser_file_path(std::string_view path);
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

    ApiResult set_feed(Simulation::Feed::FeedDescriptor feed_variant, Phase phase = Phase::Liquid);

    ApiResult set_feed_constant(double _flow,
                                double _concentration,
                                std::size_t _species,
                                std::size_t _position,
                                bool gas = false,
                                bool fed_batch = false);

    ApiResult set_feed_constant_different_output(double _flow,
                                                 double _concentration,
                                                 std::size_t _species,
                                                 std::size_t input_position,
                                                 std::size_t output_position,
                                                 bool gas = false);

    // /**
    //  * @brief Configure feed constants for the simulation.
    //  *
    //  * @param _flow A constant flow value.
    //  * @param concentrations A span of concentrations values.
    //  * @param _position A span of position indices.
    //  * @param _species A span of species indices.
    //  * @param gas Flag indicating whether the feed is gas-phase.
    //  * @return True if the feed was successfully configured; false otherwise.
    //  */
    // [[deprecated]]bool set_feed_constant(double _flow,
    //                        std::span<double> _concentration,
    //                        std::span<std::size_t> _position,
    //                        std::span<std::size_t> _species,
    //                        bool gas = false,
    //                        bool fed_batch = false);

    // [[deprecated]]bool set_feed_constant_from_position(double _flow,
    //                                      std::span<double> _concentration,
    //                                      std::span<std::size_t> _position,
    //                                      std::span<std::size_t> _species,
    //                                      std::optional<std::vector<std::size_t>>
    //                                      _output_position, bool gas = false, bool fed_batch =
    //                                      false);

    // [[deprecated]]bool set_feed_constant_from_rvalue(double _f,
    //                                    std::vector<double>&& _target,
    //                                    std::vector<std::size_t>&& _position,
    //                                    std::vector<std::size_t>&& _species,
    //                                    bool gas = false,
    //                                    bool fed_batch = false);

    ApiResult register_scalar_initiazer(Core::ScalarFactory::ScalarVariant&& var);

    /**
     * @brief Retrieve the simulation instance's unique identifier.
     *
     * @return The ID of the simulation instance.
     */
    [[nodiscard]] int get_id() const;

    [[nodiscard]] const ExecInfo& get_exec_info() const;

    /**
     * @brief Execute the simulation.
     *
     * @return An ApiResult indicating the success or failure of the operation.
     */
    ApiResult exec() noexcept;

    void set_logger(std::shared_ptr<IO::Logger> _logger)
    {

      logger = std::move(_logger);
    }

    [[nodiscard]] auto& get_logger() const
    {
      return logger;
    }

  private:
    int id{}; ///< The unique identifier to connect with c api.

    /**
     * @brief simulation instance constructor
     *
     * @param argc Number of runtime argument.
     * @param argv Tuntime arguments.
     * @param id A unique identifier for the simulation instance (optional).

     */
    SimulationInstance(int argc, char** argv, std::optional<std::size_t> run_id);

    std::shared_ptr<IO::Logger> logger;

    std::optional<Core::ScalarFactory::ScalarVariant> scalar_initializer_variant = std::nullopt;
    Core::CaseData _data;               ///< Case data for the simulation.
    Core::UserControlParameters params; ///< User-defined control parameters.
    bool loaded = false;                ///< Flag indicating if the instance is loaded.
    bool applied = false;               ///< Flag indicating if the configuration is applied.
    bool registered = false;            ///< Flag indicating if resources are registered.
    std::optional<Simulation::Feed::SimulationFeed> feed =
        std::nullopt; ///< Optional feed configuration.
    std::optional<Simulation::MassTransfer::Type::MtrTypeVariant> mtr_type = std::nullopt;
  };

} // namespace Api
#endif
