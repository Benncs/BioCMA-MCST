#ifndef __CORE_SCALAR_FACTORY_HPP__
#define __CORE_SCALAR_FACTORY_HPP__

#include <cstddef>
#include <optional>
#include <simulation/scalar_initializer.hpp>
#include <span>
#include <string_view>
#include <variant>
#include <vector>

/**
 * @namespace ScalarFactory
 * @brief Provides structures and functions to initialize scalar data for simulations, supporting multiple data sources.
 *
 * The ScalarFactory namespace contains several structures representing different methods of initializing scalar values
 * (liquid and gas concentrations) for simulations.
 */
namespace Core::ScalarFactory
{
  /**
   * @struct Uniform
   * @brief Represents uniform concentration data for all compartments.
   */
  struct Uniform
  {
    /**
     * @brief Vector of liquid phase concentrations for all compartments.
     */
    std::vector<double> liquid_concentrations;

    /**
     * @brief Optional vector of gas phase concentrations, applied uniformly if provided.
     */
    std::optional<std::vector<double>> gas_concentration = std::nullopt;
  };

  /**
   * @struct Local
   * @brief Represents localized concentration data for specific compartments.
   */
  struct Local
  {
    /**
     * @brief Vector of liquid phase concentrations for specified compartments.
     */
    std::vector<double> liquid_concentrations;

    /**
     * @brief Vector of indices corresponding to compartments with specified liquid concentrations.
     */
    std::vector<size_t> liquid_indices;

    /**
     * @brief Optional vector of gas phase concentrations, applied to specific compartments if provided.
     */
    std::optional<std::vector<double>> gas_concentration = std::nullopt;

    /**
     * @brief Optional vector of indices corresponding to compartments with specified gas concentrations.
     */
    std::optional<std::vector<size_t>> gas_indices = std::nullopt;
  };

  /**
   * @struct File
   * @brief Represents scalar data sourced from an external file.
   */
  struct File
  {
    /**
     * @brief Number of compartments represented in the file.
     */
    size_t n_compartment;

    /**
     * @brief Path to the file containing scalar data for compartments.
     */
    std::string_view path;
  };

  /**
   * @struct CustomScript
   * @brief Represents scalar data initialized through a custom script.
   */
  struct CustomScript
  {
    /**
     * @brief Path to the custom script used for initializing scalar data.
     */
    std::string_view path;
  };

  /**
   * @typedef ScalarVariant
   * @brief A variant type representing different scalar data sources.
   */
  using ScalarVariant = std::variant<Uniform, Local, File, CustomScript>;

  /**
   * @struct Visitor
   * @brief Functor struct to visit each type in ScalarVariant and initialize scalar data.
   */
  struct Visitor
  {

    Simulation::ScalarInitializer operator()(Uniform args);

    Simulation::ScalarInitializer operator()(Local args);

    Simulation::ScalarInitializer operator()(File filepath);

    Simulation::ScalarInitializer operator()(CustomScript path);
  };

  /**
   * @brief Factory function to initialize scalar data based on the specified input source.
   *
   * The `scalar_factory` function uses a provided variant type to initialize scalar data for both liquid
   * and optionally gas phases. This function accommodates a variety of data sources for flexibility in setup.
   *
   * @param f_init_gas_flow Boolean flag indicating whether gas flow initialization is required.
   * @param gas_volume Span of gas volume data for each compartment.
   * @param liquid_volume Span of liquid volume data for each compartment.
   * @param arg_liq Scalar variant specifying the source of scalar data for liquid concentration.
   * @return A `ScalarInitializer` configured with the provided data.
   */
  Simulation::ScalarInitializer scalar_factory(bool f_init_gas_flow,
                                               std::span<double> gas_volume,
                                               std::span<double> liquid_volume,
                                               ScalarVariant arg_liq);

  bool sanitize(const Simulation::ScalarInitializer &res);
} // namespace Core::ScalarFactory

#endif