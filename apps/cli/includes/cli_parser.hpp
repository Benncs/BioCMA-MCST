#ifndef __CLI_PARSER_HPP__
#define __CLI_PARSER_HPP__

#include <common/common.hpp>
#include <iostream>
#include <optional>

/**
 * @brief Parses command-line arguments to extract simulation parameters.
 *
 * This function processes the command-line arguments provided to the
 * application, extracting and validating the necessary simulation parameters.
 * If the parameters are successfully parsed and valid, they are returned as a
 * `std::optional<SimulationParameters>`. If the parsing fails or the parameters
 * are invalid, an empty `std::optional` is returned.
 *
 * @param argc The number of command-line arguments.
 * @param argv The array of command-line arguments.
 * @return A `std::optional<SimulationParameters>` containing the parsed
 * parameters if successful, or an empty `std::optional` if parsing fails or
 * parameters are invalid.
 * @exception noexcept This function does not throw exceptions.
 */
std::optional<SimulationParameters> parse_cli(int argc, char **argv) noexcept;

/**
 * @brief Print Help message to specified buffer
 */
void showHelp(std::ostream &os) noexcept;

#endif //__CLI_PARSER_HPP__
