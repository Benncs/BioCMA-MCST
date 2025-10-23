#ifndef __CLI_PARSER_HPP__
#define __CLI_PARSER_HPP__

#include <common/logger.hpp>
#include <common/results.hpp>
#include <core/simulation_parameters.hpp>
#include <memory>
#include <utility>

// TODO WIP
template <typename S> struct CliResults : Result<S, std::string>
{
  explicit CliResults(std::string_view t) noexcept
      : Result<S, std::string>(std::string(t))
  {
  }

  explicit CliResults(S&& value) noexcept
      : Result<S, std::string>(std::move(value))
  {
  }

  explicit constexpr CliResults() noexcept = default;

  explicit operator CliResults<S>() &&
  {
    return CliResults<S>(std::move(*this));
  }

  explicit operator CliResults<S>() const&
  {
    return CliResults<S>(this->get());
  }
};

/**
 * @brief Parses command-line arguments to extract simulation parameters.
 *
 * This function processes the command-line arguments provided to the
 * application, extracting and validating the necessary simulation parameters.
 * If the parameters are successfully parsed and valid, they are returned as a
 * `std::optional<UserControlParameters>`. If the parsing fails or the
 * parameters are invalid, an empty `std::optional` is returned.
 *
 * @param argc The number of command-line arguments.
 * @param argv The array of command-line arguments.
 * @return A `std::optional<UserControlParameters>` containing the parsed
 * parameters if successful, or an empty `std::optional` if parsing fails or
 * parameters are invalid.
 * @exception noexcept This function does not throw exceptions.
 */
CliResults<Core::UserControlParameters> parse_cli(
    const std::shared_ptr<IO::Logger>& logger, int argc, char** argv) noexcept;

/**
 * @brief Print Help message to specified buffer
 */
std::string get_help_message() noexcept;

#endif //__CLI_PARSER_HPP__
