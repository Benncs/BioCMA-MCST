#ifndef __RUNTIME_INIT_HPP__
#define __RUNTIME_INIT_HPP__

#include <chrono>
#include <common/execinfo.hpp>
#include <core/simulation_parameters.hpp>
#include <ctime>
#include <iomanip>
#include <string>
#include <string_view>

// [[deprecated]]void set_n_thread_current_rank(int rank,
//                                int size,
//                                ExecInfo &info,
//                                const Core::UserControlParameters &params) noexcept;

/**
 * @brief Appends the current date and time to the given stream.
 *

 * @tparam Stream The type of the stream (e.g., `std::ostream`,
 `std::stringstream`).
 * @param stream The output stream to which the date and time will be appended.
 */
template <typename Stream> void append_date_time(Stream& stream) noexcept
{
  auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  stream << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S");
}

// /**
//  * @brief Appends the current date and time to a given string.
//  *
//  * @param string The string view to which the date and time will be appended.
//  * @return A `std::string` containing the original string with the appended date
//  * and time.
//  */
// std::string sappend_date_time(std::string_view string) noexcept;

// /**
//  * @brief Initializes the runtime environment based on command-line arguments
//  * and simulation parameters.
//  *
//  * This function sets up the necessary runtime environment for the simulation
//  * by:
//  * - Initializing MPI (Message Passing Interface) if applicable.
//  * - Setting up Kokkos for parallel programming.
//  * - Configuring functions to be executed upon program exit.
//  * - Handling signals to ensure proper shutdown and resource cleanup.
//  *
//  * The function uses the provided command-line arguments and simulation
//  * parameters to configure the runtime environment accordingly.
//  *
//  * @param argc The number of command-line arguments.
//  * @param argv The array of command-line arguments.
//  * @param params The `UserControlParameters` object containing configuration
//  * settings for the simulation.
//  * @return An `ExecInfo` object containing details about the initialized runtime
//  * environment, including execution context and other relevant metadata.
//  */
// ExecInfo
// runtime_init(int argc, char **argv, Core::UserControlParameters &params) noexcept;

/**
 * @brief Print run metadata to log file before running */
void register_run(const ExecInfo& exec, const Core::UserControlParameters& params) noexcept;

#endif //__RUNTIME_INIT_HPP__
