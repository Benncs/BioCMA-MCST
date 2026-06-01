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
//                                const Core::UserControlParameters &params)
//                                noexcept;

/**
 * @brief Appends the current date and time to the given stream.
 *

 * @tparam Stream The type of the stream (e.g., `std::ostream`,
 `std::stringstream`).
 * @param stream The output stream to which the date and time will be appended.
 */
template <typename Stream>
void
append_date_time(Stream& stream) noexcept
{
  auto now
      = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  stream << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S");
}

/**
 * @brief Print run metadata to log file before running */
void register_run(const ExecInfo& exec,
                  const Core::UserControlParameters& params) noexcept;

#endif //__RUNTIME_INIT_HPP__
