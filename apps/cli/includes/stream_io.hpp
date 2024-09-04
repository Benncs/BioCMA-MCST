#ifndef __STREAM_IO_HPP__
#define __STREAM_IO_HPP__

#include <streambuf>
/**
 * @namespace FlagCompileTIme
 * @brief Namespace containing compile-time flags for application configuration.
 *
 * This namespace defines compile-time constants that control specific behavior
 * of the application, such as verbosity and output redirection. These flags are
 * set based on predefined macros and allow for conditional compilation.
 */
namespace FlagCompileTIme
{
  /**
   * @var verbose
   * @brief Controls whether verbose logging is enabled.
   *
   * The `verbose` flag is set to `true` if `__FLAG_APP_VERBOSE__` is defined at
   * compile-time; otherwise, it is set to `false`.
   *
   * @var __f__redirect
   * @brief Controls whether standard output should be redirected.
   *
   * The `__f__redirect` flag is set to `true` if `__FLAG_APP_REDIRECT_STDOUT__`
   * is defined at compile-time; otherwise, it is set to `false`.
   */



#ifndef __FLAG_APP_REDIRECT_STDOUT__
  constexpr bool __f__redirect = false;
#else
  constexpr bool __f__redirect = true;
#endif
} // namespace FlagCompileTIme

/**
 * @brief Redirects the standard output to a custom stream.
 *
 * This function redirects the standard output (stdout) to a provided string
 * stream. It stores the original buffer of stdout so that it can be restored
 * later.
 *
 * @param original_buffer A reference to a pointer that will store the original
 *                        stdout buffer. This allows the buffer to be restored
 * later.
 * @param variable_stream The string stream to which stdout will be redirected.
 * @return The file descriptor of the original stdout, which can be used to
 * restore it later.
 */
int redirect_stdout(std::streambuf *&original_buffer,
                    std::stringstream &variable_stream);

/**
 * @brief Restores the standard output to its original state.
 *
 * This function restores the standard output (stdout) to its original buffer,
 * undoing the redirection performed by `redirect_stdout`.
 *
 * @param original_stdout The file descriptor of the original stdout, obtained
 * from `redirect_stdout`.
 * @param original_buffer A reference to the pointer to the original stdout
 * buffer. This will be set back to stdout.
 */
void restore_stdout(int original_stdout, std::streambuf *&original_buffer);

// C Style Macro to hide implementation detail in Main
// If redirect, capture the stdout during block execution in a buffer and print
// it at the end of the block
#define REDIRECT_BLOCK(__block__)                                              \
  if constexpr ((FlagCompileTIme::__f__redirect))                              \
  {                                                                            \
    std::stringstream output_variable;                                         \
    output_variable.str("");                                                   \
    std::streambuf *original_buffer = nullptr;                                 \
    auto original_stdout_fd =                                                  \
        redirect_stdout(original_buffer, output_variable);                     \
    {__block__};                                                               \
    restore_stdout(original_stdout_fd, original_buffer);                       \
    if constexpr ((FlagCompileTIme::verbose))                                  \
    {                                                                          \
      std::cout << output_variable.str() << std::endl;                         \
    }                                                                          \
  }                                                                            \
  else                                                                         \
  {                                                                            \
    __block__                                                                  \
  }
#endif