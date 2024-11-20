#ifndef __SIGNAL_HANDLING_HPP__
#define __SIGNAL_HANDLING_HPP__

#include <stdexcept>

namespace Core
{

  /**
   * @class SignalHandler
   * @brief Manages POSIX signals for the application.
   *
   * Singleton that captures specific POSIX signals
   * (e.g., SIGUSR1 and SIGUSR2) during the runtime of the application. It provides mechanisms
   * to check if signals have been raised and ensures proper cleanup and handling.
   */
  class SignalHandler
  {
  public:
    /**
     * @brief Constructor for SignalHandler.
     *
     * Sets up the signal handlers for SIGUSR1 and SIGUSR2. Throws a runtime error
     * if more than one instance of SignalHandler is created.
     */
    explicit SignalHandler();

    /**
     * @brief Checks if the SIGUSR1 signal has been raised.
     *
     * This function is used to query the state of the SIGUSR1 signal. If the signal has
     * been raised, the function returns `true` and resets the internal flag to `false`.
     *
     * @throws std::runtime_error If the SignalHandler instance has not been initialized.
     * @return `true` if SIGUSR1 was raised since the last query; `false` otherwise.
     */
    [[nodiscard]] static inline bool is_usr1_raised() noexcept(false)
    {
      if (instance == nullptr)
      {
        throw std::runtime_error("SignalHandler not initialized before use");
      }

      const auto ret = instance->f_usr1_raised;
      instance->f_usr1_raised = false;
      return ret;
    }

  private:
    /**
     * @brief Handles the SIGUSR1 signal.
     *
     * This static function is called when the SIGUSR1 signal is raised. It sets the
     * internal flag to indicate the signal was received.
     *
     * @param signal The signal number (ignored in this implementation).
     */
    static void handle_SIGUSR1(int signal) noexcept;

    /**
     * @brief Handles the SIGUSR2 signal.
     *
     * This static function is intended to handle the SIGUSR2 signal. (Implementation pending.)
     *
     * @param signal The signal number (ignored in this implementation).
     */
    static void handle_SIGUSR2(int signal) noexcept;

    bool f_usr1_raised;             ///< Flag to indicate if the SIGUSR1 signal was triggered.
    static SignalHandler* instance; ///< Singleton instance of the SignalHandler.
  };

} // namespace Core

#endif //__SIGNAL_HANDLING_HPP__