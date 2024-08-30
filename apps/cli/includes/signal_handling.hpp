#ifndef __SIGNAL_HANDLING_HPP__
#define __SIGNAL_HANDLING_HPP__

#include <stdexcept>

/**
 * @class SignalHandler
 * @brief Manages POSIX SIGNALS.
 */
class SignalHandler
{
public:
  explicit SignalHandler();
  [[nodiscard]] static inline bool is_usr1_raised()
  {
    if (instance == nullptr)
    {
      throw std::runtime_error(
          "Signal not handler not initialized before before used");
    }

    const auto ret = instance->f_usr1_raised;
    instance->f_usr1_raised = false;
    return ret;
  }

private:
  static void handle_SIGUSR1(int) noexcept;
  bool f_usr1_raised; //flag to check if ursr1 signal has been triggered 
  static SignalHandler *instance;
};

#endif //__SIGNAL_HANDLING_HPP__