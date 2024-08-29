#ifndef __SIGNAL_HANDLING_HPP__
#define __SIGNAL_HANDLING_HPP__

class SignalHandler
{
public:
  explicit SignalHandler();
  [[nodiscard]] static inline bool is_usr1_raised()
  {
    const auto ret= instance->do_save;
    instance->do_save = false;
    return ret;
  }

private:
  static void handle_SIGUSR1(int);
  bool do_save;
  static SignalHandler *instance;
};

#endif //__SIGNAL_HANDLING_HPP__