#include <signal_handling.hpp>

#include <csignal>
#include <iostream>
namespace Core
{

  SignalHandler* SignalHandler::instance = nullptr;

  SignalHandler::SignalHandler() : f_usr1_raised(false)
  {
    instance = this;

    std::signal(SIGUSR1, &SignalHandler::handle_SIGUSR1);
    std::signal(SIGUSR2, &SignalHandler::handle_SIGUSR2); // Delete usr2 for not implemented signal

    instance->f_usr1_raised = false;
  }

  void SignalHandler::handle_SIGUSR1(int /*unused*/) noexcept
  {
    if (instance != nullptr)
    {
      instance->f_usr1_raised = true;
    }
  }

  void SignalHandler::handle_SIGUSR2(int /*unused*/) noexcept
  {
    std::cerr << "SIGUSR2 not implemented" << std::endl;
  }

} // namespace Core