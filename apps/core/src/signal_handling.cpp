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
    std::signal(SIGINT, &SignalHandler::handle_SIGINT); // Delete usr2 for not implemented signal

    instance->f_usr1_raised = false;
    instance->f_usr2_raised = false;
    instance->f_sigint_raised = false;
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
    if (instance != nullptr)
    {
      instance->f_usr2_raised = true;
    }
  }

  void SignalHandler::handle_SIGINT(int /*unused*/) noexcept
  {
    if (instance != nullptr)
    {
      instance->f_sigint_raised = true;
    }
  }

} // namespace Core