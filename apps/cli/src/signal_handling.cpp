#include <signal_handling.hpp>

#include <csignal>

SignalHandler *SignalHandler::instance = nullptr;

SignalHandler::SignalHandler() : f_usr1_raised(false)
{
  instance = this;

  std::signal(SIGUSR1, &SignalHandler::handle_SIGUSR1);

  instance->f_usr1_raised = false;
}

void SignalHandler::handle_SIGUSR1(int /*unused*/) noexcept
{
  if (instance != nullptr)
  {
    instance->f_usr1_raised = true;
  }
}