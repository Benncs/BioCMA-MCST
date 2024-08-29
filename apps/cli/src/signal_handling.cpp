#include <signal_handling.hpp>

#include <csignal>

SignalHandler *SignalHandler::instance = nullptr;

SignalHandler::SignalHandler() : do_save(false)
{
  instance = this;

  std::signal(SIGUSR1, &SignalHandler::handle_SIGUSR1);

  instance->do_save=false;
}

void SignalHandler::handle_SIGUSR1(int /*unused*/)
{
  if (instance != nullptr)
  {
    instance->do_save = true;
  }
}