#include <cassert>
#include <csignal>
#include <iostream>
#include <signal_handling.hpp>
#include <unistd.h>

void check_sigus1()
{
  assert(!Core::SignalHandler::is_usr1_raised() &&
         "No signal should be raised initially");
  std::raise(SIGUSR1);
  assert(Core::SignalHandler::is_usr1_raised() &&
         "Signal should have been raised");
  assert(!Core::SignalHandler::is_usr1_raised() &&
         "Signal flag should reset after querying");
}

void check_sigus2()
{
  assert(!Core::SignalHandler::is_usr1_raised() &&
         "No signal should be raised initially");
  std::raise(SIGUSR2);
  assert(!Core::SignalHandler::is_usr1_raised() &&
         "Signal should not have been raised");
}

int main()
{
  try
  {
    auto _ = Core::SignalHandler::is_usr1_raised();
    assert(false);
  }
  catch (std::runtime_error& e)
  {
    // pass
  }
  Core::SignalHandler handler;
  Core::SignalHandler _handler;

  check_sigus1();
  check_sigus2();
  std::cout << "All tests passed!" << std::endl;

  return 0;
}