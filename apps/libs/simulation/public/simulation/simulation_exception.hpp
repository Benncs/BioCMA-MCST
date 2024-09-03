#ifndef __SIMULATION_EXCEPTION_HPP__
#define __SIMULATION_EXCEPTION_HPP__

#include <exception>
#include <stdexcept>

namespace Simulation
{

  enum class ErrorCodes
  {
    BadConcentrationInitLiq,
    BadConcentrationInitGas
  };

  class SimulationException : public std::exception
  {
  public:
    explicit SimulationException(ErrorCodes code);

    [[nodiscard]] const char *what() const noexcept override
    {

      return msg.c_str();
    }

  private:
    ErrorCodes code_;
    std::string msg;
    [[nodiscard]] std::string getMessage() const;
  };
} // namespace Simulation

#endif