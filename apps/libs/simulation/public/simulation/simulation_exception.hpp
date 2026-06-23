#ifndef __SIMULATION_EXCEPTION_HPP__
#define __SIMULATION_EXCEPTION_HPP__

#include <exception>
#include <source_location>

namespace Simulation
{

  enum class ErrorCodes
  {
    Unknown,
    BadConcentrationInitLiq,
    BadConcentrationInitGas,
    BadInitialiser,
    MismatchSize,
    MismatchSizeVolume,
  };

  enum class FeedExceptionError
  {
    NegativeFlow,
    NegativeConcentration,
    OutofRange
  };

  class SimulationException : public std::exception
  {
  public:
    [[nodiscard]] virtual const std::source_location& where() const noexcept
        = 0;
  };

  class BaseSimulationException : public SimulationException
  {
  public:
    explicit BaseSimulationException(ErrorCodes code,
                                     std::source_location loc
                                     = std::source_location::current());
    [[nodiscard]] const char* what() const noexcept override;
    [[nodiscard]] const std::source_location& where() const noexcept override;

  private:
    ErrorCodes code_;
    std::source_location loc_;
  };

  class FeedException : public SimulationException
  {
  public:
    explicit FeedException(FeedExceptionError code,
                           std::source_location loc
                           = std::source_location::current());
    [[nodiscard]] const char* what() const noexcept override;
    [[nodiscard]] const std::source_location& where() const noexcept override;

  private:
    FeedExceptionError code_;
    std::source_location loc_;
  };
} // namespace Simulation

#endif
