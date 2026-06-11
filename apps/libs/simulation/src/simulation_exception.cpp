#include <simulation/simulation_exception.hpp>
#if defined(_MSC_VER) && !defined(__clang__)
#  define _unreachable __assume(false);
#else
#  define _unreachable __builtin_unreachable();
#endif
namespace Simulation
{
  BaseSimulationException::BaseSimulationException(ErrorCodes code,
                                                   std::source_location loc)
      : code_(code), loc_(loc)
  {
  }

  [[nodiscard]] const char*
  BaseSimulationException::what() const noexcept
  {
    switch (code_)
    {
    case ErrorCodes::BadConcentrationInitLiq:
      return "Simulation: UDF liquid concentration returned negative values";
    case ErrorCodes::BadConcentrationInitGas:
      return "Simulation: UDF gas concentration returned negative values";
    case ErrorCodes::MismatchSizeVolume:
      return "Simulation: volume size does not match compartment count";
    case ErrorCodes::BadInitialiser:
      return "Simulation: bad initialiser";
    case ErrorCodes::MismatchSize:
      return "Simulation: concentration array size mismatch";
    case ErrorCodes::Unknown:
      return "Simulation: unknown error";
    }
    _unreachable
  }

  [[nodiscard]] const std::source_location&
  BaseSimulationException::where() const noexcept
  {
    return loc_;
  }

  FeedException::FeedException(FeedExceptionError code,
                               std::source_location loc)
      : code_(code), loc_(loc)
  {
  }

  [[nodiscard]] const char*
  FeedException::what() const noexcept
  {
    switch (code_)
    {
    case FeedExceptionError::NegativeFlow:
      return "Feed: negative flow value";
    case FeedExceptionError::NegativeConcentration:
      return "Feed: negative concentration value";
    case FeedExceptionError::OutofRange:
      return "Feed: position out of range";
    }
    _unreachable
  }

  [[nodiscard]] const std::source_location&
  FeedException::where() const noexcept
  {
    return loc_;
  }

} // namespace Simulation
