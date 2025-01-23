#include <simulation/simulation_exception.hpp>

namespace Simulation
{
  SimulationException::SimulationException(ErrorCodes code) : code_(code)
  {
    msg = this->getMessage();
  }

  [[nodiscard]] std::string SimulationException::getMessage() const
  {
    switch (code_)
    {
    case ErrorCodes::BadConcentrationInitLiq:
      return "Simulation::post_init_concentration User Defined Function for "
             "liquid concentration leads to negative values\r\nPlease check "
             "UDF "
             "and restart";

    case ErrorCodes::BadConcentrationInitGas:
      return "Simulation::post_init_concentration User Defined Function for "
             "gas concentration leads to negative values\r\nPlease check "
             "UDF "
             "and restart";

    // default:
    //   return "Error: Unknown error.";
    case ErrorCodes::BadInitialiser:
      return "Bad initialiser";
    case ErrorCodes::MismatchSize:
      return "post_init_concentration: Mismatch Size";
    }
#if defined(_MSC_VER) && !defined(__clang__) // MSVC
    __assume(false);
#else // GCC, Clang
    __builtin_unreachable();
#endif
  }
} // namespace Simulation