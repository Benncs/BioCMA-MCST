#ifndef __SIMULATION_MIXTURE_HPP__
#define __SIMULATION_MIXTURE_HPP__
#include <optional>
#include <string>

namespace Simulation
{

  struct Phase
  {
    std::string name_phase_1;
    double rho;
  };

  class Mixture
  {
    Phase continuous;
    std::optional<Phase> dispersed;

    /** @brief New two phase flow mixture */
    static Mixture new_tpf();
    /** @brief New momophasic mixture */
    static Mixture new_mf();

    double alpha();
    double rho();

  private:
    Mixture() = default;
  };

} // namespace Simulation

#endif //__SIMULATION_MIXTURE_HPP__