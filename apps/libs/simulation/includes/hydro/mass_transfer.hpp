#ifndef __SIMULATION_HYDRO_MASS_TRANSFER_HPP__
#define __SIMULATION_HYDRO_MASS_TRANSFER_HPP__

#include <cma_read/reactorstate.hpp>
#include <mc/domain.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Simulation{
  class ScalarSimulation;
};

namespace MassTransfer
{

  struct KlFonctor
  {
    double operator()(double kl);
    Eigen::ArrayXXd operator()(const Eigen::ArrayXd& energy_dissipation_array,
                               const Eigen::ArrayXd& diffusion);
    double operator()(double kla, bool _a);
  };

  struct InterfacialAreaFonctor
  {
    double operator()()
    {
      return 1;
    }
    double operator()(double a);
    Eigen::ArrayXXd operator()(const Eigen::ArrayXd& a);
    Eigen::ArrayXXd operator()(double diameter, const Eigen::ArrayXd& a);
  };

  struct EquilibriumFunctor
  {
    double operator()(double cstar);
    Eigen::ArrayXXd operator()(const Eigen::ArrayXXd& cstar);
    Eigen::ArrayXXd operator()(const Eigen::ArrayXXd& henry, bool _henry);
  };

  /** @brief Model class combining a Coefficient and Equilibrium*/
  Eigen::MatrixXd mass_transfer_rate(Eigen::ArrayXXd& res_kla,
                                     const Eigen::MatrixXd& Vliq,
                                     const Eigen::ArrayXXd& liq_scalar_as_array,
                                     const Eigen::ArrayXXd& gas_scalar_as_array,
                                     const CmaRead::ReactorState& state);

} // namespace MassTransfer



Eigen::MatrixXd gas_liquid_mass_transfer(Eigen::ArrayXXd &res_kla,
                                         const Eigen::MatrixXd &Vliq,
                                         const Eigen::ArrayXXd &liq_scalar_as_array,
                                         const Eigen::ArrayXXd &gas_scalar_as_array,
                                         const CmaRead::ReactorState &state);

void gas_liquid_mass_transfer(Simulation::ScalarSimulation* liquid_scalar,
                                             Simulation::ScalarSimulation* gas_scalar,
                                             const CmaRead::ReactorState& state);

#endif //__SIMULATION_HYDRO_MASS_TRANSFER_HPP__ __SIM_HYDRO_MASS_TRANSFER__
