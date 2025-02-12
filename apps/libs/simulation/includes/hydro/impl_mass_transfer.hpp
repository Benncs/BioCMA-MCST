#ifndef __IMPL_HYDRO_MASS_TRANSFER_HPP__
#define __IMPL_HYDRO_MASS_TRANSFER_HPP__

#include "scalar_simulation.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cma_read/reactorstate.hpp>
#include <mc/domain.hpp>
#include <memory>
#include <simulation/mass_transfer.hpp>

namespace Simulation
{
  class ScalarSimulation;
};

namespace Simulation::MassTransfer
{
  struct MassTransferProxy
  {
    MatrixType mtr;
    Eigen::ArrayXXd kla;
  };
  
  namespace Impl
  {

    void flowmap_gas_liquid_mass_transfer(MassTransferProxy& mtr,
                                          const std::shared_ptr<ScalarSimulation>& liquid_scalar,
                                          const std::shared_ptr<ScalarSimulation>& gas_scalar,
                                          const CmaRead::ReactorState& state);
  }; // namespace Impl

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

} // namespace Simulation::MassTransfer

#endif //__SIMULATION_HYDRO_MASS_TRANSFER_HPP__ __SIM_HYDRO_MASS_TRANSFER__
