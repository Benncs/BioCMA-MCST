#ifndef __IMPL_HYDRO_MASS_TRANSFER_HPP__
#define __IMPL_HYDRO_MASS_TRANSFER_HPP__

#include "cma_utils/iteration_state.hpp"
#include "eigen_kokkos.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <mc/domain.hpp>
#include <simulation/mass_transfer.hpp>

namespace Simulation
{
  class ScalarSimulation;
};

namespace Simulation::MassTransfer
{
  struct MassTransferProxy
  {
    ColMajorMatrixtype mtr;
    Eigen::ArrayXXd kla;
    Eigen::ArrayXd Henry;
    double db;
  };

  namespace Impl
  {

    void flowmap_gas_liquid_mass_transfer(MassTransferProxy& mtr,
                                          const Eigen::ArrayXXd& liquid_concentration,
                                          const Eigen::ArrayXXd& gas_concentration,
                                          const Eigen::MatrixXd& liquid_volume,
                                          const CmaUtils::IterationState& state);

    void fixed_kla_gas_liquid_mass_transfer(MassTransferProxy& mtr,
                                            const Eigen::ArrayXXd& liquid_concentration,
                                            const Eigen::ArrayXXd& gas_concentration,
                                            const Eigen::MatrixXd& liquid_volume,
                                            const CmaUtils::IterationState& state);

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
