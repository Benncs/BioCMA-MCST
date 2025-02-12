#include "simulation/mass_transfer.hpp"
#include <cassert>
#include <hydro/impl_mass_transfer.hpp>
#include <scalar_simulation.hpp>
#include <utility>

namespace
{

  constexpr double c_kinematic_viscosity(double temp);

  Eigen::ArrayXXd get_gas_fraction(const CmaRead::ReactorState& state)
  {
    const auto gas_array = Eigen::Map<Eigen::ArrayXd>(const_cast<double*>(state.gasVolume.data()),
                                                      state.gasVolume.size());

    const auto liq_array = Eigen::Map<Eigen::ArrayXd>(
        const_cast<double*>(state.liquidVolume.data()), state.liquidVolume.size());

    return (gas_array / (liq_array + gas_array));
  }

  // Code from
  // https://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html
  constexpr double c_kinematic_viscosity(double temp)
  {

    assert(temp > 0 && temp < 370);

    if (temp < 85)
    {
      return std::round((0.00000000000282244333 * std::pow(temp, 6) -
                         0.00000000126441088087 * std::pow(temp, 5) +
                         0.00000023336659710795 * std::pow(temp, 4) -
                         0.0000234079044336466 * std::pow(temp, 3) +
                         0.00144686943485654 * std::pow(temp, 2) - 0.0607310297913931 * temp +
                         1.79194000343777) *
                        0.000001 * 10000000000) /
             10000000000;
    }

    return std::round((0.00000000000000178038 * std::pow(temp, 6) -
                       0.00000000000277495333 * std::pow(temp, 5) +
                       0.00000000181964246491 * std::pow(temp, 4) -
                       0.00000064995487357883 * std::pow(temp, 3) +
                       0.000136367622445752 * std::pow(temp, 2) - 0.0166081298727911 * temp +
                       1.08486933174497) *
                      0.000001 * 10000000000) /
           10000000000;
  }

} // namespace

namespace Simulation::MassTransfer
{
  template <typename T, typename G>
  Eigen::MatrixXd
  impl_mtr(const T& kla, const Eigen::ArrayXd& Cl, const G& Cs, const Eigen::MatrixXd& Vliq)
  {
    return (kla * (Cs - Cl)).matrix() * Vliq;
  }

  namespace Impl
  {
    constexpr double oxygen_diffusion_constant = 1e-9;
    constexpr double temperature = 20.;

    void flowmap_gas_liquid_mass_transfer(MassTransferProxy& mtr,
                                          const std::shared_ptr<ScalarSimulation>& liquid_scalar,
                                          const std::shared_ptr<ScalarSimulation>& gas_scalar,
                                          const CmaRead::ReactorState& state)
    {
      const double kinematic_viscosity = c_kinematic_viscosity(temperature);

      const double schmidtnumber = kinematic_viscosity / oxygen_diffusion_constant;

      constexpr double db = 5e-3;
      const auto energy_dissipation_array = Eigen::Map<Eigen::ArrayXd>(
          const_cast<double*>(state.energy_dissipation.data()), state.energy_dissipation.size());

#define kl_array                                                                                   \
  (0.3 * (energy_dissipation_array * kinematic_viscosity).pow(0.25) * std::pow(schmidtnumber, -0.5))

#define interfacial_area (6 * get_gas_fraction(state) / db)

#define res_kla_array (kl_array * interfacial_area).transpose()

      mtr.kla.row(1) = res_kla_array.transpose();

      
#define c_star (3.181e-2 * gas_scalar->getConcentrationArray())

      mtr.mtr = Simulation::MassTransfer::impl_mtr(
          mtr.kla, liquid_scalar->getConcentrationArray(), c_star, liquid_scalar->getVolume());
    }

  } // namespace Impl

} // namespace Simulation::MassTransfer
