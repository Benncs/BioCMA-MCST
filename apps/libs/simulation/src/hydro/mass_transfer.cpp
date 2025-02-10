#include "simulation/mass_transfer.hpp"
#include <cassert>
#include <hydro/impl_mass_transfer.hpp>
#include <scalar_simulation.hpp>

static constexpr double c_kinematic_viscosity(double temp);
static constexpr double oxygen_diffusion_constant = 1e-9;
static constexpr double temperature = 20.;

namespace MassTransfer
{
  template <typename T, typename G>
  Eigen::MatrixXd
  impl_mtr(const T& kla, const Eigen::ArrayXd& Cl, const G& Cs, const Eigen::MatrixXd& Vliq)
  {
    return (kla * (Cs - Cl)).matrix() * Vliq;
  }

  Eigen::MatrixXd mass_transfer_rate(Eigen::ArrayXXd& res_kla,
                                     const Eigen::MatrixXd& Vliq,
                                     const Eigen::ArrayXXd& liq_scalar_as_array,
                                     const Eigen::ArrayXXd& gas_scalar_as_array,
                                     const CmaRead::ReactorState& state)
  {
    return impl_mtr(res_kla, liq_scalar_as_array, gas_scalar_as_array, Vliq);
  }

}; // namespace MassTransfer

static Eigen::ArrayXXd get_gas_fraction(const CmaRead::ReactorState& state)
{
  const auto gas_array = Eigen::Map<Eigen::ArrayXd>(const_cast<double*>(state.gasVolume.data()),
                                                    state.gasVolume.size());

  const auto liq_array = Eigen::Map<Eigen::ArrayXd>(const_cast<double*>(state.liquidVolume.data()),
                                                    state.liquidVolume.size());

                                                    
  return (gas_array / (liq_array + gas_array));
}

namespace
{
  void flowmap_gas_liquid_mass_transfer(Simulation::ScalarSimulation* liquid_scalar,
                                        Simulation::ScalarSimulation* gas_scalar,
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

    liquid_scalar->vec_kla.row(1) = res_kla_array.transpose();
#define c_star (3.181e-2 * gas_scalar->getConcentrationArray())

    liquid_scalar->set_mass_transfer(MassTransfer::impl_mtr(liquid_scalar->vec_kla,
                                                            liquid_scalar->getConcentrationArray(),
                                                            c_star,
                                                            liquid_scalar->getVolume()));
  }



} // namespace

void gas_liquid_mass_transfer(MassTransfer::MTRType type,Simulation::ScalarSimulation* liquid_scalar,
                              Simulation::ScalarSimulation* gas_scalar,
                              const CmaRead::ReactorState& state)
{
  switch (type) {
    case MassTransfer::MTRType::FixedKla:
    {
      
    };
    case MassTransfer::MTRType::Flowmap:
    {
        flowmap_gas_liquid_mass_transfer(liquid_scalar,gas_scalar,state);
    }
    default:
    {

    };
  }
  
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
