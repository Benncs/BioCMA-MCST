#include <cassert>
#include <hydro/mass_transfer.hpp>
#include <stdexcept>

static constexpr double c_kinematic_viscosity(double temp);
static constexpr double oxygen_diffusion_constant = 1e-9;
static constexpr double temperature = 20.;

Eigen::MatrixXd
gas_liquid_mass_transfer(Eigen::ArrayXXd &res_kla,
                         const Eigen::MatrixXd &Vliq,
                         const Eigen::ArrayXXd &liq_scalar_as_array,
                         const Eigen::ArrayXXd &gas_scalar_as_array,
                         const CmaRead::ReactorState &state)
{

  const double kinematic_viscosity = c_kinematic_viscosity(temperature);

  const double schmidtnumber = kinematic_viscosity / oxygen_diffusion_constant;

  constexpr double db = 1e-3;
  // const auto energy_dissipation_array =
  //     Eigen::Map<Eigen::ArrayXd>(const_cast<double*>(state.energy_dissipation.data()),
  //                                 state.energy_dissipation.size());

  // const auto gas_array = Eigen::Map<Eigen::ArrayXd>(const_cast<double*>(state.gasVolume.data()),
  //                                 state.gasVolume.size());

  // const auto liq_array = Eigen::Map<Eigen::ArrayXd>(const_cast<double*>(state.liquidVolume.data()),
  //                                 state.liquidVolume.size());

  // //LAZY 
  // #define kl_array  (0.3 * (energy_dissipation_array * kinematic_viscosity).pow(0.25) * \
  //           std::pow(schmidtnumber, -0.5))

  // #define gas_fraction_array  ((gas_array)/(liq_array+gas_array))

  // #define interfacial_area  (6*gas_fraction_array/db)

  // #define res_kla kl_array*interfacial_area

  for (int i_c = 0; i_c < static_cast<int>(state.n_compartments); ++i_c)
  {

    double eps_turb{};
    // TODO FIXME: just to use flowmap without turbulence
    try
    {
      eps_turb = state.energy_dissipation.at(i_c);
    }
    catch (const std::out_of_range &e)
    {
      eps_turb = 0.5;
    }
    const double kl = 0.3 * std::pow(eps_turb * kinematic_viscosity, 0.25) *
                      std::pow(schmidtnumber, -0.5);

    const double gas_fraction =
        state.gasVolume[i_c] / (state.liquidVolume[i_c] + state.gasVolume[i_c]);

    const double a = 6 * gas_fraction / db;

    const double kla = kl * a;

    res_kla.coeffRef(1, i_c) = kla;
  }

// LAZY EVALUATION
#define c_star (1.3e-5 * gas_scalar_as_array / 32e-3 * 8.314 * (273.15 + 30))
#define transfer_g_liq (res_kla * (c_star - liq_scalar_as_array))

  return (transfer_g_liq.matrix() * Vliq).eval();
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
                       0.00144686943485654 * std::pow(temp, 2) -
                       0.0607310297913931 * temp + 1.79194000343777) *
                      0.000001 * 10000000000) /
           10000000000;
  }

  return std::round((0.00000000000000178038 * std::pow(temp, 6) -
                     0.00000000000277495333 * std::pow(temp, 5) +
                     0.00000000181964246491 * std::pow(temp, 4) -
                     0.00000064995487357883 * std::pow(temp, 3) +
                     0.000136367622445752 * std::pow(temp, 2) -
                     0.0166081298727911 * temp + 1.08486933174497) *
                    0.000001 * 10000000000) /
         10000000000;
}
