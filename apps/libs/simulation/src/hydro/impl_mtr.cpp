#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>
#include <hydro/impl_mass_transfer.hpp>
#include <scalar_simulation.hpp>
#include <simulation/mass_transfer.hpp>

namespace
{

  constexpr double c_kinematic_viscosity(double temp);

  inline auto get_gas_fraction(const CmaUtils::IterationState& state)
  {
    const auto gas_array = Eigen::Map<Eigen::ArrayXd>(const_cast<double*>(state.gas->volume.data()),
                                                      state.gas->volume.size());

    const auto liq_array = Eigen::Map<Eigen::ArrayXd>(const_cast<double*>(state.liq->volume.data()),
                                                      state.liq->volume.size());

    return (gas_array / (liq_array + gas_array));
  }

  inline auto get_interfacial_area(const double db, const CmaUtils::IterationState& state)
  {
    return 6. * get_gas_fraction(state) / db;
  }

  template <typename G, typename K>
  inline auto kl_correlation(const double schmidtnumber,
                             const G& kinematic_viscosity,
                             const K& energy_dissipation_array)
  {
    return 0.3 * (energy_dissipation_array * kinematic_viscosity).pow(0.25) *
           std::pow(schmidtnumber, -0.5);
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
                                          const Eigen::ArrayXXd& liquid_concentration,
                                          const Eigen::ArrayXXd& gas_concentration,
                                          const Eigen::MatrixXd& liquid_volume,
                                          const CmaUtils::IterationState& state)
    {

      const double kinematic_viscosity = c_kinematic_viscosity(temperature);

      const double schmidtnumber = kinematic_viscosity / oxygen_diffusion_constant;
      const auto eps = state.infos.at("energy_dissipation");
      const Eigen::Map<Eigen::ArrayXd> energy_dissipation_array =
          Eigen::Map<Eigen::ArrayXd>(const_cast<double*>(eps.data()), EIGEN_INDEX(eps.size()));


      assert(gas_concentration.stride() == liquid_concentration.stride());
      assert(energy_dissipation_array.rows() == liquid_concentration.cols());

      mtr.kla.row(1) =
          (kl_correlation(schmidtnumber, kinematic_viscosity, energy_dissipation_array) *
           get_interfacial_area(mtr.db, state))
              .transpose();

      mtr.mtr = impl_mtr(
          mtr.kla, liquid_concentration, gas_concentration.colwise() * mtr.Henry, liquid_volume);
    }

    void fixed_kla_gas_liquid_mass_transfer(MassTransferProxy& mtr,
                                            const Eigen::ArrayXXd& liquid_concentration,
                                            const Eigen::ArrayXXd& gas_concentration,
                                            const Eigen::MatrixXd& liquid_volume,
                                            const CmaUtils::IterationState& state)
    {
      (void)state;
      mtr.mtr = impl_mtr(
          mtr.kla, liquid_concentration, gas_concentration.colwise() * mtr.Henry, liquid_volume);
    }

  } // namespace Impl

} // namespace Simulation::MassTransfer
