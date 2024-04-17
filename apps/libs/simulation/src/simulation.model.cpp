#include <cma_read/reactorstate.hpp>
#include <mc/domain.hpp>
#include <simulation/simulation.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <scalar_simulation.hpp>

static Eigen::MatrixXd mock_transfer(Eigen::ArrayXXd &res_kla,
                                     const Eigen::MatrixXd &Vliq,
                                     const Eigen::ArrayXXd &liq_scalar_as_array,
                                     const Eigen::ArrayXXd &gas_scalar_as_array,
                                     ReactorState *const state)
{

  double kinematic_viscosity = 1.0023e-06;
  double oxygen_diffusion_constant = 1e-9;

  double schmidtnumber = kinematic_viscosity / oxygen_diffusion_constant;

  double db = 1e-3;

  for (int i_c = 0; i_c < static_cast<int>(state->n_compartments); ++i_c)
  {

    const double eps_turb = state->energy_dissipation[i_c];
    const double kl = 0.3 * std::pow(eps_turb * kinematic_viscosity, 0.25) *
                      std::pow(schmidtnumber, -0.5);

    const double gas_fraction =
        state->gasVolume[i_c] /
        (state->liquidVolume[i_c] + state->gasVolume[i_c]);
    const double a = 6 * gas_fraction / db;

    const double kla = kl * a;

    res_kla.coeffRef(1, i_c) = kla;
  }

//LAZY EVALUATION
#define c_star (1.3e-5 * gas_scalar_as_array / 32e-3 * 8.314 * (273.15 + 30))
#define transfer_g_liq (res_kla * (c_star - liq_scalar_as_array))

  return (transfer_g_liq.matrix() * Vliq).eval();
}
// double kla = 1e-5; //MOCK USING WRONG FLOWMAP

namespace Simulation
{

  void SimulationUnit::post_process_reducing()
  {
    for (size_t i_thread = 0; i_thread < n_thread; ++i_thread)
    {
      this->liquid_scalar->merge(i_thread);
      this->container->merge(i_thread);
    }
  }

  std::span<double> SimulationUnit::get_contributionData()
  {
    return std::span(this->liquid_scalar->biomass_contribution.data(),
                     this->liquid_scalar->biomass_contribution.size());
  }

  std::span<double> SimulationUnit::getCliqData()
  {
    auto *data = this->liquid_scalar->C.data();

    return std::span<double>(data, this->liquid_scalar->C.size());
  }

  void SimulationUnit::reduceContribs(std::span<double> data, size_t n_rank)
  {

    size_t nr = this->liquid_scalar->C.rows();
    size_t nc = this->gas_scalar->C.cols();

    this->liquid_scalar->biomass_contribution.setZero();

    for (int i = 0; i < static_cast<int>(n_rank); ++i)
    {
      this->liquid_scalar->biomass_contribution += Eigen::Map<Eigen::MatrixXd>(
          &data[i * nr * nc], static_cast<int>(nr), static_cast<int>(nc));
    }
  }

  void SimulationUnit::clearContribution()
  {
    for (size_t i = 0; i < n_thread; ++i)
    {
      this->liquid_scalar->contribs[i].setZero();
      container->extras.clear();
      this->liquid_scalar->vec_kla.setZero();
    }

    this->liquid_scalar->biomass_contribution.setZero();
  }

  void SimulationUnit::step(double d_t,ReactorState& state)
  {

    // if (state == nullptr)
    // {
    //   throw std::runtime_error("Error no given reactor state");
    // }

    auto mat_transfer_g_liq = mock_transfer(this->liquid_scalar->vec_kla,
                                            liquid_scalar->m_volumes,
                                            liquid_scalar->C.array(),
                                            gas_scalar->C.array(),
                                            &state);

    this->liquid_scalar->performStep(
        d_t, flow_liquid->transition_matrix, mat_transfer_g_liq);
    this->gas_scalar->performStep(
        d_t, flow_gas->transition_matrix, -1 * mat_transfer_g_liq);
  }
} // namespace Simulation
