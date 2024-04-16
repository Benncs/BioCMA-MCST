#include <cma_read/reactorstate.hpp>
#include <mc/domain.hpp>
#include <simulation/simulation.hpp>
#include <stdexcept>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <scalar_simulation.hpp>

static void mock_transfer(Eigen::ArrayXXd &vec, ReactorState *const state)
{

  double kinematic_viscosity = 1.0023e-06;
  double oxygen_diffusion_constant = 1e-9;

  double schmidtnumber = kinematic_viscosity / oxygen_diffusion_constant;

  double db = 1e-3;

  for (int i_c = 0; i_c < static_cast<int>(state->n_compartments); ++i_c)
  {
    // double eps_turb = state->energy_dissipation[i_c];
    // double kl =
    //     0.3 *
    //     std::pow(eps_turb * kinematic_viscosity, 0.25) *
    //     std::pow(schmidtnumber, -0.5);

    // double gas_fraction = domain[i_c].volume_gas /
    //              (domain[i_c].volume_gas + domain[i_c].volume_liq);
    // double a = 6 * gas_fraction / db;

    double kla = 1e-5; // kl*a;

    vec.coeffRef(1, i_c) = kla;
  }
}

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

  std::span<double>  SimulationUnit::get_contributionData()
  {
    return std::span(this->liquid_scalar->biomass_contribution.data(), this->liquid_scalar->biomass_contribution.size());
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
    }

    this->liquid_scalar->biomass_contribution.setZero();
  }

  void SimulationUnit::step(double d_t)
  {

    if (state == nullptr)
    {
      throw std::runtime_error("Error no given reactor state");
    }
    int n_species = this->liquid_scalar->n_species();
    int n_compartments = this->mc_unit->domain.n_compartments();

    Eigen::ArrayXXd vec_kla;
    vec_kla.resize(n_species, n_compartments);
    vec_kla.setZero();
    mock_transfer(vec_kla, state);

    auto c_star =
        1.3e-5 * gas_scalar->C.array() / 32e-3 * 8.314 * (273.15 + 30);

    auto transfer_g_liq = vec_kla * (c_star.array() - liquid_scalar->C.array());
    auto mat_transfer_g_liq = transfer_g_liq.matrix() * liquid_scalar->V;

    this->liquid_scalar->performStep(
        d_t, flow_liquid->transition_matrix, mat_transfer_g_liq);
    this->gas_scalar->performStep(
        d_t, flow_gas->transition_matrix, -1 * mat_transfer_g_liq);
  }
} // namespace Simulation
