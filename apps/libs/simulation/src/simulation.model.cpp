#include "mc/domain.hpp"
#include "reactorstate.hpp"
#include <simulation/simulation.hpp>
#include <stdexcept>

#include <iostream>

static void mock_transfer(const MC::ReactorDomain &domain,
                          Eigen::ArrayXXd &vec,
                          ReactorState *const state)
{

//https://docs.mstarcfd.com/1b_HowToGuides/predicting-mass-transfer-kLA.html
  double kinematic_viscosity = 1.0023e-06;
  double oxygen_diffusion_constant = 1e-9;

  double schmidtnumber = kinematic_viscosity / oxygen_diffusion_constant;

  double db = 1e-3;

  for (int i_c = 0; i_c < static_cast<int>(state->n_compartments); ++i_c)
  {

    double kl =
        0.3 *
        std::pow(state->energy_dissipation[i_c] * kinematic_viscosity, 0.25) *
        std::pow(schmidtnumber, -0.5);

    double eps = domain[i_c].volume_gas /
                 (domain[i_c].volume_gas + domain[i_c].volume_liq);
    double a = 6 * eps / db;
    vec.coeffRef(1, i_c) = kl * a;
  }
}

namespace Simulation
{
  void SimulationUnit::clear_contribution()
  {
    for (auto &&i : contribs)
    {
      i.setZero();
    }
    this->liquid_scalar->biomass_contribution.setZero();

    // FIXME
    for (auto i = 0; i < liquid_scalar->C.cols(); ++i)
    {
      liquid_scalar->C.coeffRef(0, i) = 10;
      
    }
    if(gas_scalar)
        gas_scalar->C.coeffRef(1, 10) = 3;
  }

  void SimulationUnit::step(double d_t)
  {

    if (!state)
    {
      throw std::runtime_error("Error no given reactor state");
    }
    int n_species = this->liquid_scalar->n_species();
    int n_compartments = this->unit->domain.n_compartments();

    Eigen::ArrayXXd vec_kla;
    vec_kla.resize(n_species, n_compartments);
    vec_kla.setZero();
    mock_transfer(unit->domain, vec_kla, state);


    auto c_star = 1.3e-5*gas_scalar->C.array()/32e-3*8.314*(273.15+30);



    auto transfer_g_liq = vec_kla * (c_star.array() - liquid_scalar->C.array());
    auto mat_transfer_g_liq = transfer_g_liq.matrix()* liquid_scalar->V;

    this->liquid_scalar->performStep(
        d_t, flow_liquid.transition_matrix, mat_transfer_g_liq);
    this->gas_scalar->performStep(
        d_t, flow_gas.transition_matrix, -1 * mat_transfer_g_liq);
  }
} // namespace Simulation
