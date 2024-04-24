#include <cma_read/reactorstate.hpp>
#include <mc/domain.hpp>
#include <simulation/simulation.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <hydro/mass_transfer.hpp>
#include <scalar_simulation.hpp>

// double kla = 1e-5; //MOCK USING WRONG FLOWMAP

namespace Simulation
{

  void SimulationUnit::post_process_reducing()
  {
    for (size_t i_thread = 0; i_thread < n_thread; ++i_thread)
    {
      this->liquid_scalar->merge(i_thread);
      this->mc_container->merge(i_thread);
    }
  }

  std::span<double> SimulationUnit::getContributionData()
  {

    return liquid_scalar->getContributionData();
  }

  std::span<double> SimulationUnit::getCliqData()
  {
    return this->liquid_scalar->getCliqData();
  }

  void SimulationUnit::reduceContribs(std::span<double> data, size_t n_rank)
  {

    size_t nr = this->liquid_scalar->C.rows();
    size_t nc = this->liquid_scalar->C.cols();

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
      this->liquid_scalar->getThreadContribs()[i].setZero();
      this->liquid_scalar->vec_kla.setZero();
    }

    this->liquid_scalar->biomass_contribution.setZero();
  }

  void SimulationUnit::step(double d_t, ReactorState &state)
  {

    // if (state == nullptr)
    // {
    //   throw std::runtime_error("Error no given reactor state");
    // }

    auto mat_transfer_g_liq =
        gas_liquid_mass_transfer(this->liquid_scalar->vec_kla,
                                 liquid_scalar->getVolume(),
                                 liquid_scalar->C.array(),
                                 gas_scalar->C.array(),
                                 state);

    this->liquid_scalar->performStep(
        d_t, flow_liquid->transition_matrix, mat_transfer_g_liq);
    this->gas_scalar->performStep(
        d_t, flow_gas->transition_matrix, -1 * mat_transfer_g_liq);
  }
} // namespace Simulation
