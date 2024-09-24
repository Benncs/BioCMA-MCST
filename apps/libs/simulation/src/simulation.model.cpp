#include <cma_read/reactorstate.hpp>
#include <iterator>
#include <mc/domain.hpp>
#include <simulation/simulation.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <hydro/mass_transfer.hpp>
#include <iterator>
#include <scalar_simulation.hpp>
#include <stdexcept>

namespace Simulation
{

  std::span<double> SimulationUnit::getContributionData() const
  {

    return liquid_scalar->getContributionData();
  }

  std::span<double> SimulationUnit::getCliqData() const
  {
    return this->liquid_scalar->getConcentrationData();
  }

  [[nodiscard]] std::span<double> SimulationUnit::getCgasData() const
  {
    if (!gas_scalar)
    {
      return {};
    }
    return this->gas_scalar->getConcentrationData();
  }

  [[nodiscard]] std::tuple<size_t, size_t>
  SimulationUnit::getDim() const noexcept
  {
    return {this->liquid_scalar->concentration.rows(),
            this->liquid_scalar->concentration.cols()};
  }

  void SimulationUnit::reduceContribs(std::span<double> data,
                                      size_t n_rank) const
  {

    size_t nr = this->liquid_scalar->concentration.rows();
    size_t nc = this->liquid_scalar->concentration.cols();

    this->liquid_scalar->biomass_contribution.setZero();

    for (int i = 0; i < static_cast<int>(n_rank); ++i)
    {
      this->liquid_scalar->biomass_contribution.noalias() +=
          Eigen::Map<Eigen::MatrixXd>(
              &data[i * nr * nc], static_cast<int>(nr), static_cast<int>(nc));
    }
  }

  void SimulationUnit::clearContribution() const noexcept
  {

    this->liquid_scalar->vec_kla.setZero();

    this->liquid_scalar->biomass_contribution.setZero();
  }

  void SimulationUnit::update_feed(const double t,const double d_t) const
  {
    constexpr double v = 0.09;
    constexpr double aim = 0.5/3600.;
    constexpr double flow =0.;// 0.035;//aim*v; //mu/2*V
    // constexpr double tau = 1/(flow/v);
    constexpr uint64_t i_exit = 0;
    constexpr double s_feed = 1.;//0.25; //g/l
   
    this->liquid_scalar->feed.coeffRef(0, 0) = s_feed * flow;
    
    
    this->liquid_scalar->sink.coeffRef(i_exit,i_exit) = flow;
    // for (int i = 1; i < this->liquid_scalar->concentration.cols() - 2; ++i)
    // {
    //   this->liquid_scalar->feed.coeffRef(0, i) = 50 * 10 / 3600;
    // }

    index_leaving_flow(0) = i_exit;
    leaving_flow(0) = flow;
  }

  void SimulationUnit::step(double d_t,
                            const CmaRead::ReactorState &state) const
  {

    auto mat_transfer_g_liq =
        (is_two_phase_flow)
            ? gas_liquid_mass_transfer(this->liquid_scalar->vec_kla,
                                       liquid_scalar->getVolume(),
                                       liquid_scalar->concentration.array(),
                                       gas_scalar->concentration.array(),
                                       state)
            : Eigen::MatrixXd::Zero(this->liquid_scalar->concentration.rows(),
                                    this->liquid_scalar->concentration.cols());

    this->liquid_scalar->performStep(
        d_t, flow_liquid->transition_matrix, mat_transfer_g_liq);

    if (is_two_phase_flow)
    {

      this->gas_scalar->performStep(
          d_t, flow_gas->transition_matrix, -1 * mat_transfer_g_liq);
    }
  }
} // namespace Simulation
