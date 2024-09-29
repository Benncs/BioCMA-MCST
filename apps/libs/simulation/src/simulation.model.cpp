#include "common/common.hpp"
#include <cma_read/reactorstate.hpp>
#include <iterator>
#include <mc/domain.hpp>
#include <optional>
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

  [[nodiscard]] std::optional<std::span<double>>
  SimulationUnit::_getCgasData() const
  {
    if (!gas_scalar)
    {
      return std::nullopt;
    }
    return this->gas_scalar->getConcentrationData();
  }

  [[nodiscard]] std::tuple<size_t, size_t>
  SimulationUnit::getDim() const noexcept
  {
    return {this->liquid_scalar->n_row(), this->liquid_scalar->n_col()};
  }

  void SimulationUnit::reduceContribs(std::span<double> data,
                                      size_t n_rank) const
  {

    auto [nr, nc] = getDim();

    this->liquid_scalar->biomass_contribution.setZero();

    for (int i = 0; i < static_cast<int>(n_rank); ++i)
    {
      this->liquid_scalar->biomass_contribution.noalias() +=
          Eigen::Map<Eigen::MatrixXd>(
              &data[i * nr * nc], EIGEN_INDEX(nr), EIGEN_INDEX(nc));
    }
  }

  void SimulationUnit::clearContribution() const noexcept
  {

    this->liquid_scalar->vec_kla.setZero();

    this->liquid_scalar->biomass_contribution.setZero();
  }


  void SimulationUnit::update_feed(const double t, const double d_t,bool update_scalar)
  {
    // Get references to the index_leaving_flow and leaving_flow data members
    auto &_index_leaving_flow = this->index_leaving_flow;
    auto &_leaving_flow = this->leaving_flow;

    // Get the index of the exit compartment
    const uint64_t i_exit = mc_unit->domain.getNumberCompartments() - 1;

    // Define the set_feed lambda function
    auto set_feed = [t, d_t, i_exit, &_index_leaving_flow, &_leaving_flow,update_scalar](
                        const pimp_ptr_t &scalar,
                        auto &&descritor,
                        bool mc_f = false)
    {
      double flow = 0.; // Initialize the flow variable

      // Iterate through each current_feed in the descriptor
      for (auto &&current_feed : descritor)
      {
        current_feed.update(t, d_t); // Update the current_feed
        flow =
            current_feed.flow_value; // Get the flow_value of the current_feed

        if (update_scalar)
        {
          // Iterate through the species, positions, and values of the
          // current_feed
          for (std::size_t i_f = 0; i_f < current_feed.n_v; ++i_f)
          {
            const std::size_t i_species = current_feed.species[i_f];
            scalar->set_feed(i_species,
                             current_feed.position[i_f],
                             flow * current_feed.value[i_f]);
          }
        }
      }

      if (update_scalar)
      {
        // Set the sink for the exit compartment
        scalar->set_sink(i_exit, flow);
      }

      // Update Flow for mc particle
      if (mc_f)
      {
        _index_leaving_flow(0) = i_exit;
        _leaving_flow(0) = flow;
      }
    };

    set_feed(this->liquid_scalar, feed.liquid, true);

    if (is_two_phase_flow && feed.gas.has_value())
    {
      set_feed(this->gas_scalar, *feed.gas);
    }

    // this->liquid_scalar->set_feed(0, 0, s_feed * flow);
    // constexpr double v = 20e-3;
    // constexpr double aim = 0.5 / 3600.;
    //  constexpr double flow =aim*v;
    //  constexpr double tau = 1/(flow/v);
    // constexpr double flow =
    //     0.03813511651379644; // 1.250000e-05; // 0.035;//aim*v; //mu/2*V
    // // constexpr double tau = 1/(flow/v);
    // const uint64_t i_exit = mc_unit->domain.getNumberCompartments() - 1;
    // constexpr double s_feed = 0.25; // g/l
    // // this->liquid_scalar->set_feed(0, 0, s_feed * flow);

    // this->liquid_scalar->set_sink(i_exit, i_exit, flow);

    // if (is_two_phase_flow)
    // {
    //   // this->gas_scalar->set_feed(0, 0, 0.21 * flow);
    //   this->gas_scalar->set_sink(i_exit, i_exit, flow);
    // }

    // index_leaving_flow(0) = i_exit;
    // leaving_flow(0) = flow;
  }

  void SimulationUnit::step(double d_t,
                            const CmaRead::ReactorState &state) const
  {

    if (is_two_phase_flow)
    {
      this->liquid_scalar->mass_transfer =
          gas_liquid_mass_transfer(this->liquid_scalar->vec_kla,
                                   liquid_scalar->getVolume(),
                                   liquid_scalar->getConcentrationArray(),
                                   gas_scalar->getConcentrationArray(),
                                   state);

      this->gas_scalar->performStep(d_t,
                                    flow_gas->transition_matrix,
                                    -1 * this->liquid_scalar->mass_transfer);
    }

    this->liquid_scalar->performStep(d_t,
                                     flow_liquid->transition_matrix,
                                     this->liquid_scalar->mass_transfer);
  }
} // namespace Simulation
