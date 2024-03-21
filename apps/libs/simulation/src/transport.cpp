#include <get_cumulative_proba.hpp>
#include <simulation/transport.hpp>

#include <iostream>
namespace Simulation
{

  Eigen::MatrixXd get_CP(const std::vector<std::vector<size_t>> &neighbors,
                         int n_c,
                         const Simulation::MatFlow &flows)
  {

    Eigen::MatrixXd probability_ij(n_c, n_c);

    for (int i = 0; i < n_c; ++i)
    {
      for (int j = 0; j < n_c; ++j)
      {
        if (flows.transition_matrix.coeff(i, i) == 0)
        {
          probability_ij(i, j) = 0.0; 
        }
        else
        {
          probability_ij.coeffRef(i, j) =
              flows.flows.coeff(i, j) / (-flows.transition_matrix.coeff(i, i));
        }
      }
    }
    auto n_neighbors = static_cast<int>(neighbors[0].size());

    Eigen::MatrixXd cumulative_probability(n_c, n_neighbors);

    for (int i_c = 0; i_c < static_cast<int>(neighbors.size()); ++i_c)
    {
      double cumsum = 0;
      for (int i_n = 0; i_n < n_neighbors; ++i_n)
      {
        cumsum +=
            probability_ij.coeff(i_c, static_cast<int>(neighbors[i_c][i_n]));
        cumulative_probability.coeffRef(i_c, i_n) = cumsum;
      } 
      //TO CHECK
       if (cumsum != 0) {
            cumulative_probability.row(i_c) /= cumsum;
        }
    }

    return cumulative_probability;
  }

  Eigen::SparseMatrix<double>
  get_transition_matrix(const Eigen::SparseMatrix<double> &flows)
  {
    auto n_compartiments = flows.rows(); // It's SHOULD be square
    Eigen::SparseMatrix<double> m_transition(n_compartiments, n_compartiments);

    // TODO OPTI
    for (int i = 0; i < n_compartiments; ++i)
    {
      for (int j = 0; j < n_compartiments; ++j)
      {
        if (i != j)
        {
          m_transition.coeffRef(i, j) = flows.coeff(i, j);
        }
        else
        {
          Eigen::RowVectorXd row = flows.row(i);

          auto sum_i =
              std::accumulate(row.begin(), row.end(), 0.) - flows.coeff(i, i);

          m_transition.coeffRef(i, i) = -sum_i;
        }
      }
    }
    m_transition.makeCompressed();
    return m_transition;
  }

  move_kernel pbf(MC::MonteCarloUnit &unit,
                  MC::ParticlesContainer &container,
                  const MatFlow &flows)
  {
    auto n_c = static_cast<int>(unit.domain.n_compartments());

    auto cumulative_probability =
        get_CP(unit.domain.getNeighbors(), n_c, flows);


    auto &m_transition = flows.transition_matrix;

    auto move_kernel = [&unit,
                        &m_transition,
                        cumulative_probability = std::move(
                            cumulative_probability)](double random_number,
                                                     double random_number2,
                                                     MC::Particles &particle,
                                                     double d_t) -> void
    {
      const size_t i_compartment = particle.current_container;
      auto &current_container = unit.domain[i_compartment];
      const auto i_neighbor = unit.domain.getNeighbors();
      const int max_neighbor = static_cast<int>(i_neighbor[0].size());

      const double &v_p = current_container.volume;

      const double leaving_flow = m_transition.coeff(
          static_cast<int>(i_compartment), static_cast<int>(i_compartment));

      double theta_p = -v_p / leaving_flow;

      double probability = std::exp(-d_t / theta_p);

      bool leaving = false;

      if (random_number > probability)
      {
        return;
      }
      else
      {
        leaving = true;
        size_t next = i_neighbor[i_compartment][0];
        for (int k = 1; k < max_neighbor - 1; ++k)
        {

          auto prec = cumulative_probability.coeff(
                          static_cast<int>(i_compartment), k) < random_number2;

          auto nxt =
              random_number2 < cumulative_probability.coeff(
                                   static_cast<int>(i_compartment), k + 1);

          if (prec && nxt)
          {
            next = i_neighbor[i_compartment][k];
          }
        }

#pragma omp atomic
        current_container.n_cells -= static_cast<size_t>(leaving);

        particle.current_container = next;
#pragma omp atomic
        unit.domain[next].n_cells += 1;
      }
    };

    return move_kernel;
  }

} // namespace Simulation
