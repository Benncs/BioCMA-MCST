#include "simulation/matflows.hpp"
#include <get_cumulative_proba.hpp>
#include <numeric>
#include <simulation/transport.hpp>

namespace Simulation
{
  template <typename T>
  using RowMajorDynMatrix =
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  template <typename T>
  using ColMajorDynMatrix =
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

  FlowMatrixType FlowmapToMat(std::span<double> data, size_t n_row)
  {
    if (n_row * n_row != data.size())
    {
      throw std::invalid_argument("FlowMap should be square");
    }
    Eigen::Map<ColMajorDynMatrix<double>> matrix_view(
        data.data(), static_cast<int>(n_row), static_cast<int>(n_row));

    // return matrix_view.sparseView();

    return matrix_view;
  }

  Eigen::MatrixXd get_CP(const std::vector<std::vector<size_t>> &neighbors,
                         int nb_zone,
                         const Simulation::MatFlow &flows)
  {
    // Initialize probability matrix
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(nb_zone, static_cast<int>(neighbors[0].size()));

    // Calculate cumulative sum and probability matrix
    for (int k = 0; k < nb_zone; ++k)
    {
      // Calculate cumulative sum
      double cumsum = 0;
      size_t last_n = 1e5; // Assuming a large enough initial value
      int count_nei = 0;
      for (auto &&i_neighbor : neighbors[k])
      {
        if (i_neighbor != last_n)
        {
          cumsum += flows.flows.coeff(k, (int)i_neighbor);
        }
        P.coeffRef(k, count_nei) = cumsum;
        count_nei++;
        last_n = i_neighbor;
      }

      // Normalize probabilities if cumsum is non-zero
      if (cumsum != 0)
      {
        P.row(k) /= cumsum;
      }
      else
      {
        // Set all probabilities to zero if cumsum is zero
        P.row(k).setZero();
      }
    }

    return P;
  }

  FlowMatrixType get_transition_matrix(const FlowMatrixType &flows)
  {
    auto n_compartiments = flows.rows(); // It's SHOULD be square
    FlowMatrixType m_transition(n_compartiments, n_compartiments);

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

    // m_transition.makeCompressed();
    return m_transition;
  }

  move_kernel pbf(MC::MonteCarloUnit &unit,
                  MC::ParticlesContainer & /*container*/,
                  const MatFlow *flows)
  {
    auto n_c = static_cast<int>(unit.domain.n_compartments());

    auto cumulative_probability =
        get_CP(unit.domain.getNeighbors(), n_c, *flows);

    // std::cout<<cumulative_probability<<std::endl;
    const auto &m_transition = flows->transition_matrix;

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
      const auto &i_neighbor = unit.domain.getNeighbors();
      const int max_neighbor = static_cast<int>(i_neighbor[0].size());

      const double &v_p = current_container.volume_liq;

      const double leaving_flow = m_transition.coeff(
          static_cast<int>(i_compartment), static_cast<int>(i_compartment));

      double theta_p = -v_p / leaving_flow;

      double probability = std::exp(-d_t / theta_p);

      bool leaving = false;

      if (random_number > probability)
      {
        return;
      }

      leaving = true;
      size_t next = i_neighbor[i_compartment][0];
      for (int k = 1; k < max_neighbor - 1; ++k)
      {

        auto prec = cumulative_probability.coeff(
                        static_cast<int>(i_compartment), k) < random_number2;

        auto nxt = random_number2 < cumulative_probability.coeff(
                                        static_cast<int>(i_compartment), k + 1);

        if (prec && nxt)
        {
          next = i_neighbor[i_compartment][k];
        }
      }

#pragma omp atomic
      current_container.n_cells -= static_cast<size_t>(leaving);
#pragma omp atomic
      unit.domain[next].n_cells += 1;

      particle.current_container = next;
    };

    return move_kernel;
  }

} // namespace Simulation
