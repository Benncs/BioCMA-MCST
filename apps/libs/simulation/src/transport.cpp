#include "common/thread_safe_container.hpp"
#include "mc/domain.hpp"
#include "simulation/matflows.hpp"
#include <climits>
#include <get_cumulative_proba.hpp>
#include <numeric>
#include <transport.hpp>

// TODO REMOVE
#include <iostream>

//TODO REMOVE 
template <typename T>
using RowMajorDynMatrix =
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <typename T>
using ColMajorDynMatrix =
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

namespace Simulation
{

  static size_t
  find_next_compartment(size_t i_compartment,
                        double random_number,
                        std::span<const size_t> i_neighbor,
                        const Eigen::MatrixXd &cumulative_probability);

  move_kernel population_balance_flow(MC::ReactorDomain &domain,
            
                                      const MatFlow *flows)
  {
    const int n_c = static_cast<int>(domain.getNumberCompartments());

    auto __cumulative_probability = get_CP(domain.getNeighbors(), n_c, *flows);

    const auto &m_transition = flows->transition_matrix;

    auto move_kernel =
        [&m_transition,
         cumulative_probability =
             std::move(__cumulative_probability)](double random_number,
                                                  double random_number2,
                                                  MC::ReactorDomain &domain,
                                                  MC::Particles &particle,
                                                  double d_t) -> void
    {
      const size_t i_compartment = particle.current_container;
      auto &current_container = domain[i_compartment];
      const std::span<const size_t> i_neighbor =
          domain.getNeighbors(i_compartment);

      const double &v_p = current_container.volume_liq;

      const double leaving_flow = m_transition.coeff(
          static_cast<int>(i_compartment), static_cast<int>(i_compartment));

      const double theta_p = -v_p / leaving_flow;

      const double probability = std::exp(-d_t / theta_p);

      if (random_number > probability)
      {
        return;
      }

      size_t next = find_next_compartment(
          i_compartment, random_number2, i_neighbor, cumulative_probability);

      __ATOM_DECR__(current_container.n_cells) // Cell leaves current
                                               // compartment
      __ATOM_INCR__(domain[next].n_cells); // Cell go to new compartment

      particle.current_container = next;
    };

    return move_kernel;
  }

  Eigen::MatrixXd get_CP(const std::vector<std::vector<size_t>> &neighbors,
                         int nb_zone,
                         const Simulation::MatFlow &flows)
  {
    // Initialize probability matrix
    static Eigen::MatrixXd P =
        Eigen::MatrixXd::Zero(nb_zone, static_cast<int>(neighbors[0].size()));

    // Calculate cumulative sum and probability matrix
    for (int k = 0; k < nb_zone; ++k)
    {
      // Calculate cumulative sum
      double cumsum = 0;
      size_t last_n = UINT_MAX; // Assuming a large enough initial value
      int count_nei = 0;

      for (auto &&i_neighbor : neighbors[k])
      {
        if (i_neighbor != last_n)
        {
          cumsum += flows.flows.coeff(k, static_cast<int>(i_neighbor));
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

  size_t find_next_compartment(size_t i_compartment,
                               double random_number,
                               std::span<const size_t> i_neighbor,
                               const Eigen::MatrixXd &cumulative_probability)
  {
    const int max_neighbor = static_cast<int>(i_neighbor.size());
    size_t next = i_neighbor[0];
    for (int k = 1; k < max_neighbor - 1; ++k)
    {

      auto prec = cumulative_probability.coeff(static_cast<int>(i_compartment),
                                               k) < random_number;

      auto nxt = random_number < cumulative_probability.coeff(
                                     static_cast<int>(i_compartment), k + 1);

      if (prec && nxt)
      {
        next = i_neighbor[k];
      }
    }
    return next;
  }

  FlowMatrixType flowmap_to_matrix(std::span<double> data, size_t n_row)
  {
    if (n_row * n_row != data.size())
    {
      throw std::invalid_argument("FlowMap should be square");
    }
    Eigen::Map<ColMajorDynMatrix<double>> matrix_view(
        data.data(), static_cast<int>(n_row), static_cast<int>(n_row));

    return matrix_view;
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
          // Eigen::RowVectorXd col = flows.col(i);

          // auto sum_i =
          //     std::accumulate(col.begin(), col.end(), 0.) - flows.coeff(i,
          //     i);

          // m_transition.coeffRef(i, i) = -sum_i;

          Eigen::VectorXd row = flows.row(i); // Get the i-th row

          auto sum_i =
              row.sum() - flows(i, i); // Sum all elements in the row, excluding
                                       // the diagonal element

          m_transition.coeffRef(i, i) =
              -sum_i; // Set the diagonal element of the transition matrix
        }
      }
    }

    // m_transition.makeCompressed();
    return m_transition;
  }

} // namespace Simulation
