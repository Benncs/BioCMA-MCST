#ifndef __GET_CUMULATIVE_PROBABILITY_HPP__
#define __GET_CUMULATIVE_PROBABILITY_HPP__

#include <Eigen/Sparse>
#include <simulation/transport.hpp>
#include <vector>

namespace Simulation
{

  static Eigen::MatrixXd
  get_CP(const std::vector<std::vector<size_t>> &neighbors,
         int n_c,
         const MatFlow &flows)
  {

    Eigen::MatrixXd probability_ij(n_c, n_c);

    for (int i = 0; i < n_c; ++i)
    {
      for (int j = 0; j < n_c; ++j)
      {
        probability_ij.coeffRef(i, j) =
            flows.flows.coeff(i, j) / (-flows.transition_matrix.coeff(i, i));
      }
    }
    auto n_neighbors = static_cast<int>(neighbors[0].size());

    Eigen::MatrixXd cumulative_propability(n_c, n_neighbors);

    for (int i_c = 0; i_c < neighbors.size(); ++i_c)
    {
      double cumsum = 0;
      for (int i_n = 0; i_n < n_neighbors; ++i_n)
      {
        cumsum +=
            probability_ij.coeff(i_c, static_cast<int>(neighbors[i_c][i_n]));
        cumulative_propability.coeffRef(i_c, i_n) = cumsum;
      }
    }

    return cumulative_propability;
  }

} // namespace Simulation
#endif //__GET_CUMULATIVE_PROBABILITY_HPP__