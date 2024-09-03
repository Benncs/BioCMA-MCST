#include <get_cumulative_proba.hpp>

namespace Simulation
{
  Eigen::MatrixXd get_cumulative_probabilities(
      CmaRead::Neighbors::Neighbors_const_view_t neighbors,
      const FlowMatrixType &m_transition)
  {
    const size_t n_compartment = neighbors.getNRow();

    // Initialize the cumulative probability matrix
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(
        static_cast<int>(n_compartment), static_cast<int>(neighbors.getNCol()));

    // Iterate through each compartment to compute its cumulative probabilities
    // Use of int indices to be avoid casting when using Eigen
    for (int k_compartment = 0; k_compartment < static_cast<int>(n_compartment);
         ++k_compartment)
    {
      double cumsum = 0.0;    // Cumulative sum for the current compartment
      int count_neighbor = 0; // Counter for the number of neighbors processed

      // Get the list of neighborsfor the current compartment
      const auto row = neighbors.getRow(k_compartment);

      // Iterate through each neighbor of the current compartment
      for (auto &&i_neighbor : row)
      {
        // Convert neighbor index to integer
        const auto colId = static_cast<int>(i_neighbor);

        // Skip if the neighbor is the same as the current compartment
        if (colId == k_compartment)
        {
          break;
        }

        // Calculate the transition probability from the current compartment to
        // the neighbor
        const double proba_out =
            m_transition.coeff(k_compartment, colId) /
            std::abs(m_transition.coeff(k_compartment, k_compartment));

        // Compute the cumulative probability for the current neighbor
        const double p_cp = proba_out + cumsum;
        P.coeffRef(k_compartment, count_neighbor) =
            p_cp;            // Store the cumulative probability
        count_neighbor++;    // Increment the neighbor count
        cumsum += proba_out; // Update the cumulative sum
      }
    }

    return P; // Return the computed cumulative probability matrix
  }
} // namespace Simulation