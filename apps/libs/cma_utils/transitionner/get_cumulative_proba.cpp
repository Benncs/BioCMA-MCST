#ifndef NDEBUG
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif 
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#ifndef NDEBUG
#pragma GCC diagnostic pop
#endif 
#include <common/common.hpp>
#include <transitionner/transport.hpp>

namespace CmaUtils
{
  CumulativeProbaType
  get_cumulative_probabilities(CmaRead::Neighbors::Neighbors_const_view_t neighbors,
                               const FlowMatrixType& m_transition)
  {
    PROFILE_SECTION("host:get_cumulative_probabilities")
    const size_t n_compartment = neighbors.getNRow();

    // Initialize the cumulative probability matrix
    CumulativeProbaType cumsum_proba = CumulativeProbaType::Zero(
        static_cast<int>(n_compartment), static_cast<int>(neighbors.getNCol()));

    // Iterate through each compartment to compute its cumulative probabilities
    // Use of int indices to be avoid casting when using Eigen
    for (int k_compartment = 0; k_compartment < static_cast<int>(n_compartment); ++k_compartment)
    {
      double cumsum = 0.0;    // Cumulative sum for the current compartment
      int count_neighbor = 0; // Counter for the number of neighbors processed

      // Get the list of neighborsfor the current compartment
      const auto row = neighbors.getRow(k_compartment);

      // Iterate through each neighbor of the current compartment
      for (auto&& i_neighbor : row)
      {
        // Convert neighbor index to integer
        const auto colId = static_cast<int>(i_neighbor);

        // Skip if the neighbor is the same as the current compartment
        if (colId == k_compartment)
        {
          /*0.5.2: replace break by continue and ++ to allow to have neighbors ordering like:
          for compartment i [j,k,i,i,l]
          This ordering should be avoided because it's a flowmap misconception but using continue
          instead ok break fix this issue Note that continue instead of break doesn´t change normal
          behavior where ordering is [j,k,l,i,i] with i padding we just perform more iterations
          */
          count_neighbor++; // Increment the neighbor count
          continue;
        }

        const double out_flow = m_transition.coeff(k_compartment, k_compartment);
        // Calculate the transition probability from the current compartment to
        // the neighbor
        const double proba_out =
            (out_flow != 0) ? m_transition.coeff(k_compartment, colId) / std::abs(out_flow) : 0;

        // Compute the cumulative probability for the current neighbor
        const double p_cp = proba_out + cumsum;
        cumsum_proba.coeffRef(k_compartment, count_neighbor) =
            p_cp;            // Store the cumulative probability
        count_neighbor++;    // Increment the neighbor count
        cumsum += proba_out; // Update the cumulative sum
      }
    }
    return cumsum_proba; // Return the computed cumulative probability matrix
  }

} // namespace CmaUtils