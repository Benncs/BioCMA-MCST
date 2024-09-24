#include <Eigen/Dense>
#include <Kokkos_Core.hpp>
#include <cmt_common/zip.hpp>
#include <get_cumulative_proba.hpp>
#include <simulation/pc_hydro.hpp>
#include <transport.hpp>

namespace Simulation
{

  size_t
  find_next_compartment(int i_compartment,
                        double random_number,
                        std::span<const size_t> i_neighbor,
                        CmaRead::L2DView<const double> cumulative_probability)
  {
    const int max_neighbor = static_cast<int>(i_neighbor.size());
    size_t next = i_neighbor[0]; // Default to the first neighbor

    // Iterate through the neighbors to find the appropriate next compartment
    for (int k_neighbor = 0; k_neighbor < max_neighbor - 1; ++k_neighbor)
    {
      // Get the cumulative probability range for the current neighbor
      const auto pi = cumulative_probability(i_compartment, k_neighbor);
      const auto pn = cumulative_probability(i_compartment, k_neighbor + 1);

      // Check if the random number falls within the probability range
      if (random_number <= pn && pi <= random_number)
      {
        next = i_neighbor[k_neighbor + 1]; // Update to the next neighbor
        // No need to break, as we're looking for the last valid neighbor in the
        // range
      }
    }

    return next; // Return the index of the chosen next compartment
  }

 
  std::vector<double> get_diag_transition(const FlowMatrixType &m_transition)
  {
    std::vector<double> res(m_transition.rows());
    for (int i = 0; i < m_transition.rows(); ++i)
    {
      res[i] = -m_transition.coeff(i, i);
    }
    return res;
  }

  FlowMatrixType
  get_transition_matrix(const CmaRead::FlowMap::FlowMap_const_view_t &flows)
  {
    const int n_compartments = 
        static_cast<int>(flows.getNRow()); // It SHOULD be square

    // Uncomment with dense matrix
    //  auto rd = flows.data();
    //  std::vector<double> flow_copy =
    //  std::vector<double>(rd.begin(),rd.end()); FlowMatrixType m_transition =
    //  Eigen::Map<Eigen::MatrixXd>(flow_copy.data(),n_compartments,n_compartments);

    FlowMatrixType m_transition =
        FlowMatrixType(n_compartments, n_compartments);

    for (int i = 0; i < n_compartments; ++i)
    {
      for (int j = 0; j < n_compartments; ++j)
      {
        if(i!=j)
        {
                  const double val = flows(i, j);
        m_transition.coeffRef(i, j) += val;  
        }

      }
    }

    for (int i = 0; i < n_compartments; ++i)
    {
      m_transition.coeffRef(i, i) = -1. * m_transition.row(i).sum();
    }
    m_transition.makeCompressed();

    return m_transition;
  }

} // namespace Simulation
