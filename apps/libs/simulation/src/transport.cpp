#include <Eigen/Dense>
#include <Kokkos_Core.hpp>
#include <cmt_common/zip.hpp>
#include <get_cumulative_proba.hpp>
#include <simulation/pc_hydro.hpp>
#include <transport.hpp>
#include <mc/domain.hpp>
#include <mc/particles/particles_container.hpp>
#include <mc/unit.hpp>
namespace Simulation
{


  std::vector<double> get_diag_transition(const FlowMatrixType &m_transition)
  {
    std::vector<double> res(m_transition.rows());
    for (int i = 0; i < m_transition.rows(); ++i)
    {
      res[i] = -m_transition.coeff(i, i);
    }
    return res;
  }

  // FlowMatrixType
  // get_transition_matrix(const CmaRead::FlowMap::FlowMap_const_view_t &flows)
  // {
  //   const int n_compartments =
  //       static_cast<int>(flows.getNRow()); // It SHOULD be square

  //   // Uncomment with dense matrix
  //   //  auto rd = flows.data();
  //   //  std::vector<double> flow_copy =
  //   //  std::vector<double>(rd.begin(),rd.end()); FlowMatrixType m_transition =
  //   //  Eigen::Map<Eigen::MatrixXd>(flow_copy.data(),n_compartments,n_compartments);

  //   FlowMatrixType m_transition =
  //       FlowMatrixType(n_compartments, n_compartments);

  //   for (int i = 0; i < n_compartments; ++i)
  //   {
  //     for (int j = 0; j < n_compartments; ++j)
  //     {
  //       if(i!=j)
  //       {
  //                 const double val = flows(i, j);
  //       m_transition.coeffRef(i, j) += val;
  //       }

  //     }
  //   }

  //   for (int i = 0; i < n_compartments; ++i)
  //   {
  //     m_transition.coeffRef(i, i) = -1. * m_transition.row(i).sum();
  //   }
  //   m_transition.makeCompressed();

  //   return m_transition;
  // }

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

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(n_compartments * n_compartments); // Reserve space for all elements

    // Temporary vector to keep track of the row sums for diagonal elements
    std::vector<double> rowSums(n_compartments, 0.0);

    for (int i = 0; i < n_compartments; ++i)
    {
      for (int j = 0; j < n_compartments; ++j)
      {
        if (i != j)
        {
          const double val = flows(i, j);
          tripletList.emplace_back(i, j, val);
          rowSums[i] += val; // Accumulate the row sum for the diagonal element
        }
      }
    }

    // Set the diagonal elements to the negative row sums
    for (int i = 0; i < n_compartments; ++i)
    {
      tripletList.emplace_back(i, i, -rowSums[i]);
    }

    // Construct the sparse matrix from the triplets
    m_transition.setFromTriplets(tripletList.begin(), tripletList.end());
    return m_transition;
  }

} // namespace Simulation
