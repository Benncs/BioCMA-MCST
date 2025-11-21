#ifndef __CMA_UTILS_D_TRANSITION_HPP__
#define __CMA_UTILS_D_TRANSITION_HPP__

#include "Kokkos_Core.hpp"
#include "cma_utils/cache_hydro_state.hpp"
#include "common/common.hpp"
#include <cma_utils/iteration_state.hpp>
#include <lib.rs.h>
namespace CmaUtils
{
  using TransitionnerPtrType = rust::Box<TransionnerWrapper>;
  using IterationStatePtrType = ::rust::Box<::IterationStateWrapper>;

  using StateCooMatrixType = ::rust::Box<::CooMatrixWrap>;

  void sparse_from_coo(Eigen::SparseMatrix<double>& res,
                       std::size_t n_c,
                       std::span<const size_t> rows,
                       std::span<const size_t> cols,
                       std::span<const double> vals);

  DiagonalView<ComputeSpace>
  diagonal_from_sparse(const Eigen::SparseMatrix<double>& transition);

  // CumulativeProbabilityView<ComputeSpace>
  // get_cumulative_probabilities(HostNeighsView neighbors,
  //                              const FlowMatrixType& m_transition);

} // namespace CmaUtils

#endif