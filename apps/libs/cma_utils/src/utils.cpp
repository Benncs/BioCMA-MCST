#include "Kokkos_Array.hpp"
#include "Kokkos_Assert.hpp"
#include "common/common.hpp"
#include <cma_utils/d_transitionner.hpp>

#define UNROLL_CHECK(i, s)                                                     \
  if (i < s)                                                                   \
  {                                                                            \
    triplets[i] = T(rows[i], cols[i], vals[i]);                                \
  }

#define UNROLL(idx) triplets[idx] = T(rows[idx], cols[idx], vals[idx]);
#include <iostream>
namespace CmaUtils
{
  void sparse_from_coo(Eigen::SparseMatrix<double>& res,
                       std::size_t n_c,
                       std::span<const size_t> rows,
                       std::span<const size_t> cols,
                       std::span<const double> vals)
  {

    // KOKKOS_ASSERT(res.cols() == EIGEN_INDEX(n_c));
    // using T = Eigen::Triplet<double>;
    // std::vector<T> triplets(vals.size());

    // triplets.reserve(vals.size());

    // for (std::size_t i = 0; i < vals.size(); ++i)
    // {
    //   triplets[i] = T(rows[i], cols[i], vals[i]);
    // }

    // // res.
    // res.setFromTriplets(triplets.begin(), triplets.end());
    // res.makeCompressed();

    KOKKOS_ASSERT(res.cols() == EIGEN_INDEX(n_c));

    using T = Eigen::Triplet<double>;
    const auto n = vals.size();
    std::vector<T> triplets(n);
    triplets.reserve(n);

    // #pragma omp simd
    for (std::size_t i = 0; i < n; i += 8)
    {
      std::size_t idx0 = i;
      std::size_t idx1 = i + 1;
      std::size_t idx2 = i + 2;
      std::size_t idx3 = i + 3;
      std::size_t idx4 = i + 4;
      std::size_t idx5 = i + 5;
      std::size_t idx6 = i + 6;
      std::size_t idx7 = i + 7;

      if (idx7 < n)
      {
        UNROLL(idx0);
        UNROLL(idx1);
        UNROLL(idx2);
        UNROLL(idx3);
        UNROLL(idx4);
        UNROLL(idx5);
        UNROLL(idx6);
        UNROLL(idx7);
      }
      else
      {
        UNROLL_CHECK(idx0, n)
        UNROLL_CHECK(idx1, n)
        UNROLL_CHECK(idx2, n)
        UNROLL_CHECK(idx3, n)
        UNROLL_CHECK(idx4, n)
        UNROLL_CHECK(idx5, n)
        UNROLL_CHECK(idx6, n)
      }
    }

    // After filling the triplets, create the sparse matrix
    res.setFromTriplets(triplets.begin(), triplets.end());
    res.makeCompressed();
  }

  DiagonalView<ComputeSpace>
  diagonal_from_sparse(const Eigen::SparseMatrix<double>& transition)
  {
    DiagonalView<HostSpace> diag("host_diag", transition.rows());

    for (int i = 0; i < transition.rows(); ++i)
    {
      diag[i] = -transition.coeff(i, i);
    }
    return Kokkos::create_mirror_view_and_copy(ComputeSpace(), diag);
    ;
  }
} // namespace CmaUtils