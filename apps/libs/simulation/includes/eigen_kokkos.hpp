#ifndef __SIMULATION_EIGEN_KOKKOS_HPP__
#define __SIMULATION_EIGEN_KOKKOS_HPP__

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Kokkos_Core.hpp>
#include <common/common.hpp>
#include <common/traits.hpp>
#include <span>

// Layout definitions
constexpr auto EigenLayoutRight = Eigen::RowMajor;
constexpr auto EigenLayoutLeft = Eigen::ColMajor;
constexpr auto CompileMatrixSizeEigen = -1;

// Templated Eigen Matrix Type
template <int Layout, FloatingPointType ftype>
using MatrixType = Eigen::
    Matrix<ftype, CompileMatrixSizeEigen, CompileMatrixSizeEigen, Layout>;

template <FloatingPointType ftype>
using ColMajorMatrixtype = MatrixType<EigenLayoutLeft, ftype>;

template <int Layout, FloatingPointType ftype>
using SparseMatrixType = Eigen::SparseMatrix<ftype, Layout>;

using DiagonalType = Eigen::DiagonalMatrix<double, CompileMatrixSizeEigen>;

// template <typename ExecSpace, typename... MemoryTraits>
// using KokkosScalarMatrix = Kokkos::View<double**, Kokkos::LayoutLeft,
// ExecSpace, MemoryTraits...>;

template <int EigenLayout> struct KokkosLayoutMapper;

template <> struct KokkosLayoutMapper<Eigen::ColMajor>
{
  using type = Kokkos::LayoutLeft;
};

template <> struct KokkosLayoutMapper<Eigen::RowMajor>
{
  using type = Kokkos::LayoutRight;
};

template <typename ExecSpace,
          int EigenLayout,
          FloatingPointType ftype,
          typename... MemoryTrait>
using KokkosScalarMatrix
    = Kokkos::View<ftype**,
                   typename KokkosLayoutMapper<EigenLayout>::type,
                   ExecSpace,
                   MemoryTrait...>;

template <FloatingPointType ftype>
using RowMajorKokkosScalarMatrix
    = KokkosScalarMatrix<ComputeSpace, Eigen::RowMajor, ftype>;
template <FloatingPointType ftype>
using ColMajorKokkosScalarMatrix
    = KokkosScalarMatrix<ComputeSpace, Eigen::ColMajor, ftype>;

template <int EigenLayout, FloatingPointType float_type> struct EigenKokkosBase
{
  using float_t = float_type;
  using EigenMatrix = MatrixType<EigenLayout, float_type>;
  using HostView = KokkosScalarMatrix<HostSpace, EigenLayout, float_type>;
  using ComputeView
      = KokkosScalarMatrix<ComputeSpace,
                           EigenLayout,
                           float_type,
                           Kokkos::MemoryTraits<Kokkos::RandomAccess> >;

  HostView host;
  ComputeView compute;
  EigenMatrix eigen_data;

  EigenKokkosBase(std::size_t n_row, std::size_t n_col)
      : eigen_data(MatrixType<EigenLayout, float_type>(n_row, n_col))
  {
    eigen_data.setZero();
    host = HostView(eigen_data.data(), n_row, n_col);
    compute = Kokkos::create_mirror_view_and_copy(ComputeSpace(), host);
  }

  [[nodiscard]] std::span<const double>
  get_span() const
  {
    return { eigen_data.data(), static_cast<size_t>(eigen_data.size()) };
  }

  std::span<double>
  get_span()
  {
    return { eigen_data.data(), static_cast<size_t>(eigen_data.size()) };
  }

  void
  update_host_to_compute() const
  {
    Kokkos::deep_copy(compute, host);
  }

  void
  update_compute_to_host() const
  {
    Kokkos::deep_copy(host, compute);
  }
};

// ColMajor version (default)
template <FloatingPointType floatype>
using EigenKokkos = EigenKokkosBase<EigenLayoutLeft, floatype>;

// RowMajor version
template <FloatingPointType floatype>
using RowMajorEigenKokkos = EigenKokkosBase<EigenLayoutRight, floatype>;

#endif
