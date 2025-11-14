#ifndef __SIMULATION_EIGEN_KOKKOS_HPP__
#define __SIMULATION_EIGEN_KOKKOS_HPP__

#include "common/traits.hpp"
#ifndef NDEBUG
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#  pragma GCC diagnostic ignored "-Wnan-infinity-disabled"
#endif
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#ifndef NDEBUG
#  pragma GCC diagnostic pop
#endif
#include <Kokkos_Core.hpp>
#include <common/common.hpp>

// constexpr auto EigenLayoutRight = Eigen::RowMajor;

// constexpr auto DataLayoutEigen = Eigen::ColMajor;
// constexpr auto CompileMatrixSizeEigen = -1;
// using MatrixType =
//     Eigen::Matrix<double, CompileMatrixSizeEigen, CompileMatrixSizeEigen,
//     DataLayoutEigen>;

// using SparseMatrixType = Eigen::SparseMatrix<double, DataLayoutEigen>;
// using DiagonalType = Eigen::DiagonalMatrix<double, CompileMatrixSizeEigen>;

// template <typename ExecSpace,typename ...Memorytrait>
// using KokkosScalarMatrix = Kokkos::View<double**, Kokkos::LayoutLeft,
// ExecSpace,Memorytrait...>;

// struct RowMajorEigenKokkos
// {

// };

// struct EigenKokkos
// {
//   KokkosScalarMatrix<HostSpace> host;
//   KokkosScalarMatrix<ComputeSpace,Kokkos::MemoryTraits<Kokkos::RandomAccess>>
//   compute; MatrixType eigen_data;

//   EigenKokkos(std::size_t n_row, std::size_t n_col);
//   [[nodiscard]] std::span<const double> get_span() const
//   {
//     return {eigen_data.data(), static_cast<size_t>(eigen_data.size())};
//   }
//   std::span<double> get_span()
//   {
//     return {eigen_data.data(), static_cast<size_t>(eigen_data.size())};
//   }

//   void update_host_to_compute() const;
//   void update_compute_to_host() const;
// };

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
using KokkosScalarMatrix =
    Kokkos::View<ftype**,
                 typename KokkosLayoutMapper<EigenLayout>::type,
                 ExecSpace,
                 MemoryTrait...>;

template <FloatingPointType ftype>
using RowMajorKokkosScalarMatrix =
    KokkosScalarMatrix<ComputeSpace, Eigen::RowMajor, ftype>;
template <FloatingPointType ftype>
using ColMajorKokkosScalarMatrix =
    KokkosScalarMatrix<ComputeSpace, Eigen::ColMajor, ftype>;

template <int EigenLayout, FloatingPointType ftype> struct EigenKokkosBase
{
  using EigenMatrix = MatrixType<EigenLayout, ftype>;
  using HostView = KokkosScalarMatrix<HostSpace, EigenLayout, ftype>;
  using ComputeView =
      KokkosScalarMatrix<ComputeSpace,
                         EigenLayout,
                         ftype,
                         Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  HostView host;
  ComputeView compute;
  EigenMatrix eigen_data;

  EigenKokkosBase(std::size_t n_row, std::size_t n_col)
      : eigen_data(MatrixType<EigenLayout, ftype>(n_row, n_col))
  {
    eigen_data.setZero();
    host = HostView(eigen_data.data(), n_row, n_col);
    compute = Kokkos::create_mirror_view_and_copy(ComputeSpace(), host);
  }

  [[nodiscard]] std::span<const double> get_span() const
  {
    return {eigen_data.data(), static_cast<size_t>(eigen_data.size())};
  }

  std::span<double> get_span()
  {
    return {eigen_data.data(), static_cast<size_t>(eigen_data.size())};
  }

  void update_host_to_compute() const
  {
    Kokkos::deep_copy(compute, host);
  }

  void update_compute_to_host() const
  {
    Kokkos::deep_copy(host, compute);
  }
};

// ColMajor version (default)
template <FloatingPointType ftype>
using EigenKokkos = EigenKokkosBase<EigenLayoutLeft, ftype>;

// RowMajor version
template <FloatingPointType ftype>
using RowMajorEigenKokkos = EigenKokkosBase<EigenLayoutRight, ftype>;

#endif
