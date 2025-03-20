#ifndef __SIMULATION_EIGEN_KOKKOS_HPP__
#define __SIMULATION_EIGEN_KOKKOS_HPP__

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Kokkos_Core.hpp>
#include <common/common.hpp>
#include <common/kokkos_vector.hpp>

constexpr auto DataLayoutEigen = Eigen::ColMajor;
constexpr auto CompileMatrixSizeEigen = -1;
using MatrixType =
    Eigen::Matrix<double, CompileMatrixSizeEigen, CompileMatrixSizeEigen, DataLayoutEigen>;
using SparseMatrixType = Eigen::SparseMatrix<double, DataLayoutEigen>;
using DiagonalType = Eigen::DiagonalMatrix<double, CompileMatrixSizeEigen>;

template <typename ExecSpace,typename ...Memorytrait>
using KokkosScalarMatrix = Kokkos::View<double**, Kokkos::LayoutLeft, ExecSpace,Memorytrait...>;

struct EigenKokkos
{
  KokkosScalarMatrix<HostSpace> host;
  KokkosScalarMatrix<ComputeSpace,Kokkos::MemoryTraits<Kokkos::RandomAccess>> compute;
  MatrixType eigen_data;

  EigenKokkos(std::size_t n_row, std::size_t n_col);
  [[nodiscard]] std::span<const double> get_span() const
  {
    return {eigen_data.data(), static_cast<size_t>(eigen_data.size())};
  }
  std::span<double> get_span()
  {
    return {eigen_data.data(), static_cast<size_t>(eigen_data.size())};
  }

  void update_host_to_compute() const;
  void update_compute_to_host() const;
};


#endif 