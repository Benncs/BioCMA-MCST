#ifndef __CACHE_HYDRO_STATE__
#define __CACHE_HYDRO_STATE__

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <common/kokkos_vector.hpp>
#include <vector>

template <typename ExecSpace>
using DiagonalView = Kokkos::
    View<double*, Kokkos::LayoutLeft, ExecSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

template <typename Space>
using CumulativeProbabilityView =
    Kokkos::View<double**, Kokkos::LayoutRight, Space, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

using FlowMatrixType = Eigen::SparseMatrix<double>;

namespace CmaUtils
{
  class ProxyPreCalculatedHydroState;

  class PreCalculatedHydroState
  {
  public:
    friend class ProxyCalculatedHydroState;

    PreCalculatedHydroState();
    ~PreCalculatedHydroState() = default;
    explicit PreCalculatedHydroState(const FlowMatrixType& _tm);
    PreCalculatedHydroState& operator=(PreCalculatedHydroState&& rhs) = default;
    PreCalculatedHydroState& operator=(const PreCalculatedHydroState& rhs) = delete;
    PreCalculatedHydroState(PreCalculatedHydroState&& other) = default;
    PreCalculatedHydroState(const PreCalculatedHydroState& other) = delete;

    FlowMatrixType transition_matrix;
    Eigen::Matrix<double, -1, -1, Eigen::RowMajor> cumulative_probability;
    std::vector<double> inverse_volume;
    std::vector<double> volume;
    DiagonalView<ComputeSpace> diagonal_compute;
    CumulativeProbabilityView<ComputeSpace> compute_cumulative_probability;


    [[nodiscard]] const FlowMatrixType& get_transition() const;
    [[nodiscard]] DiagonalView<ComputeSpace> get_kernel_diagonal() const;
  

  private:
    
  };

  inline DiagonalView<ComputeSpace> PreCalculatedHydroState::get_kernel_diagonal() const
  {
    return diagonal_compute;
  }



} // namespace CmaUtils

#endif //__CACHE_HYDRO_STATE__