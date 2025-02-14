#ifndef __CACHE_HYDRO_STATE__
#define __CACHE_HYDRO_STATE__

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cma_read/light_2d_view.hpp>
#include <cma_read/reactorstate.hpp>
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

  class PreCalculatedHydroState
  {
  public:
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

    void set_transition_matrix(const CmaRead::FlowMap::FlowMap_const_view_t& flows_view);
    void set_transition_matrix(FlowMatrixType&& matrix);
    void set_cumulative_probability(const CmaRead::Neighbors::Neighbors_const_view_t& neighbors);
    void set_diag_transition(std::vector<double>&& diag);

    [[nodiscard]] const FlowMatrixType& get_transition() const;
    DiagonalView<ComputeSpace> get_kernel_diagonal();

  private:
    DiagonalView<ComputeSpace> diagonal_compute;
    CumulativeProbabilityView<ComputeSpace> compute_cumulative_probability;
  };

  inline DiagonalView<ComputeSpace> PreCalculatedHydroState::get_kernel_diagonal()
  {
    return diagonal_compute;
  }

  

} // namespace CmaUtils

#endif //__CACHE_HYDRO_STATE__