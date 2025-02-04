#ifndef __SIMUALTION__PC_HYDRO_HPP__
#define __SIMUALTION__PC_HYDRO_HPP__

#include "common/kokkos_vector.hpp"
#include "simulation/alias.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cma_read/light_2d_view.hpp>
#include <cma_read/reactorstate.hpp>
#include <vector>

using FlowMatrixType = Eigen::SparseMatrix<double>;

namespace Simulation
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

  struct TransitionState
  {
    CmaRead::ReactorState state;
    PreCalculatedHydroState liquid_pc;
    PreCalculatedHydroState gas_pc;
  };

} // namespace Simulation

#endif //__SIMUALTION__PC_HYDRO_HPP__