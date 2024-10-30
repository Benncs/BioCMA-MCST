#ifndef __SIMUALTION__PC_HYDRO_HPP__
#define __SIMUALTION__PC_HYDRO_HPP__

#include <cma_read/light_2d_view.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cma_read/reactorstate.hpp>

// using FlowMatrixType = Eigen::MatrixXd;

using FlowMatrixType = Eigen::SparseMatrix<double>;


//TODO MOVE ELSEWHERE
inline CmaRead::L2DView<double> get_eigen_view(Eigen::MatrixXd& matrix)
{
  return CmaRead::L2DView<double>(
        std::span<double>(matrix.data(),
                                matrix.size()),
        matrix.rows(),
        matrix.cols(),
        false);
}

inline CmaRead::L2DView<const double> get_eigen_view(const Eigen::MatrixXd& matrix)
{
  return CmaRead::L2DView<const double>(
        std::span<const double>(matrix.data(),
                                matrix.size()),
        matrix.rows(),
        matrix.cols(),
        false);
}


namespace Simulation
{

  class PreCalculatedHydroState
  {
  public:
    // FlowMatrixType flows;
    FlowMatrixType transition_matrix;
    Eigen::Matrix<double, -1, -1, Eigen::RowMajor> cumulative_probability;

    [[nodiscard]] CmaRead::L2DView<const double> get_view_cum_prob() const;

    [[nodiscard]] std::span<const double> get_diag_transition() const;

    std::vector<double> inverse_volume;
    std::vector<double> diag_transition;

    PreCalculatedHydroState() = default;
    ~PreCalculatedHydroState() = default;

    explicit PreCalculatedHydroState(const FlowMatrixType &_tm);

    PreCalculatedHydroState &operator=(PreCalculatedHydroState &&rhs) = default;
    PreCalculatedHydroState &
    operator=(const PreCalculatedHydroState &rhs) = delete;

    PreCalculatedHydroState(PreCalculatedHydroState &&other) = default;
    PreCalculatedHydroState(const PreCalculatedHydroState &other) = delete;
  };

  [[nodiscard]] inline CmaRead::L2DView<const double>
  PreCalculatedHydroState::get_view_cum_prob() const
  {
    return get_eigen_view(cumulative_probability);
  }

  [[nodiscard]] inline std::span<const double> PreCalculatedHydroState::get_diag_transition() const
  {
    return diag_transition;
  }

  struct TransitionState
  {
    CmaRead::ReactorState state;
    PreCalculatedHydroState liquid_pc;
    PreCalculatedHydroState gas_pc;
  };

} // namespace Simulation

#endif //__SIMUALTION__PC_HYDRO_HPP__