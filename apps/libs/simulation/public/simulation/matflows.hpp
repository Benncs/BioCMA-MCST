#ifndef __SIMUALTION__MATFLOWS_HPP__
#define __SIMUALTION__MATFLOWS_HPP__

#include <Eigen/Dense>
#include <Eigen/Sparse>

using FlowMatrixType = Eigen::MatrixXd;

namespace Simulation
{

  class MatFlow
  {
    public:
    FlowMatrixType flows;
    FlowMatrixType transition_matrix;
    std::vector<double> inverse_volume;

    MatFlow() = default;
    ~MatFlow() = default;

    MatFlow(const FlowMatrixType &_flows, const FlowMatrixType &_tm);

    MatFlow &operator=(MatFlow &&rhs) noexcept;
    MatFlow &operator=(const MatFlow &rhs) = delete;

    MatFlow(MatFlow &&other) noexcept;
    MatFlow(const MatFlow &other) = delete;
  };

} // namespace Simulation

#endif //__SIMUALTION__MATFLOWS_HPP__