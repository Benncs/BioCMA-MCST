#ifndef __SIMUALTION__MATFLOWS_HPP__
#define __SIMUALTION__MATFLOWS_HPP__

#include <Eigen/Dense>
#include <Eigen/Sparse>

using FlowMatrixType = Eigen::MatrixXd;

namespace Simulation
{


  struct MatFlow
  {
    FlowMatrixType flows;
    FlowMatrixType transition_matrix;

    MatFlow() = default;

    MatFlow(const FlowMatrixType &_flows, const FlowMatrixType &_tm);

    MatFlow(const FlowMatrixType &&_flows, const FlowMatrixType &&_tm);

    MatFlow &operator=(MatFlow &&rhs)noexcept;

    MatFlow(MatFlow &&other)noexcept;
  };

  struct VecMatFlows
  {
    std::vector<MatFlow> data;

    explicit VecMatFlows(size_t n)noexcept
    {
      data.resize(n);
    }
  };
} // namespace Simulation

#endif //__SIMUALTION__MATFLOWS_HPP__