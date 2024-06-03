#ifndef __SIMUALTION__PC_HYDRO_HPP__
#define __SIMUALTION__PC_HYDRO_HPP__

#include "cma_read/reactorstate.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

using FlowMatrixType = Eigen::MatrixXd;

namespace Simulation
{

  class PreCalculatedHydroState
  {
    public:
    // FlowMatrixType flows;
    FlowMatrixType transition_matrix;
    FlowMatrixType cumulative_probability;
    std::vector<double> inverse_volume;

    PreCalculatedHydroState() = default;
    ~PreCalculatedHydroState() = default;

    explicit PreCalculatedHydroState(const FlowMatrixType &_tm);

    PreCalculatedHydroState &operator=(PreCalculatedHydroState &&rhs) noexcept;
    PreCalculatedHydroState &operator=(const PreCalculatedHydroState &rhs) = delete;

    PreCalculatedHydroState(PreCalculatedHydroState &&other) noexcept;
    PreCalculatedHydroState(const PreCalculatedHydroState &other) = delete;
  };

  struct TransitionState
  {
    ReactorState state;
    PreCalculatedHydroState pc;
  };



} // namespace Simulation

#endif //__SIMUALTION__PC_HYDRO_HPP__