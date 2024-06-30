#ifndef __SIMULATION_HYDRO_MAsS_TRANSFER_HPP__
#define __SIMULATION_HYDRO_MAsS_TRANSFER_HPP__

#include <cma_read/reactorstate.hpp>
#include <mc/domain.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

Eigen::MatrixXd
gas_liquid_mass_transfer(Eigen::ArrayXXd &res_kla,
                         const Eigen::MatrixXd &Vliq,
                         const Eigen::ArrayXXd &liq_scalar_as_array,
                         const Eigen::ArrayXXd &gas_scalar_as_array,
                         const CmaRead::ReactorState &state);

#endif //__SIMULATION_HYDRO_MAsS_TRANSFER_HPP__ __SIM_HYDRO_MASS_TRANSFER__
