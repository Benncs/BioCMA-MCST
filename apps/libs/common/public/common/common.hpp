#ifndef __COMMON_EXPORT_HPP__
#define __COMMON_EXPORT_HPP__

#include <common/common_types.hpp>
#include <common/execinfo.hpp>
#include <common/simulation_parameters.hpp>
#include <common/alg.hpp>

#include <Eigen/Sparse>
Eigen::SparseMatrix<double>  FlowmapToMat(std::span<double> data,size_t n_row);

void test_common();

#endif //__COMMON_EXPORT_HPP__