#ifndef __COMMON_EXPORT_HPP__
#define __COMMON_EXPORT_HPP__

#include <common/alg.hpp>
#include <common/common_types.hpp>
#include <common/execinfo.hpp>
#include <common/simulation_parameters.hpp>

void test_common();

#define EIGEN_INDEX(__VALUE__) static_cast<int>(__VALUE__)

// consteval int eigen_i(size_t v){return static_cast<int>(v);}

// #define TEST_ASSERT __asse

#endif //__COMMON_EXPORT_HPP__