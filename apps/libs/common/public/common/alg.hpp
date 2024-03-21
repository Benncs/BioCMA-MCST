#ifndef __COMMON_ALG_HPP__
#define __COMMON_ALG_HPP__
#include <functional>

double naive_newton(const std::function<double(double)>& f, double x_init, bool *sucess,double tolerance) ;
#endif //__COMMON_ALG_HPP__