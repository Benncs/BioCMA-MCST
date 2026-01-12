#ifndef __COMMON_ALG_HPP__
#define __COMMON_ALG_HPP__
#include <functional>
#include <limits>

double naive_newton(const std::function<double(double)>& f,
                    double x_init,
                    bool* success,
                    double tolerance = std::numeric_limits<double>::epsilon());
#endif //__COMMON_ALG_HPP__