#ifndef __MC_CDF_HPP__
#define __MC_CDF_HPP__

#include <cmath>
#include <iostream> 
namespace MC
{

  inline double pdf_truncated_normal(double x, double mean, double stddev)
  {
    double sqrt2 = sqrt(2.0);
    return (1.0 / (stddev * sqrt2 * sqrt(M_PI))) *
           exp(-0.5 * pow((x - mean) / stddev, 2));
  }

  inline double cdf_truncated_normal(double x, double mean, double stddev)
  {
    return 0.5 * (1 + erf((x - mean) / (stddev * sqrt(2.0))));
  }

  inline bool cdf_truncated_normal_check(double random,
                                         double x,
                                         double mean,
                                         double stddev)
  {
    
    return random < (cdf_truncated_normal(x, mean, stddev));
  }

  inline double cdf_exponential(double x)
  {
    return 1-std::exp(-x);
  }

  inline bool cdf_exponential_check(double random, double x)
  {
    return random < cdf_exponential(x);
  }

} // namespace MC

#endif //__MC_CDF_HPP__