#ifndef __SIMPLE_MODEL_HPP__
#define __SIMPLE_MODEL_HPP__

struct SimpleModel
{
  double mu;
  double phi;
  double phi_a;
  static constexpr double YSX = 0.5;
  static constexpr double qSmax = 2;
  static constexpr double ks = 0.01;
  static constexpr double taum = 50;
  static constexpr double tauphi = 0.1;
  static constexpr double tauphi_ex = 0.1;
  static constexpr double tau_mu = 1;
};

#endif //__SIMPLE_MODEL_HPP__
