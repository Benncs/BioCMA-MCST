#ifndef __BACKWARD_EULER_HPP__
#define __BACKWARD_EULER_HPP__

#include <Eigen/Core>
#include <Eigen/Dense>
template <size_t dim> using vec_type = Eigen::Vector<double, dim>;

namespace
{
  template <size_t dim, typename F>
  void jacobian(double d_t,
                double s,
                Eigen::Matrix<double, dim, dim> &J,
                const vec_type<dim> &f_val,
                const vec_type<dim> &xi,
                vec_type<dim> &xi_next,
                const F &f,
                double epsilon)
  {
    vec_type<dim> xi_next_eps;
    vec_type<dim> xi_dot_eps;
    vec_type<dim> f_val_eps;
    for (size_t i = 0; i < xi.size(); ++i)
    {
      xi_next_eps = xi_next;
      xi_next_eps[i] += epsilon;
      f(s, xi_next_eps, xi_dot_eps);
      f_val_eps = xi_next_eps - (xi + d_t * xi_dot_eps);
      J.col(i) = (f_val_eps - f_val) / epsilon;
    }
  }
} // namespace

template <size_t dim, typename F>
vec_type<dim>
backward_euler_update(double d_t, double s, vec_type<dim> &xi, const F &f)
{
  vec_type<dim> xi_next = xi;
  vec_type<dim> xi_dot = xi;

  vec_type<dim> f_val;
  vec_type<dim> delta;
  Eigen::Matrix<double, dim, dim> J = Eigen::Matrix<double, dim, dim>::Zero();

  constexpr double tol = 1e-12;
  constexpr double max_iter = 500;
  constexpr double jacobian_epsilon = 1e-10;

  for (int iter = 0; iter < max_iter; ++iter)
  {
    f_val = xi_next - (xi + d_t * xi_dot);

    ::jacobian<dim>(d_t, s, J, f_val, xi, xi_next, f, jacobian_epsilon);
    delta = J.fullPivLu().solve(-f_val);
    xi_next += delta;

    if (delta.norm() < tol)
    {
      break;
    }

    f(s, xi_next, xi_dot);
  }

  return xi_next;
}

#endif