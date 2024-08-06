#include "mc/prng/prng.hpp"
#include <Eigen/Dense>
#include <any>
#include <backward_euler.hpp>
#include <models/monod.hpp>

#include <random>
#include <stdexcept>
namespace monod
{

  constexpr bool implicit = true;
  static std::random_device rd{};
  static std::mt19937 gen{rd()};
  static std::normal_distribution<double> n_dist(1, 0.1);

  enum XI_N
  {
    lenght = 0,
    mu_eff = 1,
  };

  constexpr size_t xi_t_size = 2;
  using xi_t = Eigen::Vector<double, xi_t_size>;

  struct alignas(ExecInfo::cache_line_size) MonodSimple
  {

    xi_t xi;
    Eigen::Vector<double, 1> contrib;

    static constexpr double mu_max = 0.46 / 3600;

    static constexpr double YXS = 0.5;
    static constexpr double critcal_division_length = 11e-6;
    static constexpr double minimal_length = 7e-6;
    static constexpr double maximal_length = 18e-6;
    static constexpr double tau_metabolism = 1.25 / mu_max;
    static constexpr double KS = 0.01;
    static MC::KPRNG rng;
  };

  void init_monod_model_(MC::Particles &p)
  {

    p.data = MonodSimple();
    auto &model = std::any_cast<MonodSimple &>(p.data);
    model.contrib.setZero();

    model.xi[XI_N::lenght] = p.rng->double_unfiform(
        MonodSimple::minimal_length, MonodSimple::maximal_length);

    // model.xi[XI_N::mu_eff] =
    //     MonodSimple::mu_max * (1. - p.rng->double_unfiform());

    model.xi[XI_N::mu_eff]  = p.rng->double_unfiform(
        1e-8, 0.90*MonodSimple::mu_max);
  }

  inline double division_gamma(const auto &xi)
  {

    const auto length = xi[XI_N::lenght];
    if (length <= MonodSimple::critcal_division_length)
    {
      return 0.;
    }
    constexpr double kappa = 5.;

    return std::pow(length / MonodSimple::maximal_length, kappa);
  }

  inline void u_xi_dot(double mu_p, double S, const xi_t &xi, xi_t &xi_dot)
  {
    constexpr double cell_lenghtening = 6.51e-6;
    const auto tau_metabolism = 3600;//1.25 / MonodSimple::mu_max;
    // std::cout<<(mu_p - xi[XI_N::mu_eff])<<std::endl;
    xi_dot = {
        MonodSimple::YXS *  xi[XI_N::mu_eff] * cell_lenghtening,
        (1.0 / tau_metabolism) * (mu_p - xi[XI_N::mu_eff]),
    };
  }

  void f(double s, const xi_t &xi, xi_t &xi_dot)
  {
    const double mu_p = xi[XI_N::mu_eff] * s / (MonodSimple::KS + s);
    

    u_xi_dot(mu_p, s, xi, xi_dot);
  }

  void update_monod_model(double d_t,
                          MC::Particles &p,
                          std::span<double const> concentrations)
  {
    auto &model = std::any_cast<MonodSimple &>(p.data);
    auto &xi = model.xi;
    double S = std::max(concentrations[0], 0.);

    if (S < 1e-20)
    {
      S = 0.;
    }

    if constexpr (implicit)
    {
      xi = backward_euler_update<xi_t_size>(d_t, S, model.xi, &f);
    }
    else
    {
      xi_t xi_dot;
      xi_dot.setZero();

      f(S, xi, xi_dot);
      xi += d_t * xi_dot;
    }

    model.contrib(0) = xi[XI_N::mu_eff] * S / (MonodSimple::KS + S);

    auto proba_div = (1. - std::exp(-division_gamma(xi) * d_t));

    if (p.rng->double_unfiform() < proba_div)
    {
      p.status = MC::CellStatus::CYTOKINESIS;
    }
  }

  MC::Particles division_monod_model(MC::Particles &p)
  {

    auto &parent_model = std::any_cast<MonodSimple &>(p.data);
    const double div_length = p.rng->double_unfiform();

    const auto current_length = parent_model.xi[XI_N::lenght];

    MC::Particles child(p);

    auto child_model = MonodSimple(parent_model);
    child_model.xi = parent_model.xi;

    parent_model.xi[XI_N::lenght] = div_length * current_length;
    child_model.xi[XI_N::lenght] = (1 - div_length) * current_length;

    child.status = MC::CellStatus::IDLE;
    p.status = MC::CellStatus::IDLE;

    child.data.emplace<MonodSimple>(child_model);

    return child;

    // // parent_model.xi[XI_N::lenght]/=2.;
    // // p.weight*=2;
  }

  void contribution_monod_model(MC::Particles &p, Eigen::MatrixXd &contribution)
  {
    auto &model = std::any_cast<MonodSimple &>(p.data);

    // int ic = static_cast<int>(p.current_container);

    // contribution.col(0) -= (p.weight * model.contrib);

    contribution(0, 0) -= (p.weight * model.contrib(0));
  }

  model_properties_detail_t properties(const MC::Particles &p)
  {
    const auto &model = std::any_cast<const MonodSimple &>(p.data);
    const auto &xi = model.xi;
    return {
        {"lenght", xi[XI_N::lenght]},
        {"mu_eff", xi[XI_N::mu_eff]},
    };
  }
} // namespace monod

KModel get_model_monod()
{
  return {&monod::init_monod_model_,
          &monod::update_monod_model,
          &monod::division_monod_model,
          &monod::contribution_monod_model,
          &monod::properties};
}