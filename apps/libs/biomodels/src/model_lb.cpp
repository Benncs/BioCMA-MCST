#include "mc/prng/prng.hpp"
#include <Eigen/Dense>
#include <any>
#include <models/lb_model.hpp>

#include <random>
namespace lb
{

  constexpr bool implicit = false;

  static std::random_device rd{};
  static std::mt19937 gen{rd()};
  static std::normal_distribution<double> n_dist(1, 0.1);

  enum XI_N
  {
    lenght = 0,
    mu_eff = 1,
    a_pts = 2,
    a_permease = 3,
    n_permease = 4,
  };

  constexpr size_t xi_t_size = 5;
  using xi_t = Eigen::RowVector<double, xi_t_size>;

  struct LenghtBasedModel
  {

    MC::PRNG prng;
    xi_t xi;
    Eigen::RowVector<double, 3> contrib;

    static constexpr double length_to_mass_c = 1.;
    static constexpr double cell_diamater = 1e-6;

    static constexpr double tauPTS = 25.0;
    static constexpr double tau_f = 5.0;
    static constexpr double tau_d = 5.0;
    static constexpr double NPermease_max = 5e4;
    static constexpr double NPermease_init = 1e4;
    static constexpr double tauAu = 5.0;
    static constexpr double tauAd = 5.0;
    static constexpr double tau_metabolisme = 0.1;
    static constexpr double phi_pts_max = 4.454e-12;
    static constexpr double kpts = 1e-2;
    static constexpr double kppermease = 1e-3;
    static constexpr double psi_permease = 1.25e-13;
    static constexpr double YXS = 0.5;
    static constexpr double YXO = 1e-4;
    static constexpr double YSA = 1e-4;
    static constexpr double psi_o_meta = 20e-3 / 3600 * 32e-3; // 20mmmol O2 /h;

    static constexpr double critcal_division_length = 11e-6;

    static constexpr double minimal_length = 7e-6;
    static constexpr double maximal_length = 18e-6;
  };

  void init_lb_model_(MC::Particles &p)
  {

    // auto model = LenghtBasedModel();

    p.data = LenghtBasedModel();
    auto &model = std::any_cast<LenghtBasedModel &>(p.data);
    model.contrib.setZero();
    auto rd_length = n_dist(gen) * LenghtBasedModel::critcal_division_length;

    model.xi[XI_N::lenght] =
        std::max(LenghtBasedModel::minimal_length,
                 std::min(rd_length, LenghtBasedModel::maximal_length));

    model.xi[XI_N::a_permease] = model.prng.double_unfiform();
    model.xi[XI_N::a_pts] = model.prng.double_unfiform();
    model.xi[XI_N::n_permease] =
        model.prng.double_unfiform() * LenghtBasedModel::NPermease_init;
    model.xi[XI_N::mu_eff] = 1e-8 * model.prng.double_unfiform();

    // p.data = std::move(model);
  }

  inline double phi_pts(const auto &xi, double S)
  {
    return xi[XI_N::a_pts] *
           (LenghtBasedModel::phi_pts_max * S / (LenghtBasedModel::kpts + S));
  }

  inline double phi_permease(const auto &xi, double S)
  {
    return xi[XI_N::n_permease] * LenghtBasedModel::psi_permease *
           (xi[XI_N::lenght] * LenghtBasedModel::length_to_mass_c * 1e3) *
           xi[XI_N::a_permease] *
           (LenghtBasedModel::phi_pts_max * S /
            (LenghtBasedModel::kppermease + S));
  }

  inline double division_gamma(const auto &xi)
  {

    auto length = xi[XI_N::lenght];
    if (length <= LenghtBasedModel::minimal_length ||
        length > LenghtBasedModel::maximal_length)
    {
      return 0;
    }
    constexpr double kappa = 5;

    return std::pow(length / LenghtBasedModel::maximal_length, kappa);
  }

  inline void
  u_xi_dot(double gamma_PTS_S, double mu_p, double S, xi_t &xi, xi_t &xi_dot)
  {
    xi_dot[XI_N::lenght] = xi[XI_N::mu_eff];

    xi_dot[XI_N::a_pts] =
        (1.0 / LenghtBasedModel::tauPTS) *
        ((S / (LenghtBasedModel::kpts + S)) - xi[XI_N::a_pts]);

    xi_dot[XI_N::a_permease] =
        ((1.0 / LenghtBasedModel::tauAu) * gamma_PTS_S +
         (1.0 / LenghtBasedModel::tauAd) * (1.0 - gamma_PTS_S)) *
        (1.0 - gamma_PTS_S - xi[XI_N::a_permease]);

    xi_dot[XI_N::n_permease] =
        (1.0 / LenghtBasedModel::tau_f) * (LenghtBasedModel::kppermease /
                                           (LenghtBasedModel::kppermease + S)) +
        (1.0 / LenghtBasedModel::tau_d) * (S / (LenghtBasedModel::kpts + S)) *
            (LenghtBasedModel::NPermease_max * (1.0 - gamma_PTS_S) -
             xi[XI_N::a_permease]);
    xi_dot[XI_N::mu_eff] =
        (1.0 / LenghtBasedModel::tau_metabolisme) * (mu_p - xi[XI_N::mu_eff]);
  }

  void f(double s, xi_t &xi, xi_t &xi_dot)
  {

    const double phi_s_pts = phi_pts(xi, s);
    const double phi_s_in = phi_s_pts + phi_permease(xi, s);

    const double gamma_PTS_S = phi_s_pts / LenghtBasedModel::phi_pts_max;
    const double mu_p = LenghtBasedModel::YXS * phi_s_in;

    u_xi_dot(gamma_PTS_S, mu_p, s, xi, xi_dot);
  }

  void jacobian(double d_t,
                double s,
                Eigen::MatrixXd &J,
                xi_t &f_val,
                xi_t &xi,
                xi_t &xi_next,
                double epsilon)
  {

    for (int i = 0; i < xi.size(); ++i)
    {
      xi_t xi_next_eps = xi_next;
      xi_next_eps[i] += epsilon;
      xi_t xi_dot_eps(xi.size());
      f(s, xi_next_eps, xi_dot_eps);
      xi_t f_val_eps = xi_next_eps - xi - d_t * xi_dot_eps;
      J.col(i) = (f_val_eps - f_val) / epsilon;
    }
  }

  xi_t backward_euler_update(double d_t, double s, LenghtBasedModel &model)
  {
    auto &xi = model.xi;
    xi_t xi_dot = xi_t(xi);
    constexpr double tol = 1e-5;
    constexpr double max_iter = 100;
    xi_t xi_next = xi;
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(xi.size(), xi.size());
    for (int iter = 0; iter < max_iter; ++iter)
    {

      xi_t f_val = xi_next - xi - d_t * xi_dot;

      if (f_val.norm() < tol)
      {
        break;
      }

      jacobian(d_t, s, J, f_val, model.xi, xi_next, 1e-9);

      auto delta = J.fullPivLu().solve(-f_val.transpose());
      xi_next += delta.transpose();

      if (delta.norm() < tol)
      {
        break;
      }
    }
    return xi_next;
  }

  void update_lb_model(double d_t,
                       MC::Particles &p,
                       std::span<double const> concentrations)
  {
    auto &model = std::any_cast<LenghtBasedModel &>(p.data);
    auto &xi = model.xi;
    const double S = concentrations[0];

    if constexpr (implicit)
    {
      xi = backward_euler_update(d_t, S, model);
    }
    else
    {
      xi_t xi_dot;
      xi_dot.setZero();
      f(S, xi, xi_dot);
      xi += d_t * xi_dot;
    }

    model.contrib(0) = phi_pts(xi, S) + phi_permease(xi, S);

    auto proba_div = (1. - std::exp(-division_gamma(xi) * d_t));

    if (model.prng.double_unfiform() < proba_div)
    {
      p.status = MC::CellStatus::CYTOKINESIS;
    }
  }

  MC::Particles division_lb_model(MC::Particles &p)
  {
    auto &parent_model = std::any_cast<LenghtBasedModel &>(p.data);
    parent_model.xi[XI_N::lenght] = parent_model.xi[XI_N::lenght] / 2.;
    const auto current_perm = parent_model.xi[XI_N::n_permease];

    const double div_perm = parent_model.prng.double_unfiform();
    parent_model.xi[XI_N::n_permease] = div_perm * current_perm;

    p.status = MC::CellStatus::IDLE;
    MC::Particles child(p);

    auto child_model = LenghtBasedModel(parent_model);

    child_model.xi[XI_N::n_permease] = (1 - div_perm) * current_perm;
    child.status = MC::CellStatus::IDLE;

    child.data.emplace<LenghtBasedModel>(child_model);
    return child;
  }

  void contribution_lb_model(MC::Particles &p, Eigen::MatrixXd &contribution)
  {
    auto &model = std::any_cast<LenghtBasedModel &>(p.data);
    int ic = static_cast<int>(p.current_container);
    contribution.col(ic) -= (p.weight * model.contrib);
    // contribution.coeffRef(0,ic)-=p.weight * model.contrib[0];
  }

  model_properties_detail_t properties(const MC::Particles &p)
  {
    const auto &model = std::any_cast<const LenghtBasedModel &>(p.data);
    const auto &xi = model.xi;
    return {{"lenght", xi[XI_N::lenght]},
            {"a_pts", xi[XI_N::a_pts]},
            {"mu_eff", xi[XI_N::mu_eff]},
            {"a_permease", xi[XI_N::a_permease]},
            {"n_permease", xi[XI_N::n_permease]}};
  }
} // namespace lb

KModel get_model_lb()
{
  return {&lb::init_lb_model_,
          &lb::update_lb_model,
          &lb::division_lb_model,
          &lb::contribution_lb_model,
          &lb::properties};
}