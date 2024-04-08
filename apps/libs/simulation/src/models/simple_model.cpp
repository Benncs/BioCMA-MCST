#include "common/alg.hpp"
#include "mc/particles/mcparticles.hpp"
#include <cmath>
#include <simulation/models/simple_model.hpp>
#include <simulation/models/types.hpp>

#include <any>
#include <memory>

static double phi_pts(std::shared_ptr<SimpleModel> &model, double S);
static double
phi_permease(std::shared_ptr<SimpleModel> &model, double n_permease, double S);
static double uptake_glucose(std::shared_ptr<SimpleModel> &model,
                             unsigned n_permease,
                             double S);
static void update_xi_dot(std::shared_ptr<SimpleModel> &model,
                          double gamma_PTS_S,
                          double n_permease,
                          double S);

static double uptake_o2(std::shared_ptr<SimpleModel> &model, double O);

SimpleModel &SimpleModel::operator=(SimpleModel &&other)
{
  if (this != &other)
  {
    *xi = *other.xi;
    *xi_dot = *other.xi_dot;
  }
  return *this;
}

SimpleModel::Xi SimpleModel::Xi::operator+(const Xi &rhs) const
{
  Xi result;
  result.mass = mass + rhs.mass;
  result.mu_eff = mu_eff + rhs.mu_eff;
  result.a_pts = a_pts + rhs.a_pts;
  result.a_permease = a_permease + rhs.a_permease;
  result.n_permease = n_permease + rhs.n_permease;
  return result;
}

// Scalar multiplication operator for Xi
SimpleModel::Xi SimpleModel::Xi::operator*(double scalar) const
{
  Xi result;
  result.mass = mass * scalar;
  result.mu_eff = mu_eff * scalar;
  result.a_pts = a_pts * scalar;
  result.a_permease = a_permease * scalar;
  result.n_permease = static_cast<unsigned>(n_permease * scalar);
  return result;
}

void SimpleModel::step(double d)
{
  *this->xi = *this->xi + *this->xi_dot * d;
}

void init_simple_model_(MC::Particles &p)
{
  std::shared_ptr<SimpleModel> ptr = std::make_shared<SimpleModel>();
  ptr->xi->mass = 1e-5;
  ptr->xi->a_permease = 0.5;
  ptr->xi->a_pts = 0.5l;
  ptr->xi->mu_eff = 1e-5;
  ptr->xi->n_permease = 1e3;

  ptr->mu_p = 0;
  ptr->phi_o_in = 0;
  ptr->phi_s_in = 0;
  p.data = ptr;
}

void update_simple_model(double d_t,
                         MC::Particles &p,
                         std::span<double const> concentrations)
{
  auto &model = std::any_cast<std::shared_ptr<SimpleModel> &>(p.data);

  double S = concentrations[0];
  unsigned n_permease = static_cast<unsigned>(model->xi->n_permease);

  double gamma_PTS_S = uptake_glucose(model, n_permease, S);
  ;

  double phi_o_in = uptake_o2(model, concentrations[1]);

  // Not enough substrate to stay alive
  if (model->phi_s_in == 0 || phi_o_in == 0)
  {
    p.status = MC::CellStatus::DEAD;
    return;
  }

  double mu_p =
      std::min(SimpleModel::YXS * model->phi_s_in, SimpleModel::YXO * phi_o_in);

  model->mu_p = mu_p;

  model->phi_o_in = phi_o_in;

  update_xi_dot(model, gamma_PTS_S, n_permease, S);

  model->step(d_t);

  // Division has to occur
  if (model->xi->mass >= SimpleModel::critcal_division_mass)
  {
    p.status = MC::CellStatus::CYTOKINESIS;
  }
}

MC::Particles division_simple_model(MC::Particles &p)
{
  auto &model = std::any_cast<std::shared_ptr<SimpleModel> &>(p.data);
  model->xi->mass = model->xi->mass/2;
  p.status = MC::CellStatus::IDLE;
  auto child = MC::Particles(p);
  return child;
  // TODO
}

void contribution_simple_model(MC::Particles &p, Eigen::MatrixXd &contribution)
{
  auto &model = std::any_cast<std::shared_ptr<SimpleModel> &>(p.data);
  int ic = static_cast<int>(p.current_container);
  contribution.coeffRef(0, ic) += -model->phi_s_in * p.weight;
  contribution.coeffRef(1, ic) += -model->phi_o_in * p.weight;
  contribution.coeffRef(2, ic) +=
      SimpleModel::YSA * std::max(0., model->mu_p - model->xi->mu_eff) *
      p.weight;
}

static double phi_pts(std::shared_ptr<SimpleModel> &model, double S)
{
  return model->xi->a_pts *
         (SimpleModel::phi_pts_max * S / (SimpleModel::kpts + S));
}

static double
phi_permease(std::shared_ptr<SimpleModel> &model, double n_permease, double S)
{
  return n_permease * SimpleModel::psi_permease * (model->xi->mass*1e3) * model->xi->a_permease *
         (SimpleModel::phi_pts_max * S / (SimpleModel::kppermease + S));
}

static double uptake_glucose(std::shared_ptr<SimpleModel> &model,
                             unsigned n_permease,
                             double S)
{

  if (S == 0)
  {
    model->phi_s_in = 0.;
    return 0.;
  }

  static constexpr double tau_m = 1e-3;
  auto get_phi = [&model, n_permease, S](double Si)
  {
    double phi_s_pts = phi_pts(model, Si);
    double phi_s_permease = phi_permease(model, n_permease, Si);
    double rhs = phi_s_pts + phi_s_permease;
    double lhs = (S - Si) / tau_m;
    return std::abs(rhs - lhs);
  };

  bool success = false;
  double Si = naive_newton(get_phi, S, &success, 1e-5);
  if (Si < 0 || !success)
  {
    return 0;
  }

  double phi_s_pts = phi_pts(model, Si);

  double phi_s_in = phi_s_pts + phi_permease(model, n_permease, Si);

  if (std::isnan(phi_s_in))
  {
    phi_s_in = 0;
  }

  model->phi_s_in = phi_s_in;

  double gamma_PTS_S = phi_s_pts / SimpleModel::phi_pts_max;
  return gamma_PTS_S;
};

#include <iostream>
static void update_xi_dot(std::shared_ptr<SimpleModel> &model,
                          double gamma_PTS_S,
                          double n_permease,
                          double S)
{
  model->xi_dot->mass = 1.5e-7;//model->xi->mu_eff;
  model->xi_dot->a_pts = (1.0 / SimpleModel::tauPTS) *
                         ((S / SimpleModel::kpts + S) - model->xi->a_pts);

  model->xi_dot->a_permease =
      ((1.0 / SimpleModel::tauAu) * gamma_PTS_S +
       (1.0 / SimpleModel::tauAd) * (1.0 - gamma_PTS_S)) *
      (1.0 - gamma_PTS_S - model->xi->a_permease);

  model->xi_dot->n_permease =
      (1.0 / SimpleModel::tau_f) *
          (SimpleModel::kppermease / (SimpleModel::kppermease + S)) +
      (1.0 / SimpleModel::tau_d) * (S / (SimpleModel::kpts + S)) *
          (SimpleModel::NPermease_max * (1.0 - gamma_PTS_S) - n_permease);
  model->xi_dot->mu_eff =
      (1.0 / SimpleModel::tau_metabolisme) * (model->mu_p - model->xi->mu_eff);
}

static double uptake_o2(std::shared_ptr<SimpleModel> &model, double O)
{
  static constexpr double tau_m = 1e-3;
  int growth = 0;

  auto get_phi = [&model, O, &growth](double Oi)
  {
    double phi_o_growth = SimpleModel::YXO * model->xi->mu_eff;
    double phi_o_in = SimpleModel::psi_o_meta + growth * phi_o_growth;
    double rhs = phi_o_growth + phi_o_in;
    double lhs = (O - Oi) / tau_m;
    return std::abs(rhs - lhs);
  };
  if (O == 0)
  {
    return 0;
  }
  bool success = false;
  double Oi = naive_newton(get_phi, O, &success, 1e-5);
  if (Oi < 0 || !success)
  {
    return 0;
  }

  growth = 1; //  Don't remove this, captured in get_phi TODO avoid side effect in get_phi  
  double Oi2 = naive_newton(get_phi, O, &success, 1e-5);
  if (Oi2 < 0 || !success)
  {
    return SimpleModel::psi_o_meta;
  }
  else
  {
    double phi_o_growth = SimpleModel::YXO * model->xi->mu_eff;
    return SimpleModel::psi_o_meta + growth * phi_o_growth;
  }
}