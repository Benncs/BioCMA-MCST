#include "models/types.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <common/alg.hpp>
#include <mc/particles/mcparticles.hpp>
#define EEIGEN
#include <models/simple_model.hpp>

#include <any>


constexpr double SOLVER_TOLERANCE = 1e-5;

static double phi_pts(SimpleModel&model, double S);
static double
phi_permease(SimpleModel&model, double n_permease, double S);
static double uptake_glucose(SimpleModel&model,
                             unsigned n_permease,
                             double S);
static void update_xi_dot(SimpleModel&model,
                          double gamma_PTS_S,
                          double n_permease,
                          double S);

static double uptake_o2(SimpleModel&model, double O);

SimpleModel &SimpleModel::operator=(SimpleModel &&other) noexcept
{
  if (this != &other)
  {
    xi = other.xi;
    xi_dot = other.xi_dot;
  }
  return *this;
}

Xi Xi::operator+(const Xi &rhs) const
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
Xi Xi::operator*(double scalar) const
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
  this->xi = this->xi + this->xi_dot * d;
}

void init_simple_model_(MC::Particles &p)
{
  // std::shared_ptr<SimpleModel> ptr = std::make_shared<SimpleModel>();
  

  p.data = SimpleModel();
  auto& model = std::any_cast<SimpleModel&>(p.data);

  static double n_part = p.weight;
  static constexpr double total_volume  = 90./1000.; //TODO FIXME 
  static constexpr double avg_mass = SimpleModel::critcal_division_mass*0.1;
  static constexpr  double X_goal = 10;
  p.weight = (X_goal*(total_volume)/(n_part*avg_mass));



  // ptr->xi.mass =
  //     std::max(1e-7,
  //              MC::normal_double_rand(SimpleModel::critcal_division_mass *
  //              0.7,
  //                                     std::sqrt(1e-7)));

  model.xi.mass = SimpleModel::critcal_division_mass*0.1;
  // MC::uniform_double_rand(1e-7, SimpleModel::critcal_division_mass);
  model.xi.a_permease = 0.5;
  model.xi.a_pts = 0.5l;
  model.xi.mu_eff = 1e-5;
  model.xi.n_permease = 1e3;

  model.mu_p = 0;
  model.phi_o_in = 0;
  model.phi_s_in = 0;
  
  // std::cout<< ptr->xi->mass <<std::endl;
  assert(model.xi.mass > 0);
}
void update_simple_model(double d_t,
                         MC::Particles &p,
                         std::span<double const> concentrations)
{
  auto &model = std::any_cast<SimpleModel&> (p.data);
  double S = concentrations[0];
  auto n_permease = static_cast<unsigned>(model.xi.n_permease);

  double gamma_PTS_S = uptake_glucose(model, n_permease, S);

  double phi_o_in = uptake_o2(model, concentrations[1]);

  // Not enough substrate to stay alive
  // if (model.phi_s_in == 0 || phi_o_in == 0)
  // {
  //   p.status = MC::CellStatus::DEAD;
  //   return;
  // }

  double mu_p =
      std::min(SimpleModel::YXS * model.phi_s_in, SimpleModel::YXO * phi_o_in);

  model.mu_p = mu_p;
  model.phi_o_in = phi_o_in;

  update_xi_dot(model, gamma_PTS_S, n_permease, S);

  model.step(d_t);

  // Division has to occur
  if (model.xi.mass >= SimpleModel::critcal_division_mass)
  {
    p.status = MC::CellStatus::CYTOKINESIS;
  }
}

MC::Particles division_simple_model(MC::Particles &p)
{
  auto &model = std::any_cast<SimpleModel&> (p.data);
  model.xi.mass = model.xi.mass / 2;

  p.status = MC::CellStatus::IDLE;
  MC::Particles child(p);
  child.status = MC::CellStatus::IDLE;
  child.data.emplace<SimpleModel>(model) ;//= SimpleModel(model);//std::make_shared<SimpleModel>(*model);

  return child;
  // TODO
}
void contribution_simple_model(MC::Particles &p, Eigen::MatrixXd &contribution)
{
  auto &model = std::any_cast<SimpleModel&> (p.data);
  int ic = static_cast<int>(p.current_container);
  const double weight = p.weight;
  const auto z = static_cast<model_float_type_t>(0.);
  contribution.coeffRef(0, ic) += -model.phi_s_in * weight;
  contribution.coeffRef(1, ic) += -model.phi_o_in * weight;
  contribution.coeffRef(2, ic) +=
      SimpleModel::YSA * std::max(z, model.mu_p - model.xi.mu_eff) * weight;
}

static double phi_pts(SimpleModel&model, double S)
{
  return model.xi.a_pts *
         (SimpleModel::phi_pts_max * S / (SimpleModel::kpts + S));
}

static double
phi_permease(SimpleModel&model, double n_permease, double S)
{
  return n_permease * SimpleModel::psi_permease * (model.xi.mass * 1e3) *
         model.xi.a_permease *
         (SimpleModel::phi_pts_max * S / (SimpleModel::kppermease + S));
}

static double uptake_glucose(SimpleModel&model,
                             unsigned n_permease,
                             double S)
{

  if (S == 0)
  {
    model.phi_s_in = 0.;
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

  model.phi_s_in = phi_s_in;

  double gamma_PTS_S = phi_s_pts / SimpleModel::phi_pts_max;
  return gamma_PTS_S;
};

static void update_xi_dot(SimpleModel &model,
                          double gamma_PTS_S,
                          double n_permease,
                          double S)
{
  model.xi_dot.mass = model.xi.mu_eff;

  model.xi_dot.a_pts = (1.0 / SimpleModel::tauPTS) *
                        ((S / SimpleModel::kpts + S) - model.xi.a_pts);

  model.xi_dot.a_permease =
      ((1.0 / SimpleModel::tauAu) * gamma_PTS_S +
       (1.0 / SimpleModel::tauAd) * (1.0 - gamma_PTS_S)) *
      (1.0 - gamma_PTS_S - model.xi.a_permease);

  model.xi_dot.n_permease =
      (1.0 / SimpleModel::tau_f) *
          (SimpleModel::kppermease / (SimpleModel::kppermease + S)) +
      (1.0 / SimpleModel::tau_d) * (S / (SimpleModel::kpts + S)) *
          (SimpleModel::NPermease_max * (1.0 - gamma_PTS_S) - n_permease);
  model.xi_dot.mu_eff =
      (1.0 / SimpleModel::tau_metabolisme) * (model.mu_p - model.xi.mu_eff);
}

static double uptake_o2(SimpleModel&model, double O)
{
  static constexpr double tau_m = 1e-3;
  int growth = 0;
  double factor = 100; //TODO DELETE, JUST TO SEE OXYGENE COMSUMPTION

  auto get_phi = [&model, O, &growth,factor](double Oi)
  {
    double phi_o_growth = SimpleModel::YXO * model.xi.mu_eff;
    double phi_o_in =
        model.xi.mass * SimpleModel::psi_o_meta*factor + growth * phi_o_growth;
    double rhs = phi_o_growth + phi_o_in;
    double lhs = (O - Oi) / tau_m;
    return std::abs(rhs - lhs);
  };
  if (O == 0)
  {
    return 0;
  }
  bool success = false;
  double Oi = naive_newton(get_phi, O, &success, SOLVER_TOLERANCE);
  if (Oi < 0 || !success)
  {
    return 0;
  }

  growth = 1; //  Don't remove this, captured in get_phi TODO avoid side effect
              //  in get_phi
  double Oi2 = naive_newton(get_phi, O, &success, SOLVER_TOLERANCE);
  if (Oi2 < 0 || !success)
  {
    return model.xi.mass * SimpleModel::psi_o_meta*factor;
  }

  double phi_o_growth = SimpleModel::YXO * model.xi.mu_eff;
  return model.xi.mass * SimpleModel::psi_o_meta*factor + growth * phi_o_growth;
}

KModel get_simple_model()
{
  return {&init_simple_model_,
          &update_simple_model,
          &division_simple_model,
          &contribution_simple_model};
}