#ifndef __SIMPLE_MODEL_HPP__
#define __SIMPLE_MODEL_HPP__

#include <mc/particles/mcparticles.hpp>
#include <memory>

#include <models/types.hpp>

#include <span>

struct SimpleModel
{
  static constexpr double tauPTS = 25.0;
  static constexpr double tau_f = 5.0;
  static constexpr double tau_d = 5.0;
  static constexpr double NPermease_max = 5e4;
  static constexpr double tauAu = 5.0;
  static constexpr double tauAd = 5.0;
  static constexpr double tau_metabolisme = 0.1;
  static constexpr double phi_pts_max = 4.454e-12;
  static constexpr double kpts = 1e-3;
  static constexpr double kppermease = 1e-2;
  static constexpr double psi_permease = 1.25e-13;
  static constexpr double YXS = 0.5;
  static constexpr double YXO = 1e-4;
  static constexpr double YSA = 1e-4;
  static constexpr double psi_o_meta = 1e-12;
  static constexpr double critcal_division_mass = 1.1e-5;
  struct Xi;
  double phi_s_in =0.;
  double phi_o_in=0.;
  double mu_p;
  std::unique_ptr<Xi> xi;
  std::unique_ptr<Xi> xi_dot;
  void step(double d);
  SimpleModel() : xi(std::make_unique<Xi>()), xi_dot(std::make_unique<Xi>())
  {
  }
  SimpleModel &operator=(SimpleModel &&other) noexcept;
};

struct SimpleModel::Xi
{

  double mass;
  double mu_eff;
  double a_pts;
  double a_permease;
  double n_permease;

  Xi operator+(const Xi &rhs) const;
  Xi operator*(double scalar) const;

  Xi() : mass(0.), mu_eff(0), a_pts(0), a_permease(0), n_permease(0)
  {
  }
  Xi(const Xi &rhs)
      : mass(rhs.mass), mu_eff(rhs.mu_eff), a_pts(rhs.a_pts),
        a_permease(rhs.a_permease), n_permease(rhs.n_permease)
  {
  }
};

void init_simple_model_(MC::Particles &p);

void update_simple_model(double d_t,
                         MC::Particles &p,
                         std::span<double const> concentrations);

MC::Particles division_simple_model(MC::Particles &p);

void contribution_simple_model(MC::Particles &p, Eigen::MatrixXd &contribution);

static KModel simple_model = {&init_simple_model_,
                              &update_simple_model,
                              &division_simple_model,
                              &contribution_simple_model};

#endif //__SIMPLE_MODEL_HPP__
