#ifndef __SIMPLE_MODEL_HPP__
#define __SIMPLE_MODEL_HPP__

#include <mc/particles/mcparticles.hpp>
#include <memory>

#include <models/types.hpp>

#include <span>

struct Xi
{

  double mass{};
  double mu_eff{};
  double a_pts{};
  double a_permease{};
  double n_permease{};

  ~Xi() = default;

  Xi operator+(const Xi &rhs) const;
  Xi operator*(double scalar) const;

  Xi() = default;
  Xi(const Xi &rhs) = default;
  Xi(Xi &&rhs) = default;
  Xi &operator=(const Xi &) = default;
  Xi &operator=(Xi &&) = default;
};

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
  static constexpr double psi_o_meta = 20e-3 / 3600 * 32e-3; // 20mmmol O2 /h;
  static constexpr double critcal_division_mass = 1.7e-5;

  double phi_s_in = 0.;
  double phi_o_in = 0.;
  double mu_p = 0;
  Xi xi;
  Xi xi_dot;
  void step(double d) ;
  SimpleModel() = default;

  ~SimpleModel() = default;
  SimpleModel(SimpleModel &&) = default;
  SimpleModel(const SimpleModel &rhs)
      : xi(rhs.xi), xi_dot(rhs.xi_dot){

                    };

  SimpleModel &operator=(const SimpleModel &other) noexcept = delete;

  SimpleModel &operator=(SimpleModel &&other) noexcept;
};

void init_simple_model_(MC::Particles &p);

void update_simple_model(double d_t,
                         MC::Particles &p,
                         std::span<double const> concentrations);

MC::Particles division_simple_model(MC::Particles &p);

void contribution_simple_model(MC::Particles &p, Eigen::MatrixXd &contribution);


KModel get_simple_model();



#endif //__SIMPLE_MODEL_HPP__
