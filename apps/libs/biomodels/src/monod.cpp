#include "common/common.hpp"
#include <any>
#include <cmath>
#include <cstdlib>
#include <mc/prng/cdf.hpp>
#include <models/monod.hpp>
#include <random>
#include <Eigen/Core>
#include <iostream>

constexpr double interdivision_time = 20*60;
constexpr double ln2 = 0.6931471805599453;
constexpr double Ks = 486e-9;
constexpr double k_pts = 1e-3;
constexpr double phi_pts_max = 4.454e-12*1e-3/3600;
constexpr double YXS = 0.3;

std::random_device rd{};
std::mt19937 gen{rd()};
static std::uniform_real_distribution<double> __model_dist(0., 1.);

static std::normal_distribution<double> __model_n_dist(interdivision_time,interdivision_time/30.);

static double phi_pts(double S)
{
  return (phi_pts_max * S / (k_pts + S));
}

void init_monod_model_(MC::Particles &p)
{
  p.data = Monod{interdivision_time, 0, 0,6e-6};
}

void update_monod_model(double d_t,
                        MC::Particles &p,
                        std::span<double const> concentrations)
{
  auto &model = std::any_cast<Monod &>(p.data);
  model.age += d_t;
 
  double s = concentrations[0];
  model.phi = std::max(0.,phi_pts(s));

  model.l+=d_t*model.phi/YXS;  
  //  const auto lambda = model.age / model.interdivision_time;
  auto random = __model_dist(gen);

  // bool division = MC::cdf_truncated_normal_check(random,model.age,model.interdivision_time,model.interdivision_time/10.) &&! almost_equal(s, 0.);

  double tt = __model_n_dist(gen);
  while(tt<0)
  {
    tt = __model_n_dist(gen);
  }

  bool division = model.age>= tt && !almost_equal(s, 0.);
  if (division)
  {
    // std::cout<<model.age<<std::endl;
    p.status = MC::CellStatus::CYTOKINESIS;
  }
}

MC::Particles division_mond_model(MC::Particles &p)
{
  auto &model = std::any_cast<Monod &>(p.data);
  model.age = 0.;
  model.l = model.l/2.;
  p.status = MC::CellStatus::IDLE;
  // std::cout<<model.age<<std::endl;
  MC::Particles child(p);
  
  child.status = MC::CellStatus::IDLE;
  child.data.emplace<Monod>(
      model); //= SimpleModel(model);//std::make_shared<SimpleModel>(*model);

  return child;
}

void contribution_mond_model(MC::Particles &p, Eigen::MatrixXd &contribution)
{
  auto &model = std::any_cast<Monod &>(p.data);
  int ic = static_cast<int>(p.current_container);
  const double weight = p.weight;

  contribution.coeffRef(0, ic) -= model.phi * weight;
  // std::cout<<model.phi * weight<<std::endl;
}

KModel get_mond_model()
{

  return {init_monod_model_,
          update_monod_model,
          division_mond_model,
          contribution_mond_model};
}