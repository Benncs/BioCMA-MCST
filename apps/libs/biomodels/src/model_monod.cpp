#include "common/common.hpp"
#include <Eigen/Core>
#include <any>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <mc/prng/cdf.hpp>
#include <models/monod.hpp>
#include <random>


constexpr double interdivision_time = 1175;//20*60 ;
constexpr double ln2 = 0.6931471805599453;
constexpr double Ks = 486e-9;
constexpr double k_pts = 1e-3;
constexpr double phi_pts_max =
    4.454e-12 * 1e-3 / 3600 ; //*0.1 to divide cumption by 10;
constexpr double YXS = 0.3;

std::random_device rd{};
std::mt19937 gen{rd()};
static std::uniform_real_distribution<double> __model_dist(0., 1.);

static std::normal_distribution<double> __model_n_dist(interdivision_time, 120);

static std::normal_distribution<double> n_dist(1, 0.1);

static double phi_pts(double S)
{
  return (phi_pts_max * S / (k_pts + S));
}

double rand_division_time()
{
  return n_dist(gen) * interdivision_time;
}

void init_monod_model_(MC::Particles &p)
{
  p.data = Monod{rand_division_time(), 0., 0., 6e-2};
}

void update_monod_model(double d_t,
                        MC::Particles &p,
                        std::span<double const> concentrations)
{
  auto &model = std::any_cast<Monod &>(p.data);
  model.age += d_t;

  double s = concentrations[0];
  model.phi = std::max(0., phi_pts(s));
  model.l = model.l + d_t*model.phi / YXS; 


  // const auto lambda = (model.age / model.interdivision_time);
  //  const auto lambda = (model.age * model.interdivision_time);
  // auto random = __model_dist(gen);
  // bool division = MC::cdf_exponential_check(random,lambda)&&! almost_equal(s,
  // 0.);

  bool division = model.age >= model.interdivision_time && !almost_equal(s, 0.);
  // auto random = __model_dist(gen);
  // bool division =
  // MC::cdf_truncated_normal_check(random,model.age,model.interdivision_time,model.interdivision_time/6.)
  // &&! almost_equal(s, 0.);

  if (division)
  {
    p.status = MC::CellStatus::CYTOKINESIS;
  }
}

MC::Particles division_mond_model(MC::Particles &p)
{
  auto &model = std::any_cast<Monod &>(p.data);
  model.age = 0.;
  model.l = model.l / 2.;
  p.status = MC::CellStatus::IDLE;
  MC::Particles child(p);

  child.status = MC::CellStatus::IDLE;
  child.data.emplace<Monod>(
      model); //= SimpleModel(model);//std::make_shared<SimpleModel>(*model);

  auto &child_model = std::any_cast<Monod &>(child.data);
  
  child_model.interdivision_time = model.interdivision_time; //rand_division_time();
  return child;
}

void contribution_mond_model(MC::Particles &p, Eigen::MatrixXd &contribution)
{
  auto &model = std::any_cast<Monod &>(p.data);
  int ic = static_cast<int>(p.current_container);
  const double weight = p.weight;

  contribution.coeffRef(0, ic) -= model.phi * weight;
}

model_properties_detail_t monod_properties(const MC::Particles &p)
{
  const auto &model = std::any_cast<const Monod &>(p.data);
  
  return {{"age", model.age},
          {"interdivision_time", model.interdivision_time},
          {"lenght", model.l}};
}

KModel get_model_monod()
{

  return {init_monod_model_,
          update_monod_model,
          division_mond_model,
          contribution_mond_model,
          monod_properties};
}