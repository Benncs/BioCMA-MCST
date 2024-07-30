#include "mc/prng/prng.hpp"
#include <models/light_model.hpp>

#include <iostream>
#include <random>
// values near the mean are the most likely
// standard deviation affects the dispersion of generated values from the mean
static double avg_mass = 1e-13;
static double max_mass = 1e-11;
void init_light_model_(MC::Particles &p)
{
  std::random_device rd{};
  std::mt19937 gen{rd()};
  static double n_part = p.weight;
  static constexpr double total_volume  = 90./1000.; //TODO FIXME 
  static constexpr  double X_goal = 0.5;
  p.weight = (X_goal*total_volume/(n_part*avg_mass));

  static auto d = std::normal_distribution<>(avg_mass, avg_mass / 10.);

  double mass = std::max(0., std::min(d(gen), max_mass));
  p.data = LightModel{mass};
}

void update_light_model(double d_t,
                        MC::Particles &p,
                        std::span<double const> concentrations)
{
  // auto &model = std::any_cast<LightModel &>(p.data);
  // model.mass += avg_mass;
}

MC::Particles division_light_model(MC::Particles &p)
{
  // NOPE
}

void contribution_light_model(MC::Particles &p, Eigen::MatrixXd &contribution)
{
    int ic = static_cast<int>(p.current_container);
    
    
    // contribution(0,0) -= (1e-5);
}

KModel get_model_light()
{
  return {init_light_model_,
          update_light_model,
          division_light_model,
          contribution_light_model};
}