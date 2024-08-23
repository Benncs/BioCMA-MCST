#ifndef __LIGHT_MODEL_HPP__
#define __LIGHT_MODEL_HPP__

#include <mc/particles/mcparticles.hpp>
#include <memory>

#include <models/types.hpp>

#include <span>

struct LightModel
{
  double mass;
};

// void init_light_model_(MC::Particle &p);

// void update_light_model(double d_t,
//                          MC::Particle &p,
//                          std::span<double const> concentrations);

// MC::Particle division_light_model(MC::Particle &p);

// void contribution_light_model(MC::Particle &p, Eigen::MatrixXd &contribution);

KModel get_model_light();

#endif