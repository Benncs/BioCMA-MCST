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

void init_light_model_(MC::Particles &p);

void update_light_model(double d_t,
                         MC::Particles &p,
                         std::span<double const> concentrations);

MC::Particles division_light_model(MC::Particles &p);

void contribution_light_model(MC::Particles &p, Eigen::MatrixXd &contribution);

KModel get_light_model();

#endif