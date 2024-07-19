#ifndef __MONOD_AGE_MODEL_HPP__
#define __MONOD_AGE_MODEL_HPP__

#include <mc/particles/mcparticles.hpp>
#include <memory>

#include <models/types.hpp>

#include <span>

struct Monod
{
  double interdivision_time;
  double age;
  double phi;
  double l;
};

void init_monod_model_(MC::Particles &p);

void update_monod_model(double d_t,
                         MC::Particles &p,
                         std::span<double const> concentrations);

MC::Particles division_mond_model(MC::Particles &p);

void contribution_mond_model(MC::Particles &p, Eigen::MatrixXd &contribution);

model_properties_detail_t monod_properties(const MC::Particles &);

KModel get_model_monod_age();

#endif