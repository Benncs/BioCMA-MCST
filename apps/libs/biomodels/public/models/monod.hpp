#ifndef __MONOD_MODEL_HPP__
#define __MONOD_MODEL_HPP__

#include <mc/particles/mcparticles.hpp>
#include <memory>

#include <models/types.hpp>

#include <span>

struct Monod
{
  double interdivision_time;
  double age;
  double phi;
  double l=6e-6;
};

void init_monod_model_(MC::Particles &p);

void update_monod_model(double d_t,
                         MC::Particles &p,
                         std::span<double const> concentrations);

MC::Particles division_mond_model(MC::Particles &p);

void contribution_mond_model(MC::Particles &p, Eigen::MatrixXd &contribution);

KModel get_mond_model();

#endif