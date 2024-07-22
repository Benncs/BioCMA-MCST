#ifndef __BIO_@__model__name__@_HPP__
#define __BIO_@__model__name__@_HPP__

#include <mc/particles/mcparticles.hpp>
#include <models/types.hpp>
#include <span>

struct @__model__name__@
{
  
};

void init_@__model__name__@(MC::Particles &p);

void update_@__model__name__@(double d_t,
                         MC::Particles &p,
                         std::span<double const> concentrations);

MC::Particles division_@__model__name__@(MC::Particles &p);

void contribution_@__model__name__@(MC::Particles &p, Eigen::MatrixXd &contribution);

KModel get_@__model__name__@();



#endif