#ifndef __BIO_LB_MODEL_HPP__
#define __BIO_LB_MODEL_HPP__

#include <mc/particles/mcparticles.hpp>

#include <models/types.hpp>

#include <span>



// void init_lb_model_(MC::Particles &p);

// void update_lb_model(double d_t,
//                      MC::Particles &p,
//                      std::span<double const> concentrations);

// MC::Particles division_lb_model(MC::Particles &p);

// void contribution_lb_model(MC::Particles &p, Eigen::MatrixXd &contribution);

KModel get_model_lb();

#endif