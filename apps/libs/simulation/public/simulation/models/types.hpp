#ifndef __MODELS_TYPES_HPP__
#define __MODELS_TYPES_HPP__

#include "common/cmodel_parameters.hpp"
#include "mc/particles/mcparticles.hpp"
#include <Eigen/Dense>
#include <memory>
#include <span>

#ifdef BIO_DYNAMIC_MODULE
#include <functional>

using ModelInit = std::function<void(MC::Particles &)>;

using ModelUpdate =
    std::function<void(double, MC::Particles &m, std::span<const double>)>;

using ModelP = std::function<void(MC::Particles &)>;

using ModelContribution =
    std::function<void(MC::Particles &, Eigen::MatrixXd &)>;
#else
using ModelUpdate = void (*)(double, MC::Particles &, std::span<const double>);

using ModelP = void (*)(MC::Particles &);

using ModelInit = ModelP;

using ModelContribution = void (*)(MC::Particles &, Eigen::MatrixXd &);

#endif



struct KModel
{
  ModelInit init_kernel;
  ModelUpdate update_kernel;
  ModelP division_kernel;
  ModelContribution contribution_kernel;
};



template <ModelParameter M>
void update_kernel(double dt,
                   std::shared_ptr<M> &model,
                   std::span<const double>)
{

  model->step(dt);
}

#endif //__MODELS_TYPES_HPP__