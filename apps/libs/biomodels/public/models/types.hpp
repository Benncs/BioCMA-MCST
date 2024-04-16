#ifndef __MODELS_TYPES_HPP__
#define __MODELS_TYPES_HPP__

#include <Eigen/Dense>
#include <common/cmodel_parameters.hpp>
#include <mc/particles/mcparticles.hpp>
#include <memory>
#include <span>

#ifdef BIO_DYNAMIC_MODULE
#  include <functional>

using ModelInit = std::function<void(MC::Particles &)>;

using ModelUpdate =
    std::function<void(double, MC::Particles &m, std::span<const double>)>;

using ModelDivision = std::function<MC::Particles(MC::Particles &)>;

using ModelContribution =
    std::function<void(MC::Particles &, Eigen::MatrixXd &)>;
#else
using ModelUpdate = void (*)(double, MC::Particles &, std::span<const double>);

using ModelDivision = MC::Particles (*)(MC::Particles &);

using ModelInit = void (*)(MC::Particles &);

using ModelContribution = void (*)(MC::Particles &, Eigen::MatrixXd &);

#endif

struct KModel
{
  ModelInit init_kernel;
  ModelUpdate update_kernel;
  ModelDivision division_kernel;
  ModelContribution contribution_kernel;
};

// template <ModelParameter M>
// void update_kernel(double dt,
//                    std::shared_ptr<M> &model,
//                    std::span<const double>)
// {

//   model->step(dt);
// }

#endif //__MODELS_TYPES_HPP__