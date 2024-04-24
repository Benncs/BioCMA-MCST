#ifndef __MODELS_TYPES_HPP__
#define __MODELS_TYPES_HPP__

#include <common/cmodel_parameters.hpp>
#include <mc/particles/mcparticles.hpp>
#include <span>

#include <Eigen/Core>

#ifdef USE_PYTHON_MODULE
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

// #ifdef DEBUG
using ModelDebug = std::function<void(MC::Particles &)>;
inline void defaut_dgb(MC::Particles & /*unused*/){};
// #endif

struct KModel
{
  ModelInit init_kernel;
  ModelUpdate update_kernel;
  ModelDivision division_kernel;
  ModelContribution contribution_kernel;
  // #ifdef DEBUG
  ModelDebug f_dbg = defaut_dgb;
  // #endif
};

// template <ModelParameter M>
// void update_kernel(double dt,
//                    std::shared_ptr<M> &model,
//                    std::span<const double>)
// {

//   model->step(dt);
// }

#endif //__MODELS_TYPES_HPP__