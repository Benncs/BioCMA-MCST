#ifndef __MODELS_TYPES_HPP__
#define __MODELS_TYPES_HPP__

#include "mc/particles/mcparticles.hpp"
#include <span>
#include "common/cmodel_parameters.hpp"
#include <memory>
#include <Eigen/Dense>

using ModelUpdate = void (*)(double ,MC::Particles &,std::span<const double> );

using ModelP= void (*)(MC::Particles &);

using ModelInit = ModelP;

using ModelContribution = void(*)(MC::Particles&, Eigen::MatrixXd& contribution );

struct KModel{
  ModelInit init_kernel;
  ModelUpdate update_kernel;
  ModelP division_kernel;
  ModelContribution contribution_kernel; 
};

template<ModelParameter M>
void update_kernel(double dt,std::shared_ptr<M>& model,std::span<const double>)
{   

  


  model->step(dt);
}


#endif //__MODELS_TYPES_HPP__