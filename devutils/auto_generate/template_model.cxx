#include "mc/prng/prng.hpp"
#include <models/@__model__name__@.hpp>


void init_@__model__name__@(MC::Particles &p)
{
  p.data = __model__name__{};
}

void update_@__model__name__@(double d_t,
                        MC::Particles &p,
                        std::span<double const> concentrations)
{
  auto &model = std::any_cast<@__model__name__@ &>(p.data);
  
}

MC::Particles division_@__model__name__@(MC::Particles &p)
{
  auto &model = std::any_cast<@__model__name__@ &>(p.data);
}

void contribution_@__model__name__@(MC::Particles &p, Eigen::MatrixXd &contribution)
{
  auto &model = std::any_cast<@__model__name__@ &>(p.data);
}

KModel get_@__model__name__@()
{
  return {init_@__model__name__@,
          update_@__model__name__@,
          division_@__model__name__@,
          contribution_@__model__name__@};
}