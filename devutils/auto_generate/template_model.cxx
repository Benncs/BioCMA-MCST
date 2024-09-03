#include <models/@__model__name__@.hpp>

namespace Models
{
  KOKKOS_FUNCTION void @__model__name__@::init(MC::ParticleDataHolder &p, MC::KPRNG _rng)
  {

  }

  KOKKOS_FUNCTION void @__model__name__@::update(double d_t,
                              MC::ParticleDataHolder &p,
                              const LocalConcentrationView &concentration,
                              MC::KPRNG _rng)
                              {

                              }

  KOKKOS_FUNCTION @__model__name__@ @__model__name__@::division(MC::ParticleDataHolder &p)
  {
    return {};
  }

  KOKKOS_FUNCTION void @__model__name__@::contribution(MC::ParticleDataHolder &p,
                                    ContributionView contri)
                                    {

                                    }

  model_properties_detail_t @__model__name__@::get_properties()
  {
    return {};
  }

}; // namespace Models
