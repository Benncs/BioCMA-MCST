#ifndef __MODEL_PYTHON_HPP__
#define __MODEL_PYTHON_HPP__

#define EXPORT __attribute__((visibility("default")))
#include <mc/particles/particle_model.hpp>

namespace PythonWrap
{

  class EXPORT PimpModel
  {
  public:
    PimpModel(const PimpModel&);
 
    PimpModel& operator=(const PimpModel&);
    PimpModel(PimpModel&&)noexcept;
    PimpModel& operator=(PimpModel&&) noexcept ;
    PimpModel();
    ~PimpModel();
    KOKKOS_FUNCTION void init(MC::ParticleDataHolder& p, MC::KPRNG _rng);
    KOKKOS_FUNCTION void update(double d_t,
                                MC::ParticleDataHolder& p,
                                const LocalConcentrationView& concentration,
                                MC::KPRNG _rng);
    KOKKOS_FUNCTION PimpModel division(MC::ParticleDataHolder& p, MC::KPRNG k) noexcept;
    KOKKOS_FUNCTION void contribution(MC::ParticleDataHolder& p, const ContributionView& contrib) noexcept;

    model_properties_detail_t get_properties() noexcept;
    [[nodiscard]] KOKKOS_FUNCTION double mass() const noexcept;

  private:
    struct Impl;
    Impl* pimpl = nullptr;
  };
  static_assert(ParticleModel<PimpModel>, "Check Pimpl");
} // namespace PythonWrap

#endif 