#include <models/model_user.hpp>

#define CHECK_PIMP                                                                                 \
  if (this->pimpl == nullptr)                                                                      \
  {                                                                                                \
    return;                                                                                        \
  }

namespace Models
{
  KOKKOS_FUNCTION void User::init(MC::ParticleDataHolder& p, MC::KPRNG _rng)
  {
    CHECK_PIMP;
  }

  KOKKOS_FUNCTION void User::update(double d_t,
                                         MC::ParticleDataHolder& p,
                                         const LocalConcentrationView& concentration,
                                         MC::KPRNG _rng)
  {
    CHECK_PIMP;
  }
  KOKKOS_FUNCTION User User::division(MC::ParticleDataHolder& p, MC::KPRNG k) noexcept
  {
    return {};
  }
  KOKKOS_FUNCTION void User::contribution(MC::ParticleDataHolder& p,
                                               ContributionView contrib) noexcept
  {
    CHECK_PIMP;
  }

  model_properties_detail_t User::get_properties() noexcept
  {
    return {};
  }

  KOKKOS_FUNCTION [[nodiscard]] double User::mass() const noexcept
  {
    return 0.;
  }

} // namespace Models