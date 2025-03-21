// #ifdef DECLARE_EXPORT_UDF
#include "udf_includes.hpp"
#include <mc/particles/particle_model.hpp>
#include <models/udfmodel_user.hpp>

#define CHECK_PIMP                                                                                 \
  if (this->pimpl == nullptr)                                                                      \
  {                                                                                                \
    return;                                                                                        \
  }

namespace Models
{
  User::User()
  {
    UnsafeUDF::Loader::make_udf(*this);
  }

  KOKKOS_FUNCTION void User::init(MC::ParticleDataHolder& p, MC::KPRNG _rng)
  {
    CHECK_PIMP;
    return UnsafeUDF::Loader::init_udf(*pimpl, p,_rng);
  }

  KOKKOS_FUNCTION void User::update(double d_t,
                                    MC::ParticleDataHolder& p,
                                    const LocalConcentrationView& concentration,
                                    MC::KPRNG _rng)
  {
    CHECK_PIMP;
    return UnsafeUDF::Loader::update_udf(*pimpl, d_t, p, concentration);
  }
  KOKKOS_FUNCTION User User::division(MC::ParticleDataHolder& p, MC::KPRNG k) noexcept
  {
    if (this->pimpl == nullptr)
    {
      return {};
    }
    auto* const new_pimpl = UnsafeUDF::Loader::division_udf(*pimpl, p, k);
    if (new_pimpl != nullptr)
    {
      auto child = User(*this);
      child.pimpl = new_pimpl;
      return child;
    }

    return {};
  }
  KOKKOS_FUNCTION void User::contribution(MC::ParticleDataHolder& p,
                                          const ContributionView& contrib) noexcept
  {
    CHECK_PIMP;
    return UnsafeUDF::Loader::contribution_udf(*pimpl, p, contrib);
  }

  KOKKOS_FUNCTION [[nodiscard]] double User::mass() const noexcept
  {
    return UnsafeUDF::Loader::mass(*pimpl);
  }

  KOKKOS_FUNCTION void User::fill_properties(SubViewtype full) const
  {

    UnsafeUDF::Loader::fill_properties(*pimpl, full);
  }

  KOKKOS_FUNCTION std::size_t User::get_number()
  {
    return UnsafeUDF::Loader::get_number();
  }

  std::vector<std::string> User::names()
  {
    return UnsafeUDF::Loader::names();
  }

} // namespace Models

// #endif