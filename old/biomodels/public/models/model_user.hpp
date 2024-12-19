#ifndef __BIO_MODEL_USER_HPP__
#define __BIO_MODEL_USER_HPP__
#include <mc/particles/particle_model.hpp>


namespace Models
{
  struct UserImpl;
  struct User
  {
  public:

    User(const User&)=default;

    User& operator=(const User&)=default;
    User(User&&) noexcept=default;
    User& operator=(User&&) noexcept=default;
    User();
    ~User()=default;
    KOKKOS_FUNCTION void init(MC::ParticleDataHolder& p, MC::KPRNG _rng);
    KOKKOS_FUNCTION void update(double d_t,
                                MC::ParticleDataHolder& p,
                                const LocalConcentrationView& concentration,
                                MC::KPRNG _rng);
    KOKKOS_FUNCTION User division(MC::ParticleDataHolder& p, MC::KPRNG k) noexcept;
    KOKKOS_FUNCTION void contribution(MC::ParticleDataHolder& p, const ContributionView& contrib) noexcept;

    model_properties_detail_t get_properties() noexcept;
    [[nodiscard]] KOKKOS_FUNCTION double mass() const noexcept;
    UserImpl* pimpl = nullptr;

  private:
    
    
  };
  static_assert(ParticleModel<User>, "Check Pimpl");
} // namespace Models
#endif
