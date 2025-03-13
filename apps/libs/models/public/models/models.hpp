#ifndef __SIMPLE_MODEL_HPP__
#define __SIMPLE_MODEL_HPP__

#include <mc/prng/prng_extension.hpp>
#include <mc/traits.hpp>



namespace Models
{

  struct SimpleModel
  {

    enum class particle_var : int
    {
      mass = 0,
    };
    static constexpr std::size_t n_var = 1;
    static constexpr std::string_view name = "simple";
    using uniform_weight = std::true_type; // Using type alias
    using Self = SimpleModel;
    using FloatType = float;
    using SelfParticle = MC::ParticlesModel<Self::n_var, Self::FloatType>;

    MODEL_CONSTANT FloatType l_max_m = 5e-6;   // m
    MODEL_CONSTANT FloatType l_c_m = 3e-6;     // m
    MODEL_CONSTANT FloatType d_m = 0.6e-6;     // m
    MODEL_CONSTANT FloatType l_min_m = 0.9e-6; // m

    KOKKOS_INLINE_FUNCTION static void
    init([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
         std::size_t idx,
         const SelfParticle& arr);

    KOKKOS_INLINE_FUNCTION static MC::Status update(FloatType d_t,
                                                    std::size_t idx,
                                                    const SelfParticle& arr,
                                                    const MC::LocalConcentration& c);

    KOKKOS_INLINE_FUNCTION static void division(std::size_t idx,
                                                std::size_t idx2,
                                                const SelfParticle& arr,
                                                const SelfParticle& buffer_arr);

    KOKKOS_INLINE_FUNCTION static void contribution(std::size_t idx,
                                                    std::size_t position,
                                                    double weight,
                                                    const SelfParticle& arr,
                                                    const MC::ContributionView& contributions);

    KOKKOS_INLINE_FUNCTION static double mass(std::size_t idx, const SelfParticle& arr)
    {
      return arr(idx, static_cast<int>(particle_var::mass));
    }

    inline static std::vector<std::string> names()
    {
      return {"mass"};
    }
  };

  CHECK_MODEL(SimpleModel)

  KOKKOS_INLINE_FUNCTION void
  SimpleModel::init([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
                    std::size_t idx,
                    const SelfParticle& arr)
  {
    MODEL_CONSTANT auto l_distribution = MC::Distributions::TruncatedNormal<FloatType>(
        l_min_m, l_min_m / 5., l_min_m * 0.5, l_max_m);
    auto gen = random_pool.get_state();
    arr(idx, static_cast<int>(particle_var::mass)) = l_distribution.draw(gen);
    random_pool.free_state(gen);
  }

  KOKKOS_INLINE_FUNCTION MC::Status SimpleModel::update([[maybe_unused]] FloatType d_t,
                                                        [[maybe_unused]] std::size_t idx,
                                                        [[maybe_unused]] const SelfParticle& arr,
                                                        const MC::LocalConcentration& c)
  {
    (void)c;
    // arr(idx, static_cast<int>(particle_var::mass)) += d_t * 1e-15;
    // return (idx % 3 == 0) ? MC::Status::Division : MC::Status::Idle;
    return MC::Status::Idle;
  }

  KOKKOS_INLINE_FUNCTION void SimpleModel::division(std::size_t idx,
                                                    std::size_t idx2,
                                                    const SelfParticle& arr,
                                                    const SelfParticle& buffer_arr)
  {
    // for(size_t i=0;i<n_var;++i)
    // {
    //   buffer_arr(idx,i)=arr(idx,i); //basic copy
    // }
    const FloatType new_current_mass =
        arr(idx, static_cast<int>(particle_var::mass)) / static_cast<FloatType>(2.);
    buffer_arr(idx2, static_cast<int>(particle_var::mass)) = new_current_mass;
    arr(idx, static_cast<int>(particle_var::mass)) = new_current_mass;
  }

  KOKKOS_INLINE_FUNCTION void SimpleModel::contribution([[maybe_unused]] std::size_t idx,
                                                        std::size_t position,
                                                        double weight,
                                                        [[maybe_unused]] const SelfParticle& arr,
                                                        const MC::ContributionView& contributions)
  {
    auto access = contributions.access();
    access(position, 0) += -weight * 1e-15; //NOLINT
  }

} // namespace Models

#endif