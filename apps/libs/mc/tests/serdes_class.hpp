#ifndef __UTILS_SERDES_HPP__
#define __UTILS_SERDES_HPP__

#include <cassert>
#include <mc/traits.hpp>

// NOLINTBEGIN
struct SerdeModel
{
  using uniform_weight = std::true_type; // Using type alias
  using Self = SerdeModel;
  using FloatType = float;

  enum class particle_var : int
  {
    a = 0,
    b,
    c,
    COUNT
  };

  static constexpr std::size_t n_var = INDEX_FROM_ENUM(particle_var::COUNT);
  static constexpr std::string_view name = "simple";
  using SelfParticle = MC::ParticlesModel<Self::n_var, Self::FloatType>;

  MODEL_CONSTANT FloatType a_i = 5e-6; // m

  KOKKOS_INLINE_FUNCTION static void init([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
                                          std::size_t idx,
                                          const SelfParticle& arr)
  {
    GET_PROPERTY(SerdeModel::particle_var::a) = a_i;
  }

  KOKKOS_INLINE_FUNCTION static double mass([[maybe_unused]] std::size_t idx,
                                            [[maybe_unused]] const SelfParticle& arr)
  {
    return 1.;
  }

  KOKKOS_INLINE_FUNCTION static MC::Status
  update([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
         [[maybe_unused]] FloatType d_t,
         [[maybe_unused]] std::size_t idx,
         [[maybe_unused]] const SelfParticle& arr,
         [[maybe_unused]] const MC::LocalConcentration& c)
  {

    return MC::Status::Idle;
  }

  KOKKOS_INLINE_FUNCTION static void
  division([[maybe_unused]] const MC::KPRNG::pool_type& random_pool,
           [[maybe_unused]] std::size_t idx,
           [[maybe_unused]] std::size_t idx2,
           [[maybe_unused]] const SelfParticle& arr,
           [[maybe_unused]] const SelfParticle& buffer_arr)
  {
  }

  KOKKOS_INLINE_FUNCTION static void
  contribution([[maybe_unused]] std::size_t idx,
               [[maybe_unused]] std::size_t position,
               [[maybe_unused]] double weight,
               [[maybe_unused]] const SelfParticle& arr,
               [[maybe_unused]] const MC::ContributionView& contributions)
  {
  }
};
static_assert(ModelType<SerdeModel>, "Check non serde model");

// // NOLINTEND
// template <typename T> std::ostringstream wrap_ser(T& t)
// {
//   std::ostringstream buff(std::ios::binary);
//   cereal::BinaryOutputArchive ar(buff);
//   t.serialize(ar);
//   return buff;
// }

// template <typename T, class StreamType> void wrap_de(T& t, StreamType& buff)
// {
//   std::istringstream iss(buff.str(), std::ios::binary);
//   cereal::BinaryInputArchive ar(iss);
//   t.serialize(ar);
// }

#endif