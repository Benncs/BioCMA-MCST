#ifndef __MC_PARTICLE_MODEL_HPP__
#define __MC_PARTICLE_MODEL_HPP__

#include "common/kokkos_vector.hpp"
#include "mc/prng/prng.hpp"
#include <Kokkos_Core_fwd.hpp>
#include <mc/particles/data_holder.hpp>

using LocalConcentrationView =
    Kokkos::Subview<Kokkos::View<const double **>, int, decltype(Kokkos::ALL)>;
using ContributionView =
    Kokkos::View<double **, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace>;

using model_properties_t = std::variant<double, std::string>;

using model_properties_detail_t =
    std::unordered_map<std::string, model_properties_t>;

template <typename T>
concept ParticleModel = requires(T model,
                                 MC::ParticleDataHolder &p,
                                 double d_t,
                                 const LocalConcentrationView &concentration,
                                 ContributionView contrib,Kokkos::Random_XorShift64_Pool<> _rng,MC::KPRNG rng2) {
  { model.init(p,_rng) } -> std::same_as<void>;
  { model.update(d_t, p, concentration,rng2) } -> std::same_as<void>;
  { model.division(p) } -> std::same_as<T>;
  { model.contribution(p, contrib) } -> std::same_as<void>;
  { model.get_properties() } -> std::same_as<model_properties_detail_t>;
};



#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
class DefaultModel
{
public:
  KOKKOS_INLINE_FUNCTION void init(MC::ParticleDataHolder &p,Kokkos::Random_XorShift64_Pool<> _rng)
  {
    p.status = MC::CellStatus::IDLE;
  }

  KOKKOS_INLINE_FUNCTION void
  update(double d_t,
         MC::ParticleDataHolder &p,
         const LocalConcentrationView &concentration,MC::KPRNG _rng)
  {
  }

  KOKKOS_INLINE_FUNCTION DefaultModel division(MC::ParticleDataHolder & /*p*/)
  {
    // Division logic
    return {};
  }

  KOKKOS_INLINE_FUNCTION void contribution(MC::ParticleDataHolder &p,
                                           ContributionView contrib)
  {
    // Contribution logic
  }

  inline model_properties_detail_t get_properties(){
    return {};
  }
};

static_assert(ParticleModel<DefaultModel>, "Check default model");

#pragma clang diagnostic pop

#endif