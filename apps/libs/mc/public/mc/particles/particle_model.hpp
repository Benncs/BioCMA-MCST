#ifndef __MC_PARTICLE_MODEL_HPP__
#define __MC_PARTICLE_MODEL_HPP__

#include <Kokkos_Core.hpp>
#include <common/kokkos_vector.hpp>
#include <mc/particles/data_holder.hpp>
#include <mc/prng/prng.hpp>
#include <type_traits>
#include <variant>

#ifdef __GNUC__
// Avoid tons of warnings with root code
#  pragma GCC system_header
#endif

using LocalConcentrationView =
    Kokkos::Subview<Kokkos::View<const double **>,
                    int,
                    decltype(Kokkos::ALL)>; ///< Concentration inside the
                                            ///< current container that particle
                                            ///< can access

using ContributionView =
    Kokkos::View<double **,
                 Kokkos::LayoutLeft,
                 ComputeSpace>; ///< Contribution inside the particle's current
                                ///< container

using model_properties_t =
    std::variant<double,
                 std::string>; ///< Type of properties that model can export

using model_properties_detail_t =
    std::unordered_map<std::string,
                       model_properties_t>; ///< Type of properties data that
                                            ///< model can export

/**
 * @brief Concept that defines the requirements for a Biological Model.
 *
 * This concept specifies the interface and requirements that the class must
 * implement to be a Biological Model
 *
 * - `init(MC::ParticleDataHolder &p, MC::KPRNG rng)`: Initializes the model
 *   with the given particle data holder and random number generator.
 *
 * - `update(double d_t, MC::ParticleDataHolder &p, const LocalConcentrationView
 * &concentration, MC::KPRNG rng)`: Updates the model's state based on the
 * elapsed time, particle data, concentration view, and random number generator.
 *
 * - `division(MC::ParticleDataHolder &p)`: Returns a new instance of the model
 * representing the result of a division operation.
 *
 * - `contribution(MC::ParticleDataHolder &p, ContributionView contrib)`:
 * Calculates and applies contributions to the model based on the particle data
 * holder and contribution view.
 *
 * - `get_properties()`: Returns the properties of the model.
 *
 * @tparam T The type of the model that must satisfy the concept requirements.
 *
 * @details
 * - The `init` method must not return any value (void).
 * - The `update` method must not return any value (void).
 * - The `division` method must return a new instance of type `T`.
 * - The `contribution` method must not return any value (void).
 * - The `get_properties` method must return an object of type
 * `model_properties_detail_t`.
 *
 * @see MC::ParticleDataHolder, LocalConcentrationView, ContributionView,
 * MC::KPRNG, model_properties_detail_t
 */
template <typename T>
concept ParticleModel = requires(T model,
                                 std::is_default_constructible<T>,
                                 MC::ParticleDataHolder &p,
                                 double d_t,
                                 const LocalConcentrationView &concentration,
                                 ContributionView contrib,
                                 MC::KPRNG rng) {
  { model.init(p, rng) } -> std::same_as<void>;
  { model.update(d_t, p, concentration, rng) } -> std::same_as<void>;
  { model.division(p) } -> std::same_as<T>;
  { model.contribution(p, contrib) } -> std::same_as<void>;
  { model.get_properties() } -> std::same_as<model_properties_detail_t>;
};

/**
 * @brief Default Model that implement ParticleModel
 * This model does nothing but is required to have correct compilation even if
 * other models were not provided
 */
class DefaultModel
{
public:
  KOKKOS_INLINE_FUNCTION void init(MC::ParticleDataHolder &p, MC::KPRNG _rng)
  {
    p.status = MC::CellStatus::IDLE;
  }

  KOKKOS_INLINE_FUNCTION void
  update(double d_t,
         MC::ParticleDataHolder &p,
         const LocalConcentrationView &concentration,
         MC::KPRNG _rng)
  {
  }

  KOKKOS_INLINE_FUNCTION DefaultModel division(MC::ParticleDataHolder & /*p*/)
  {
    return {};
  }

  KOKKOS_INLINE_FUNCTION void contribution(MC::ParticleDataHolder &p,
                                           ContributionView contrib)
  {
  }

  inline model_properties_detail_t get_properties()
  {
    return {};
  }
};

static_assert(ParticleModel<DefaultModel>, "Check default model");

#endif