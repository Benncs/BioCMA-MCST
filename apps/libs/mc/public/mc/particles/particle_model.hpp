// #ifndef __MC_PARTICLE_MODEL_HPP__
// #define __MC_PARTICLE_MODEL_HPP__

// #include <common/kokkos_vector.hpp>
// #include <Kokkos_Core.hpp>
// #include <Kokkos_ScatterView.hpp>
// #include <common/common.hpp>
// #include <concepts>
// #include <mc/particles/data_holder.hpp>
// #include <mc/prng/prng.hpp>
// #include <type_traits>

// WARN_EXPERIMENTAL

// #ifdef __GNUC__
// // Avoid tons of warnings with root code
// #  pragma GCC system_header
// #endif

// using LocalConcentrationView = Kokkos::Subview<Kokkos::View<const double**>,
//                                                int,
//                                                decltype(Kokkos::ALL)>; ///< Concentration inside the
//                                                                        ///< current container that
//                                                                        ///< particle can access

// using kernelContribution = Kokkos::View<double**, Kokkos::LayoutLeft, ComputeSpace>;

// using ContributionView =
//     Kokkos::Experimental::ScatterView<double**, Kokkos::LayoutLeft>; ///< Contribution inside the
//                                                                      ///< particle's current
//                                                                      ///< container
// template <typename MemorySpace>
// using ParticlePropertyViewType = Kokkos::View<double**, Kokkos::LayoutRight, MemorySpace>;

// using SubViewtype =
//     Kokkos::Subview<ParticlePropertyViewType<ComputeSpace>, decltype(Kokkos::ALL), std::size_t>;

// template <typename T>
// concept HasMass = requires(const T obj) {
//   {
//     obj.mass()
//   } -> std::convertible_to<double>; // mass() must be callable on const object, returning something
//                                     // convertible to double
// };

// // template <typename T>
// // concept HasNumberExportProperties = requires {
// //   { T::get_number() } -> std::same_as<std::size_t>;
// //   { std::bool_constant<(T::get_number(), true)>() } -> std::same_as<std::true_type>;
// // };

// template <typename T>
// concept HasNumberExportProperties = requires {
//   { T::get_number() } -> std::same_as<std::size_t>;
//   // { std::bool_constant<(T::get_number(), true)>() } -> std::same_as<std::true_type>;
// };

// template <typename T>
// concept HasExportProperties =
//     HasNumberExportProperties<T> && requires(const T obj, SubViewtype subview) {
//       { T::names() } -> std::convertible_to<std::vector<std::string>>;
//       { obj.fill_properties(subview) } -> std::same_as<void>;
//     };

// /**
//  * @brief Concept that defines the requirements for a Biological Model.
//  *
//  * This concept specifies the interface and requirements that the class must
//  * implement to be a Biological Model
//  *
//  * - `init(MC::ParticleDataHolder &p, MC::KPRNG rng)`: Initializes the model
//  *   with the given particle data holder and random number generator.
//  *
//  * - `update(double d_t, MC::ParticleDataHolder &p, const LocalConcentrationView
//  * &concentration, MC::KPRNG rng)`: Updates the model's state based on the
//  * elapsed time, particle data, concentration view, and random number generator.
//  *
//  * - `division(MC::ParticleDataHolder &p)`: Returns a new instance of the model
//  * representing the result of a division operation.
//  *
//  * - `contribution(MC::ParticleDataHolder &p, ContributionView contrib)`:
//  * Calculates and applies contributions to the model based on the particle data
//  * holder and contribution view.
//  *
//  * - `get_properties()`: Returns the properties of the model.
//  *
//  * @tparam T The type of the model that must satisfy the concept requirements.
//  *
//  * @details
//  * - The `init` method must not return any value (void).
//  * - The `update` method must not return any value (void).
//  * - The `division` method must return a new instance of type `T`.
//  * - The `contribution` method must not return any value (void).
//  * - The `get_properties` method must return an object of type
//  * `model_properties_detail_t`.
//  *
//  * @see MC::ParticleDataHolder, LocalConcentrationView, ContributionView,
//  * MC::KPRNG, model_properties_detail_t
//  */
// template <typename T>
// concept ParticleModel = HasExportProperties<T> && HasMass<T> &&
//                         requires(T model,
//                                  std::is_default_constructible<T>,
//                                  MC::ParticleDataHolder& p,
//                                  double d_t,
//                                  const LocalConcentrationView& concentration,
//                                  const ContributionView& contrib,
//                                  MC::KPRNG rng) {
//                           { model.init(p, rng) } -> std::same_as<void>;
//                           { model.update(d_t, p, concentration, rng) } -> std::same_as<void>;
//                           { model.division(p, rng) } -> std::same_as<T>;
//                           { model.contribution(p, contrib) } -> std::same_as<void>;

//                           { model.mass() } -> std::same_as<double>;
//                         };

// /**
//  * @brief Default Model that implement ParticleModel
//  * This model does nothing but is required to have correct compilation even if
//  * other models were not provided
//  */
// class DefaultModel
// {
// public:
//   KOKKOS_INLINE_FUNCTION void init(MC::ParticleDataHolder& p, MC::KPRNG _rng) noexcept
//   {
//     (void)_rng;
//     p.status = MC::CellStatus::IDLE;
//   }

//   KOKKOS_INLINE_FUNCTION void update(const double d_t,
//                                      MC::ParticleDataHolder& p,
//                                      const LocalConcentrationView& concentration,
//                                      MC::KPRNG _rng)
//   {

//     (void)d_t;
//     (void)p;
//     (void)concentration;
//     (void)_rng;
//   }

//   KOKKOS_INLINE_FUNCTION DefaultModel division(MC::ParticleDataHolder& p, MC::KPRNG k) noexcept
//   {
//     (void)p;
//     (void)k;
//     return {};
//   }

//   KOKKOS_INLINE_FUNCTION void contribution(MC::ParticleDataHolder& p,
//                                            const ContributionView& contrib) noexcept
//   {
//     (void)p;
//     (void)contrib;
//   }

//   // TODO: make this optinnal -> separation HasExportProperties and ParticleModel
//   static std::vector<std::string> names()
//   {
//     return {"mass"};
//   }

//   KOKKOS_FORCEINLINE_FUNCTION static consteval std::size_t get_number()
//   {
//     return 1;
//   }

//   KOKKOS_FORCEINLINE_FUNCTION void fill_properties(SubViewtype full) const
//   {
//     full(0) = mass();
//   }

//   [[nodiscard]] KOKKOS_INLINE_FUNCTION double mass() const noexcept
//   {
//     return 1.;
//   }
// };

// static_assert(ParticleModel<DefaultModel>, "Check default model");

// #endif
