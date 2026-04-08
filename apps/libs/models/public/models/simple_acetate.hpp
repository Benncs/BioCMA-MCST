
// #ifndef __SIMPLE_ACETATE_MODEL_HPP__
// #define __SIMPLE_ACETATE_MODEL_HPP__

// #include "Kokkos_Assert.hpp"
// #include "Kokkos_Core_fwd.hpp"
// #include "Kokkos_Macros.hpp"
// #include "common/common.hpp"
// #include "common/traits.hpp"
// #include "mc/macros.hpp"
// #include "models/utils.hpp"
// #include <mc/prng/prng_extension.hpp>
// #include <mc/traits.hpp>
// #include <optional>
// #include <string_view>

// namespace Models
// {
//   // template <FloatingPointType F>
//   // static F consteval
//   // _get_phi_s_max(F density, F dl, F glucose_to_biomass_yield = 0.5)
//   // {
//   //   // dl and density must be same unit, dl*density -> mass and y is mass
//   //   yield return (dl * density) / glucose_to_biomass_yield;
//   // }
//   struct SimpleAcetate
//   {
//     using uniform_weight = std::true_type;
//     using Self = SimpleAcetate;
//     using FloatType = float;
//     using Config = std::nullopt_t;

//     enum class particle_var : int // NOLINT
//     {
//       length = 0,
//       l_max,
//       a_p,
//       a_max,
//       // Export
//       a_e,
//       a_e_s,
//       a_e_a,
//       // Con
//       phi_s,
//       phi_a,
//       __COUNT__
//     };

//     static constexpr std::size_t n_var
//         = INDEX_FROM_ENUM(particle_var::__COUNT__);

//     static constexpr std::string_view name = "simple_acetate";
//     using SelfParticle = MC::ParticlesModel<Self::n_var, Self::FloatType>;

//     MODEL_CONSTANT std::size_t N_N = 2;              // Number of species
//     MODEL_CONSTANT FloatType a_max_m = 2e-6 / 3600.; // m
//     MODEL_CONSTANT FloatType l_max_m = 2e-6;         // m
//     MODEL_CONSTANT FloatType l_min_m = l_max_m / 2.; // m
//     // MODEL_CONSTANT FloatType k = 1e-3;                 // m
//     MODEL_CONSTANT FloatType d_m = 0.6e-6; // m
//     MODEL_CONSTANT FloatType lin_density
//         = c_linear_density(static_cast<FloatType>(1000), d_m);

//     MODEL_CONSTANT FloatType k_s = 1e-3;
//     MODEL_CONSTANT FloatType k_a = 1e-4;
//     MODEL_CONSTANT Kokkos::Array<FloatType, N_N> k = { k_s, k_a };

//     MODEL_CONSTANT FloatType y_s = 2;
//     MODEL_CONSTANT FloatType y_a = 3;
//     MODEL_CONSTANT Kokkos::Array<FloatType, N_N> y = { y_s, y_a };

//     MODEL_CONSTANT FloatType phi_s_max
//         = _get_phi_s_max<FloatType>(lin_density, a_max_m); // kgS/s

//     MODEL_CONSTANT auto l_max_dist
//         = MC::Distributions::TruncatedNormal<FloatType>(
//             l_max_m, l_max_m / 10., l_max_m * 0.7, 1.3 * l_max_m);

//     MODEL_CONSTANT auto l_dist =
//     MC::Distributions::TruncatedNormal<FloatType>(
//         l_max_m * 0.75, l_max_m / 10., 0.7 * l_min_m, l_max_m * 1.3);

//     MC::ContribIndexBounds static get_bounds();

//     // static Self::Config get_config(std::size_t n);

//     KOKKOS_INLINE_FUNCTION static void init(const MC::pool_type& random_pool,
//                                             std::size_t idx,
//                                             const SelfParticle& arr);

//     KOKKOS_INLINE_FUNCTION static MC::Status
//     update(const MC::pool_type& random_pool,
//            FloatType d_t,
//            std::size_t idx,
//            const SelfParticle& arr,
//            const MC::LocalConcentration& c);

//     KOKKOS_INLINE_FUNCTION static void
//     division(const MC::pool_type& random_pool,
//              std::size_t idx,
//              std::size_t idx2,
//              const SelfParticle& arr,
//              const SelfParticle& buffer_arr);

//     KOKKOS_INLINE_FUNCTION static double
//     mass(std::size_t idx, const SelfParticle& arr)
//     {
//       return GET_PROPERTY(Self::particle_var::length) * lin_density;
//     }

//     static std::array<std::string_view, Self::n_var>
//     names()
//     {
//       return {
//         "length", "l_max", "a_p",   "a_max", "a_e",
//         "a_e_s",  "a_e_a", "phi_s", "phi_a",
//       };
//     }

//     // static std::vector<std::size_t>
//     // get_number()
//     // {
//     //   return { INDEX_FROM_ENUM(particle_var::length) };
//     // }
//   };

//   CHECK_MODEL(SimpleAcetate)

//   KOKKOS_INLINE_FUNCTION void
//   SimpleAcetate::init([[maybe_unused]] const MC::pool_type& random_pool,
//                       std::size_t idx,
//                       const SelfParticle& arr)
//   {

//     MODEL_CONSTANT auto a_max_dist
//         = MC::Distributions::TruncatedNormal<FloatType>(
//             a_max_m, a_max_m / 2., 0.5 * a_max_m, a_max_m * 1.5);

//     static constexpr auto ld = l_dist;
//     static constexpr auto lm = l_dist;
//     auto gen = random_pool.get_state();
//     GET_PROPERTY(particle_var::length) = ld.draw(gen);
//     GET_PROPERTY(particle_var::l_max) = lm.draw(gen);
//     GET_PROPERTY(particle_var::a_p) = a_max_m / 2.;
//     GET_PROPERTY(particle_var::a_max) = a_max_dist.mean();
//     random_pool.free_state(gen);

//     GET_PROPERTY(particle_var::phi_s) = 0.;
//   }

//   KOKKOS_INLINE_FUNCTION MC::Status
//   SimpleAcetate::update([[maybe_unused]] const MC::pool_type& random_pool,
//                         FloatType d_t,
//                         std::size_t idx,
//                         const SelfParticle& arr,
//                         const MC::LocalConcentration& c)
//   {
//     Kokkos::Array<FloatType, N_N> adm
//         = { GET_PROPERTY(particle_var::a_max),
//             GET_PROPERTY(particle_var::a_max) / 3 };

//     Kokkos::Array<FloatType, N_N> D{};

//     FloatType inv = 1. / (c[0] + k[0]);
//     D[0] = adm[0] * c[0] * inv;
//     inv = 1. / (c[1] + k[1]);
//     D[1] = adm[1] * c[1] * inv;

//     GET_PROPERTY(particle_var::a_e) = 0.;
//     Kokkos::Array<FloatType, N_N> U{};
//     U[0] = Kokkos::min(D[0], GET_PROPERTY(particle_var::a_p));
//     GET_PROPERTY(particle_var::a_e) += U[0];

//     const FloatType pa = D[0] - GET_PROPERTY(particle_var::a_p);
//     const auto mask_pa = static_cast<FloatType>(pa < 0.F);

//     U[1] = mask_pa * Kokkos::min(D[1], -pa) + (1 - mask_pa) * 0.F;
//     GET_PROPERTY(particle_var::a_e) += U[1];

//     GET_PROPERTY(particle_var::length) += d_t *
//     GET_PROPERTY(particle_var::a_e);

//     GET_PROPERTY(particle_var::a_e_s) = U[0];
//     GET_PROPERTY(particle_var::a_e_a) = U[1];

//     GET_PROPERTY(particle_var::phi_s) = -1 * D[0] * lin_density * y[0];
//     GET_PROPERTY(particle_var::phi_a)
//         = mask_pa * (-U[1] * lin_density * y[1])
//           + (1.F - mask_pa) * (pa * lin_density * y[0] / y[1]);

//     return check_div(GET_PROPERTY(Self::particle_var::length),
//                      GET_PROPERTY(Self::particle_var::l_max));
//   }

//   KOKKOS_INLINE_FUNCTION void
//   SimpleAcetate::division([[maybe_unused]] const MC::pool_type& random_pool,
//                           std::size_t idx,
//                           std::size_t idx2,
//                           const SelfParticle& arr,
//                           const SelfParticle& buffer_arr)
//   {
//     Kokkos::View<FloatType**,
//                  Kokkos::LayoutRight,
//                  Kokkos::MemoryTraits<Kokkos::MemoryTraitsFlags::Restrict>>
//         buffer_e = buffer_arr;

//     const FloatType current_l = GET_PROPERTY(particle_var::length);
//     const FloatType new_current_length = current_l / 2.F;
//     GET_PROPERTY(particle_var::length) = new_current_length;
//     constexpr auto binf = INDEX_FROM_ENUM(particle_var::length);
//     constexpr auto bsup = INDEX_FROM_ENUM(particle_var::a_e);
//     for (auto i = binf; i < bsup; ++i)
//     {
//       COPY_PROPERTY_TO(i, idx2, buffer_e);
//     }

//     const auto current_a_e = GET_PROPERTY(particle_var::a_e);
//     auto gen = random_pool.get_state();

//     const double sigma = 0.2;
//     const double average = Kokkos::log(current_a_e) - sigma * sigma / 2;

//     const auto dist = MC::Distributions::LogNormal<double>(average, sigma);
//     const auto gen1 = static_cast<FloatType>(dist.draw(gen));
//     const auto gen2 = static_cast<FloatType>(dist.draw(gen));

//     static constexpr auto local_l = l_max_dist;
//     const FloatType lmax1 = local_l.draw(gen);
//     const FloatType lmax2 = local_l.draw(gen);

//     // KOKKOS_ASSERT(gen1 != 0 && gen2 != 0);

//     GET_PROPERTY(particle_var::a_p) = gen1;
//     GET_PROPERTY(particle_var::l_max) = lmax1;
//     GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::a_p) = gen2;
//     GET_PROPERTY_FROM(idx2, buffer_arr, particle_var::l_max) = lmax2;

//     random_pool.free_state(gen);
//   }

//   inline MC::ContribIndexBounds
//   SimpleAcetate::get_bounds()
//   {
//     int begin = INDEX_FROM_ENUM(Self::particle_var::phi_s);
//     return { .begin = begin, .end = begin + 2 };
//   }

// } // namespace Models

// #endif
