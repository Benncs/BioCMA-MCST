#ifndef __BIO__UTILS_HPP__
#define __BIO__UTILS_HPP__

#include "Kokkos_Macros.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_Random.hpp>
#include <common/static_for.hpp>
#include <common/traits.hpp>
#include <cstddef>
#include <mc/alias.hpp>
#include <mc/macros.hpp>
#include <mc/prng/prng.hpp>
#include <tuple>

/**
  @brief Models definition
 */
namespace Models
{
  template <typename T>
  KOKKOS_INLINE_FUNCTION T
  get_clamped_concentration_cast(const MC::LocalConcentration& c,
                                 const std::size_t position_index,
                                 const std::size_t species_index)
  {
    return static_cast<T>(Kokkos::max(
        static_cast<typename MC::LocalConcentration::non_const_value_type>(0),
        GET_CONCENTRATION(species_index)));
  }

#define GET_CLAMPED_CONCENTRATION_CAST(__T__, __species_index__)               \
  get_clamped_concentration_cast<__T__>(c, position_index, __species_index__)

#define GET_CLAMPED_CONCENTRATION(__species_index__)                           \
  get_clamped_concentration_cast<MC::LocalConcentration::const_value_type>(    \
      c, position_index, __species_index__)

  // template<FloatingPointType T>
  // struct partiton_t
  // {
  //   T a;
  //   T b;
  // };

  template <FloatingPointType F>
  static F constexpr _get_phi_s_max(F density,
                                    F dl,
                                    F glucose_to_biomass_yield = 0.5)
  {
    // dl and density must be same unit, dl*density -> mass and y is mass yield
    return (dl * density) / glucose_to_biomass_yield;
  }

  template <FloatingPointType T>
  KOKKOS_INLINE_FUNCTION Kokkos::pair<T, T>
  g_partition(T a, T p)
  {
    const auto a1 = a * p;
    const auto a2 = (static_cast<T>(1) - p) * a;
    return Kokkos::pair(a1, a2);
  }

  template <FloatingPointType T>
  KOKKOS_INLINE_FUNCTION MC::Status
  check_div(const T l, const T lc)
  {
    return (l >= lc) ? MC::Status::Division : MC::Status::Idle;
  }

  namespace MolarMass
  {
    namespace GramPerMole
    {

      template <FloatingPointType T> constexpr T glucose = static_cast<T>(180);
      template <FloatingPointType T> constexpr T dioxygen = static_cast<T>(32);
      template <FloatingPointType T> constexpr T acetate = static_cast<T>(59);

      constexpr double co2 = 44;
      ;
    }; // namespace GramPerMole

    constexpr double glucose = GramPerMole::glucose<double> * 1e-3;
    constexpr double dioxygen = GramPerMole::dioxygen<double> * 1e-3;
    ;
    constexpr double acetate = GramPerMole::acetate<double> * 1e-3;
    ;
    constexpr double co2 = GramPerMole::co2 * 1e-3;
    ;
    [[deprecated]] constexpr double X = 113.1e-3;
  }; // namespace MolarMass

  template <FloatingPointType F>
  KOKKOS_INLINE_FUNCTION consteval F
  c_linear_density(F rho, F d)
  {
    return rho * F(Kokkos::numbers::pi) * d * d / F(4.);
  }

  static constexpr double tau_division_proba = 1e-7;

  KOKKOS_INLINE_FUNCTION bool
  check_probability_division(double d_t,
                             double gamma,
                             MC::pool_type random_pool)
  {
    (void)d_t;
    // const double proba_div =
    //     (1 - Kokkos::exp(-d_t / tau_division_proba)) * gamma;

    auto generator = random_pool.get_state();
    const double x = generator.drand(0., 1.);
    random_pool.free_state(generator);

    return x < gamma;
  }

  namespace GammaDivision
  {
    KOKKOS_INLINE_FUNCTION double
    threshold_linear(double lenght, double threshold, double upper_bound)
    {
      return (lenght <= threshold)
                 ? 0
                 : (lenght - threshold) / (upper_bound - threshold);
    }

    KOKKOS_INLINE_FUNCTION constexpr float
    logistic(float x, float xmax, float alpha)
    {
      // auto blna = alpha * std::log(xmax);
      // return 1 / (1 + std::exp(-alpha * std::log(x) + blna));
      // return x<xmax?1. / (1 + std::exp(alpha * std::log((x-xmax)/xmax ))):1.;

      const auto z = x < xmax ? std::pow(x / (xmax - x), alpha) : 1.;
      return z / (1. + z);
    }

  } // namespace GammaDivision

  KOKKOS_INLINE_FUNCTION bool
  almost_equal(double val, double val2, double tolerance = 1e-8)
  {
    return Kokkos::abs(val - val2) < tolerance;
  }

  KOKKOS_INLINE_FUNCTION constexpr bool
  almost_zero(double x, double tol = 1e-8)
  {
    return (-tol < x) && (x < tol);
  }

  template <typename... Args>
  KOKKOS_INLINE_FUNCTION double
  min_var(Args... args)
  {
    return (Kokkos::min)({ args... });
  }
  template <>
  KOKKOS_INLINE_FUNCTION double
  min_var()
  {
    return 0.;
  }

  namespace impl
  {
    /// 1/x[i] for every i with a single division
    /// (https://medium.com/eryxcoop/montgomerys-trick-for-batch-galois-field-inversion-9b6d0f399da2)
    ///
    /// 3(N-1) multiplications + 1 division instead of N divisions.
    template <typename T, std::size_t N>
    KOKKOS_INLINE_FUNCTION constexpr Kokkos::Array<T, N>
    batch_reciprocal(const Kokkos::Array<T, N>& x)
    {
      static_assert(N > 0, "batch_reciprocal needs at least one value");

      Kokkos::Array<T, N> prefix{};
      prefix[0] = x[0];
      Common::static_for<N - 1>(
          [&](auto k)
          {
            constexpr std::size_t i = decltype(k)::value + 1;
            prefix[i] = prefix[i - 1] * x[i];
          });

      T inv = static_cast<T>(1) / prefix[N - 1];

      Kokkos::Array<T, N> out{};
      Common::static_for<N - 1>(
          [&](auto k)
          {
            constexpr std::size_t i = N - 1 - decltype(k)::value;
            out[i] = inv * prefix[i - 1];
            inv *= x[i];
          });
      out[0] = inv;

      return out;
    }
  } // namespace impl

  /**
   * @brief Monod factors c_i / (c_i + k_i) of several substrates, with a
   * single division for all of them.
   *
   * Arguments are (c, k) pairs, read the result with structured bindings:
   * @code
   * const auto [m_s, m_o2] = monod_batch(c_s, k_s, c_o2, k_o);
   * @endcode
   *
   * @warning c must be clamped >= 0: every denominator is then >= its k > 0
   * and the product of all of them stays finite and non zero.
   */
  template <typename... Ts>
  KOKKOS_INLINE_FUNCTION constexpr auto
  monod_batch(Ts... args)
  {
    using T = std::common_type_t<Ts...>;
    static_assert(FloatingPointType<T>, "monod_batch: floating point only");
    static_assert(sizeof...(Ts) > 0 && sizeof...(Ts) % 2 == 0,
                  "monod_batch: arguments are (c, k) pairs");
    constexpr std::size_t n = sizeof...(Ts) / 2;

    const Kokkos::Array<T, 2 * n> ck = { { static_cast<T>(args)... } };

    Kokkos::Array<T, n> denominators{};
    Common::static_for<n>(
        [&](auto i)
        {
          constexpr std::size_t j = decltype(i)::value;
          denominators[j] = ck[2 * j] + ck[(2 * j) + 1];
        });

    const auto inv = impl::batch_reciprocal(denominators);

    Kokkos::Array<T, n> m{};
    Common::static_for<n>(
        [&](auto i)
        {
          constexpr std::size_t j = decltype(i)::value;
          m[j] = ck[2 * j] * inv[j];
        });
    return m;
  }

  /**
   * @brief Production rate of every row of a stoichiometric matrix:
   * Stoichiometry . r.
   *
   * The matrix is a template parameter (rows x reactions, Kokkos::Array of
   * Kokkos::Array), so each entry is a constant expression: the product is
   * unrolled by Common::static_for and zero entries are skipped at compile
   * time, no multiplication by 0 is emitted.
   * @code
   * const auto d = apply_stoichiometry<stoichiometry>(r);
   * @endcode
   */
  template <auto Stoichiometry, typename Rates>
  KOKKOS_INLINE_FUNCTION constexpr auto
  apply_stoichiometry(const Rates& r)
  {
    // clang sees the template parameter object as a const reference
    using matrix_type = std::remove_cvref_t<decltype(Stoichiometry)>;
    using row_type = typename matrix_type::value_type;
    using T = typename row_type::value_type;
    constexpr std::size_t n_rows = matrix_type::size();
    constexpr std::size_t n_reactions = row_type::size();
    static_assert(n_reactions == Rates::size(),
                  "apply_stoichiometry: one rate per matrix column");

    Kokkos::Array<T, n_rows> d{};
    Common::static_for<n_rows>(
        [&](auto j)
        {
          Common::static_for<n_reactions>(
              [&](auto k)
              {
                // Local copies: nvcc rejects integral_constant::value used as
                // a runtime subscript ("undefined in device code")
                constexpr std::size_t row = decltype(j)::value;
                constexpr std::size_t col = decltype(k)::value;
                constexpr T n_jk = Stoichiometry[row][col];
                if constexpr (n_jk != 0)
                {
                  d[row] += n_jk * r[col];
                }
              });
        });
    return d;
  }

} // namespace Models

#endif
