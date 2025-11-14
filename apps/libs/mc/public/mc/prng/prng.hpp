#ifndef __MC_PRNG_HPP__
#define __MC_PRNG_HPP__

#include "Kokkos_Macros.hpp"
#include <Kokkos_Random.hpp>
#include <common/common.hpp>
#include <common/traits.hpp>
#include <cstdint>
#include <mc/alias.hpp>
namespace MC
{

  using pool_type =
      Kokkos::Random_XorShift1024_Pool<Kokkos::DefaultExecutionSpace>;
  pool_type get_pool(std::size_t seed = 0);

  /** @brief Sample variable with RAI mecanism*/
  KOKKOS_INLINE_FUNCTION auto sample_random_variables(const pool_type& p,
                                                      auto functor)
  {
    auto gen = p.get_state();
    auto t = functor(gen);
    p.free_state(gen);
    return t;
  }

#define SAMPLE_RANDOM_VARIABLES(_random_pool_, ...)                            \
  auto _generator_state_ = _random_pool_.get_state();                          \
  __VA_ARGS__;                                                                 \
  _random_pool_.free_state(_generator_state_);

  /**
  @brief Utilities and wrap around kokkos random generator
   */
  class KPRNG // TODO remove deprecated
  {
  public:
    using pool_type = pool_type;
    using generator_type = pool_type::generator_type;

    explicit KPRNG(size_t _seed = 0);

    template <FloatingPointType T> KOKKOS_INLINE_FUNCTION T uniform() const
    {
      return this->uniform<T>(0, 1);
    }

    template <FloatingPointType T>
    KOKKOS_INLINE_FUNCTION T uniform(T a, T b) const
    {
      return sample_random_variables(
          random_pool,
          [a, b](auto gen)
          {
            T x;
            if constexpr (std::is_same_v<T, float>)
            {
              x = gen.frand(a, b);
            }
            else if constexpr (std::is_same_v<T, double>)
            {
              x = gen.drand(a, b);
            }
            return x;
          });
    }

    [[deprecated]] [[nodiscard]] Kokkos::View<double*, ComputeSpace>
    double_uniform(size_t n_sample, double a = 0., double b = 1.) const;

    template <FloatingPointType T, size_t n_r>
    Kokkos::View<T[n_r], ComputeSpace> random_view() const
    {
      Kokkos::View<T[n_r], ComputeSpace> A("random");
      Kokkos::fill_random(A, random_pool, 0., 1.);
      return A;
    }

    template <size_t n_r>
    [[deprecated]] KOKKOS_INLINE_FUNCTION std::array<double, n_r>
    double_uniform() const
    {
      return generate_uniform_impl<pool_type, n_r>(
          random_pool, std::make_index_sequence<n_r>{});
    }

    [[nodiscard]] KOKKOS_INLINE_FUNCTION uint64_t uniform_u(uint64_t a,
                                                            uint64_t b) const
    {

      return sample_random_variables(
          random_pool, [a, b](auto gen) { return gen.urand64(a, b); });
    }

    pool_type random_pool;

    [[nodiscard]] auto get_seed() const
    {
      return seed;
    } // TODO export this to result file

  private:
    std::size_t seed{};
    template <typename random_pool_t, size_t n_r, size_t... I>
    KOKKOS_INLINE_FUNCTION std::array<double, n_r>
    generate_uniform_impl(random_pool_t pool,
                          std::index_sequence<I...> /*unused*/) const
    {
      // Constexpr loopunrolling to fill the array
      auto generator = pool.get_state();
      std::array<double, n_r> res = {
          {(static_cast<void>(I), generator.drand(0., 1.))...}};
      pool.free_state(generator);
      return res;
    }
  };

} // namespace MC

#endif //__MC_PRNG_HPP__
