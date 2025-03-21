#ifndef __MC_PRNG_HPP__
#define __MC_PRNG_HPP__

#include <Kokkos_Random.hpp>
#include <common/common.hpp>
#include <common/traits.hpp>
#include <cstdint>

namespace MC
{
  /**
  @brief Utilities and wrap around kokkos random generator
   */
  class KPRNG //TODO remove deprecated 
  {
  public:
    using pool_type = Kokkos::Random_XorShift1024_Pool<Kokkos::DefaultExecutionSpace>;
    explicit KPRNG(size_t _seed = 0);

    [[deprecated]] [[nodiscard]] KOKKOS_INLINE_FUNCTION double double_uniform() const
    {
      auto generator = random_pool.get_state();
      double x = generator.drand(0., 1.);
      random_pool.free_state(generator);
      return x;
    }

    template <FloatingPointType T> KOKKOS_INLINE_FUNCTION T uniform() const
    {
      auto generator = random_pool.get_state();
      T x;
      if constexpr (std::is_same_v<T, float>)
      {
        x = generator.frand();
      }
      else if constexpr (std::is_same_v<T, double>)
      {
        x = generator.drand();
      }
      random_pool.free_state(generator);
      return x;
    }

    template <FloatingPointType T> KOKKOS_INLINE_FUNCTION T uniform(T a, T b) const
    {
      auto generator = random_pool.get_state();
      T x;
      if constexpr (std::is_same_v<T, float>)
      {
        x = generator.frand(a, b);
      }
      else if constexpr (std::is_same_v<T, double>)
      {
        x = generator.drand(a, b);
      }
      random_pool.free_state(generator);
      return x;
    }

    [[deprecated]] [[nodiscard]] Kokkos::View<double*, ComputeSpace>
    double_uniform(size_t n_sample, double a = 0., double b = 1.) const;

    template <size_t n_r> Kokkos::View<double[n_r], ComputeSpace> double_uniform_view() const
    {
      Kokkos::View<double[n_r], ComputeSpace> A("random");
      Kokkos::fill_random(A, random_pool, 0., 1.);
      return A;
    }

    template <size_t n_r>
    [[deprecated]] KOKKOS_INLINE_FUNCTION std::array<double, n_r> double_uniform() const
    {
      return generate_uniform_impl<pool_type, n_r>(random_pool, std::make_index_sequence<n_r>{});
    }

    [[nodiscard]] KOKKOS_INLINE_FUNCTION uint64_t uniform_u(uint64_t a, uint64_t b) const
    {

      auto generator = random_pool.get_state();
      uint64_t x = generator.urand64(a, b);
      random_pool.free_state(generator);
      return x;
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
    generate_uniform_impl(random_pool_t pool, std::index_sequence<I...> /*unused*/) const
    {
      // Constexpr loopunrolling to fill the array
      auto generator = pool.get_state();
      std::array<double, n_r> res = {{(static_cast<void>(I), generator.drand(0., 1.))...}};
      pool.free_state(generator);
      return res;
    }
  };

} // namespace MC

#endif //__MC_PRNG_HPP__