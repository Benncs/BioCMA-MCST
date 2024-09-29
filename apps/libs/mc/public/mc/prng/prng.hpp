#ifndef __MC_PRNG_HPP__
#define __MC_PRNG_HPP__

#include <Kokkos_Random.hpp>
#include <common/kokkos_vector.hpp>
#include <cstdint>

// inline unsigned tau_step(unsigned &z, int S1, int S2, int S3, unsigned M)
// {
//   unsigned b = (((z << S1) ^ z) >> S2);
//   return z = (((z & M) << S3) ^ b);
// }

// inline unsigned LCGStep(unsigned &z, unsigned A, unsigned C)
// {
//   return z = (A * z + C);
// }

namespace MC
{

  class KPRNG
  {
  public:
    using SharedKPRNG = Kokkos::View<KPRNG, ComputeSpace>;
    explicit KPRNG(size_t _seed = 0);

    [[nodiscard]] KOKKOS_INLINE_FUNCTION double double_uniform() const
    {
      auto generator = random_pool.get_state();
      double x = generator.drand(0., 1.);
      random_pool.free_state(generator);
      return x;
    }

    KOKKOS_INLINE_FUNCTION double double_uniform(double a, double b) const
    {
      auto generator = random_pool.get_state();
      double x = generator.drand(a, b);
      random_pool.free_state(generator);
      return x;
    }

    [[nodiscard]] Kokkos::View<double *, ComputeSpace>
    double_uniform(size_t n_sample, double a = 0., double b = 1.) const;

    template <size_t n_r>
    Kokkos::View<double[n_r], ComputeSpace> double_uniform_view() const
    {
      Kokkos::View<double[n_r], ComputeSpace> A("random");
      Kokkos::fill_random(A, random_pool, 0., 1.);
      return A;
    }

    template <size_t n_r>
    KOKKOS_INLINE_FUNCTION std::array<double, n_r> double_uniform() const
    {
      return generate_uniform_impl<Kokkos::Random_XorShift64_Pool<>, n_r>(
          random_pool, std::make_index_sequence<n_r>{});
    }

    [[nodiscard]] KOKKOS_INLINE_FUNCTION uint64_t uniform_u(uint64_t a,
                                                            uint64_t b) const
    {

      auto generator = random_pool.get_state();
      uint64_t x = generator.urand64(a, b);
      random_pool.free_state(generator);
      return x;
    }

    
    Kokkos::Random_XorShift64_Pool<> random_pool{};

  private:
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