#include "common/kokkos_vector.hpp"
#include <mc/prng/prng.hpp>

#ifdef NDEBUG
#  include <random>
#endif

static constexpr size_t MC_RAND_DEFAULT_SEED = 2024;
namespace MC
{
  KPRNG::KPRNG(size_t _seed)
  {
    if (_seed == 0)
    {
#ifndef NDEBUG
      const size_t seed = MC_RAND_DEFAULT_SEED;
#else
      const size_t seed = std::random_device{}();
#endif
      random_pool = Kokkos::Random_XorShift64_Pool<>(seed);
    }
    else
    {
      random_pool = Kokkos::Random_XorShift64_Pool<>(_seed);
    }
  };

  Kokkos::View<double *,ComputeSpace>
  KPRNG::double_uniform(size_t n_sample, double a, double b) const
  {
    Kokkos::View<double *,ComputeSpace> view("sample", n_sample);
    auto local_pool = random_pool;
    Kokkos::parallel_for(
        n_sample, KOKKOS_LAMBDA(auto &&i) {
          auto generator = local_pool.get_state();
          view(i) = generator.drand(a, b);
          local_pool.free_state(generator);
        });
    Kokkos::fence();

    return view;
  }
} // namespace MC