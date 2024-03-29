#ifndef __MC_PRNG_HPP__
#define __MC_PRNG_HPP__

#include <concepts>
#include <random>
#include <thread>
#include <time.h>

template <typename T>
concept IntegerType = requires(T n) {
  requires std::is_integral_v<T>;
  requires !std::is_same_v<bool, T>;
  requires std::is_arithmetic_v<decltype(n + 1)>;
  requires !std::is_pointer_v<T>;
};

namespace MC
{
  template <IntegerType T> T uniform_int_rand(const T &min, const T &max)
  {
    static thread_local std::mt19937 *generator_int = nullptr;
    if (!generator_int)
    {
      std::hash<std::thread::id> hasher;
      generator_int =
          new std::mt19937(clock() + hasher(std::this_thread::get_id()));
    }
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(*generator_int);
  }

  double uniform_double_rand(const double &min, const double &max);

  

  double double_unfiform();

} // namespace MC

#endif //__MC_PRNG_HPP__