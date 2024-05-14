#ifndef __MC_PRNG_HPP__
#define __MC_PRNG_HPP__

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
  static constexpr size_t MC_RAND_DEFAULT_SEED = 10;
  class PRNG
  {
  public:

    explicit PRNG()
    {
      gen = std::mt19937_64(MC_RAND_DEFAULT_SEED);
    }

    explicit PRNG(uint64_t seed)
    {
      gen = std::mt19937_64(seed);
    }

    template <typename T> T uniform_int_rand(const T &min, const T &max)
    {
      T rn {};
      std::uniform_int_distribution<int> int_uniform(min, max);
      rn =  int_uniform(gen);
      return rn;
    }

    double uniform_double_rand(double min, double max)
    {
      double rn =0;
      std::uniform_real_distribution<double> double_distribution(min, max);
      rn= double_distribution(gen);
      return rn;
    }

    inline double double_unfiform()
    {
      return uniform_double_rand(0., 1.);
    }

  private:
    std::mt19937_64 gen;
  };

  // template <IntegerType T> T uniform_int_rand(const T &min, const T &max)
  // {
  //   static thread_local std::mt19937 *generator_int = nullptr;
  //   if (!generator_int)
  //   {
  //     std::hash<std::thread::id> hasher;
  //     generator_int =
  //         new std::mt19937(clock() + hasher(std::this_thread::get_id()));
  //   }
  //   std::uniform_int_distribution<int> distribution(min, max);
  //   return distribution(*generator_int);
  // }

  // inline double uniform_double_rand(double min, double max)
  // {
  //   static thread_local std::mt19937 *generator_double = nullptr;
  //   if (generator_double == nullptr)
  //   {
  //     std::hash<std::thread::id> hasher;
  //     generator_double =
  //         new std::mt19937(clock() + hasher(std::this_thread::get_id()));
  //   }
  //   std::uniform_real_distribution<> distribution(min, max);
  //   return distribution(*generator_double);
  // }

  // inline double normal_double_rand(double mean, double stddev)
  // {
  //   static thread_local std::mt19937 *generator_double = nullptr;
  //   if (generator_double == nullptr)
  //   {
  //     std::hash<std::thread::id> hasher;
  //     generator_double =
  //         new std::mt19937(clock() + hasher(std::this_thread::get_id()));
  //   }
  //   std::normal_distribution<> distribution(mean, stddev);
  //   return distribution(*generator_double);
  // }

  // inline double double_unfiform()
  // {
  //   return uniform_double_rand(0., 1.);
  // }

} // namespace MC

#endif //__MC_PRNG_HPP__