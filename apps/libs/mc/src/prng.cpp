#include <mc/prng.hpp>

namespace MC {
  double uniform_double_rand(const double &min, const double &max)
  {
    static thread_local std::mt19937 *generator_double = nullptr;
    if (!generator_double)
    {
      std::hash<std::thread::id> hasher;
      generator_double =
          new std::mt19937(clock() + hasher(std::this_thread::get_id()));
    }
    std::uniform_real_distribution<> distribution(min, max);
    return distribution(*generator_double);
  }

  

  double double_unfiform()
  {
    return uniform_double_rand(0., 1.);
  }
}