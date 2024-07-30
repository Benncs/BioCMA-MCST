#ifndef __MC_DISTRIBUTION_HPP__
#define __MC_DISTRIBUTION_HPP__

#include <common/common.hpp>
#include <cstdint>
#include <functional>
#include <random>
#include <type_traits>
#include <variant>

namespace MC
{

  struct ExponentialLaw
  {
    double lambda;
  };

  struct BoundedExponentialLaw
  {
    double lambda;
    double min;
    double max;
  };

  // Struct representing Normal distribution
  struct NormalLaw
  {
    double mean;
    double stddev;
  };

  // Struct representing Uniform distribution
  struct UniformLaw
  {
    double lower;
    double upper;
  };

  struct UniformLawINT
  {
    int64_t lower;
    int64_t upper;
  };

  using DistributionVariantReal =
      std::variant<ExponentialLaw, NormalLaw, UniformLaw>;
  using DistributionVariantInt =
      std::variant<ExponentialLaw, BoundedExponentialLaw, UniformLawINT>;

  template <class... Ts> struct overloaded : Ts...
  {
    using Ts::operator()...;
  };
  
  // explicit deduction guide (not needed as of C++20)
  // template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

  template <typename Val, typename... Ts> auto match(Val &&val, Ts... ts)
  {
    return std::visit(overloaded{ts...}, std::forward<Val>(val));
  }

  template <NumberType T>
  using distribution_t = std::function<T(std::mt19937 &)>;

  template <NumberType T> T bounded_exponential(std::mt19937 &g, T min, T max);


  template <NumberType T>
  distribution_t<T> get_distribution(DistributionVariantReal params)
  {
    auto exp_func = [&](const ExponentialLaw &law) -> distribution_t<T>
    {
      std::exponential_distribution<T> dist(law.lambda);
      auto l = [_dist = dist](std::mt19937 &gen) mutable { return _dist(gen); };
      return l;
    };

    auto normal_func = [&](const NormalLaw &law) -> distribution_t<T>
    {
      std::normal_distribution<T> dist(law.mean, law.stddev);
      auto l = [_dist = dist](std::mt19937 &gen) mutable { return _dist(gen); };
      return l;
    };

    auto uniform_func = [&](const UniformLaw &law) -> distribution_t<T>
    {
      std::uniform_real_distribution<T> dist(law.lower, law.upper);
      auto l = [_dist = dist](std::mt19937 &gen) mutable { return _dist(gen); };
      return l;
    };

    return match(params, exp_func, normal_func, uniform_func);
  }

  template <typename IntegerType>
  distribution_t<IntegerType>
  get_distribution_int(DistributionVariantInt params)
  {
    auto exp_func =
        [&](const ExponentialLaw &law) -> distribution_t<IntegerType>
    {
      std::exponential_distribution<> dist(law.lambda);
      auto l = [_dist = dist](std::mt19937 &gen) mutable { return _dist(gen); };
      return l;
    };

    auto b_exp_func =
        [&](const BoundedExponentialLaw &law) -> distribution_t<IntegerType>
    {
      std::exponential_distribution<> dist(law.lambda);
      auto l = [_dist = dist, min = law.min, max = law.max](
                   std::mt19937 &gen) mutable
      {
        IntegerType val = _dist(gen);
        while (val < min || val > max)
        {
          val = _dist(gen);
        }
        return val;
      };
      return l;
    };

    auto uniform_func =
        [&](const UniformLawINT &law) -> distribution_t<IntegerType>
    {
      std::uniform_int_distribution<IntegerType> dist(law.lower, law.upper);
      auto l = [_dist = dist](std::mt19937 &gen) mutable { return _dist(gen); };
      return l;
    };

    return match(params, b_exp_func, exp_func, uniform_func);
  }
} // namespace MC

#endif //__MC_DISTRIBUTION_HPP__