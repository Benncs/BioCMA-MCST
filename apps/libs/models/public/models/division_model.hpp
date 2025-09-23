#ifndef __MODEL_DIVISION_HPP__
#define __MODEL_DIVISION_HPP__

#include <common/traits.hpp>
#include <mc/traits.hpp>
namespace Models
{

  namespace DivisionModel
  {

    template <typename T>
    concept DivisionModelModule = requires(T model) {
      { T::do_div() } -> std::same_as<bool>;
    };

    template <ModelType M> struct Adder
    {
      MODEL_CONSTANT M::FloatType v_adder = 1.5e-6; // m
    };
  }; // namespace DivisionModel

} // namespace Models

#endif
