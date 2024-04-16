#ifndef __CMODEL_PARAMETERS_HPP__
#define __CMODEL_PARAMETERS_HPP__

#include <concepts>
#include <vector>

template <typename T>
concept ModelParameter = requires(T obj, double d) {
  {
    obj.step(d)
  } -> std::same_as<void>;
};

#endif //__CMODEL_PARAMETERS_HPP__