#ifndef __CMODEL_PARAMETERS_HPP__
#define __CMODEL_PARAMETERS_HPP__

#include <concepts>
#include <vector>

template <typename T>
concept ModelParameter = requires(const T &obj) {
  {
    obj.as_array()
  } -> std::same_as<std::vector<double>>;
};

#endif //__CMODEL_PARAMETERS_HPP__