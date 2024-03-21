#ifndef __COMMON_TYPES_HPP__
#define __COMMON_TYPES_HPP__

#include <span>

#include <any>

template <typename T> using ArrayView = std::span<T>;

using NumberView = ArrayView<double>;



#endif //__COMMON_TYPES_HPP__