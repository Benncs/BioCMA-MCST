#ifndef __SPECIES_NAME_EXTRACTOR_HPP__
#define __SPECIES_NAME_EXTRACTOR_HPP__

#include <mc/traits.hpp>

namespace impl
{
  template <ModelType T>
    requires(!has_species_name<T>)
  auto
  get_species_names_impl()
  {
    std::vector<std::string> v(T::n_c);

    std::generate(
        v.begin(), v.end(), [n = 0]() mutable { return std::to_string(n++); });

    return v;
  }

  template <ModelType T>
    requires has_species_name<T>
  auto
  get_species_names_impl()
  {
    auto rd = T::species();

    if (rd.size() != T::n_c)
    {
      throw std::invalid_argument("Missing species names");
    }

    return std::vector<std::string>(rd.begin(), rd.end());
  }
} // namespace impl

#endif
