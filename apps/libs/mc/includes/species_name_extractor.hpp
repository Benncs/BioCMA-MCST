#ifndef __SPECIES_NAME_EXTRACTOR_HPP__
#define __SPECIES_NAME_EXTRACTOR_HPP__

#include <mc/traits.hpp>

namespace impl
{

  inline auto
  iota_species(std::size_t n_c)
  {
    std::vector<std::string> v(n_c);

    std::generate(
        v.begin(), v.end(), [n = 0]() mutable { return std::to_string(n++); });

    return v;
  }

  template <ModelType T>
    requires(!has_species_name<T>)
  auto
  get_species_names_impl()
  {
    return iota_species(T::n_c);
  }

  template <ModelType T>
    requires has_species_name<T>
  auto
  get_species_names_impl()
  {
    auto rd = T::species();

    if (rd.size() == 0)
    {
      return iota_species(T::n_c);
    }

    if (rd.size() != T::n_c && rd.size() != 0)
    {
      throw std::invalid_argument("Missing species names");
    }

    return std::vector<std::string>(rd.begin(), rd.end());
  }
} // namespace impl

#endif
