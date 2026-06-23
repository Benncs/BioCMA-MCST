#include <algorithm>
#include <mixture/species_descriptor.hpp>
#include <optional>

namespace
{

  constexpr Mixture::Specie Glucose = {
    "glucose",
    180.,
    std::nullopt,
  };

  constexpr Mixture::Specie Oxygen = {
    "oxygen",
    32.,
    3.181e-2,
  };

  constexpr Mixture::Specie Acetate = {
    "acetate",
    60.,
    std::nullopt,
  };

  constexpr Mixture::Specie CarbonDioxide = {
    "co2",
    44.,
    8.3e-1,
  };

  template <class R, class T>
  bool
  contains(std::initializer_list<R> r, const T& v)
  {
    return std::ranges::find(r, v) != std::ranges::end(r);
  }

  std::string
  to_lowercase(std::string s)
  {
    for (auto& c : s)
    {
      c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }

    return s;
  }

}

namespace Mixture
{

  std::optional<Specie>
  query_species(std::string_view name)
  {
    std::string n = to_lowercase(std::string(name));
    if (contains({ "g", "s", "gl", "glucose" }, n))
    {
      // overwrite name by user provided no avoid problem
      auto s = Glucose;
      s.name = name;
      return s;
    }

    if (contains({ "o2", "oxygen" }, n))
    {
      // overwrite name by user provided no avoid problem
      auto s = Oxygen;
      s.name = name;
      return s;
    }

    if (contains({ "ac", "acetate" }, n))
    {
      // overwrite name by user provided no avoid problem
      auto s = Acetate;
      s.name = name;
      return s;
    }

    if (contains({ "co2" }, n))
    {
      // overwrite name by user provided no avoid problem
      auto s = CarbonDioxide;
      s.name = name;
      return s;
    }

    return std::nullopt;
  }
} // namespace Mixture
