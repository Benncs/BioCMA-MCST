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

    return std::nullopt;
  }

  [[nodiscard]] std::size_t
  SpecieTable::n_species() const
  {
    return m_table.size();
  }

  SpecieTable::SpecieTable(std::initializer_list<Specie> values)
      : m_table(values)
  {
  }

  [[nodiscard]] std::optional<size_t>
  SpecieTable::index_of(std::string_view name) const
  {
    auto it = std::ranges::find_if(
        m_table, [name](const auto& n) { return n.name == name; });

    if (it == m_table.end())
    {
      return std::nullopt;
    }

    return it - m_table.begin();
  }

  [[nodiscard]] bool
  SpecieTable::has(std::string_view name) const
  {
    return index_of(name).has_value();
  }

  std::optional<std::reference_wrapper<const Specie>>
  SpecieTable::add(Specie&& s)
  {
    if (!has(s.name))
    {
      m_table.emplace_back(std::move(s));
      return m_table.back();
    }

    return std::nullopt;
  }

  std::optional<std::reference_wrapper<const Specie>>
  SpecieTable::get(std::string_view name) const
  {
    auto index_existing_s = index_of(name);

    if (index_existing_s.has_value())
    {
      return m_table[*index_existing_s];
    }

    return std::nullopt;
  }

}
