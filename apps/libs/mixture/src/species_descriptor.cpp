#include <algorithm>
#include <iomanip>
#include <mixture/species_descriptor.hpp>
#include <optional>
#include <ostream>
namespace Mixture
{

  std::ostream&
  operator<<(std::ostream& os, const Specie& s)
  {
    os << "{ "
       << "name=" << std::setw(10) << std::left << s.name << " "
       << "mw=" << std::setw(8) << (s.molar_weight ? *s.molar_weight : 0) << " "
       << "henry=" << std::setw(8) << (s.henry ? *s.henry : 0) << " }";

    return os;
  }

  std::ostream&
  operator<<(std::ostream& os, const SpecieTable& table)
  {
    os << "Species (" << table.m_table.size() << "):\n";
    for (std::size_t i = 0; i < table.m_table.size(); ++i)
    {
      os << "  [" << i << "] " << table.m_table[i] << '\n';
    }
    return os;
  }

  Specie
  new_specie(std::string_view _name,
             double molar_weight,
             std::optional<double> henry)
  {
    return {
      .name = std::string(_name),
      .molar_weight = molar_weight,
      .henry = henry,
    };
  }

  Specie
  new_specie(std::string_view _name)
  {
    return {
      .name = std::string(_name),
      .molar_weight = std::nullopt,
      .henry = std::nullopt,
    };
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

  SpecieTable::SpecieTable(std::size_t n)
  {
    for (std::size_t i = 0LU; i < n; ++i)
    {
      add(new_specie(std::to_string(i)));
    }
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

} // namespace Mixture
