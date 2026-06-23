#ifndef __MIXTURE_SPECIES_DESCRIPTOR_HPP__
#define __MIXTURE_SPECIES_DESCRIPTOR_HPP__
#include <initializer_list>
#include <optional>
#include <ostream>
#include <ranges>
#include <string>
#include <vector>
namespace Mixture
{

  struct Specie
  {
    std::string name;
    std::optional<double> molar_weight;
    // dimensionless
    std::optional<double> henry;

    template <class Archive>
    void
    serialize(Archive& archive)
    {
      archive(name, molar_weight, henry);
    }
  };

  Specie new_specie(std::string_view _name,
                    double molar_weight,
                    std::optional<double> henry);

  Specie new_specie(std::string_view _name);

  std::optional<Specie> query_species(std::string_view);

  class SpecieTable
  {
  private:
    std::vector<Specie> m_table;

  public:
    friend std::ostream& operator<<(std::ostream&, const SpecieTable&);
    [[nodiscard]] std::size_t n_species() const;
    SpecieTable(std::initializer_list<Specie> values);

    template <class Archive>
    void
    serialize(Archive& archive)
    {
      archive(m_table);
    }

    explicit SpecieTable(std::size_t n);

    SpecieTable() = default;
    std::optional<std::reference_wrapper<const Specie>> add(Specie&& s);
    [[nodiscard]] std::optional<std::reference_wrapper<const Specie>>
    get(std::string_view name) const;

    [[nodiscard]] std::optional<size_t> index_of(std::string_view name) const;

    [[nodiscard]] bool has(std::string_view name) const;

    [[nodiscard]] auto
    henry() const
    {
      return m_table
             | std::views::transform([](const Specie& e)
                                     { return e.henry.value_or(0.0); });
    }

    [[nodiscard]] auto
    species_list() const
    {
      return m_table
             | std::views::transform([](const Specie& e) { return e.name; });
    }
  };
  std::ostream& operator<<(std::ostream& os, const Specie& s);

  std::ostream& operator<<(std::ostream& os, const SpecieTable& table);

  struct EnvironementProperties
  {
    std::vector<Specie> species;
    double bubble_diameter;
  };

} // namespace Mixture

#endif //__SIMULATION_SPECIES_DESCRIPTOR_HPP__
