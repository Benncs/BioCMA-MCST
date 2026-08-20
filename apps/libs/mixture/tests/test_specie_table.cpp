#include <cassert>
#include <mixture/species_descriptor.hpp>
#include <optional>
#include <vector>
int
main()
{

  const auto species1
      = Mixture::Specie{ .name = "sss", .molar_weight = 1., .henry = 0. };
  const auto specie2 = Mixture::Specie{ .name = "sssaava",
                                        .molar_weight = 2.,
                                        .henry = std::nullopt };

  const auto specie3
      = Mixture::Specie{ .name = "ssssavaa", .molar_weight = 3., .henry = 45. };

  auto table = Mixture::SpecieTable({ species1, specie2, specie3 });

  assert(table.n_species() == 3);

  auto vh = table.henry();

  assert(vh[0] == 0.);
  assert(vh[1] == 0.);
  assert(vh[2] == 45.);

  table.add(*Mixture::query_species("g"));
  assert(table.n_species() == 4);
  assert(table.henry().size() == 4);

  return 0;
}
