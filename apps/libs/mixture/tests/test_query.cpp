#include <cassert>
#include <mixture/species_descriptor.hpp>
#include <optional>
int
main()
{
  // "g", "s", "gl", "glucose"
  assert(Mixture::query_species("g").has_value());
  assert(Mixture::query_species("G").has_value());
  assert(Mixture::query_species("s").has_value());
  assert(Mixture::query_species("gl").has_value());
  assert(Mixture::query_species("glucose").has_value());

  auto g = Mixture::query_species("G");
  assert(g->name == "G");
  assert(g->henry == std::nullopt);
  assert(g->molar_weight == 180.);

  // "o2", "oxygen"
  assert(Mixture::query_species("o2").has_value());
  assert(Mixture::query_species("oxygen").has_value());

  auto o = Mixture::query_species("o2");
  assert(o->name == "o2");
  assert(o->molar_weight == 32.);

  // not found -> none
  assert(!Mixture::query_species("ssssaaaavvaaa").has_value());

  return 0;
}
