#ifndef __DIMENSIONS_HPP__
#define __DIMENSIONS_HPP__

#include <cstdlib>

namespace Simulation
{
  struct Dimensions
  {
    std::size_t n_species{};
    std::size_t n_compartment{};

    template <class Archive>
    void
    serialize(Archive& archive)
    {
      archive(n_species, n_compartment);
    }
  };
} // namespace Simulation

#endif
