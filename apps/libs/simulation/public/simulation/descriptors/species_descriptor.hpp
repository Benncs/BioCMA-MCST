#ifndef __SIMULATION_SPECIES_DESCRIPTOR_HPP__
#define __SIMULATION_SPECIES_DESCRIPTOR_HPP__
#include <optional>
#include <string>
#include <vector> 
namespace Simulation
{

  struct Specie
  {
    std::string name; 
    std::optional<double> molar_weight; 
    std::optional<double> henry;
  };

  using SpecieTable = std::vector<Specie>;

  struct EnvironementProperties
  {
    std::vector<Specie> species;
    double bubble_diameter;
  };

  

} // namespace Simulation

#endif //__SIMULATION_SPECIES_DESCRIPTOR_HPP__