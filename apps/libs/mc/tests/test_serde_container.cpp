#include "Kokkos_Macros.hpp"
#include <mc/particles_container.hpp>
#include "serdes_class.hpp"
#include <Kokkos_Core.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/types/array.hpp> //MC::events use array internally
#include <cereal/types/memory.hpp>
#include <cereal/types/optional.hpp>
#include <cereal/types/string.hpp> //MC::List use vector internally
#include <cereal/types/tuple.hpp>
#include <cereal/types/variant.hpp> //MC::Unit use variant internally
#include <cereal/types/vector.hpp>  //MC::List use vector internally
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mc/prng/prng.hpp>
#include <sstream>

static constexpr std::size_t p_size = 100;

void write_to_file(const std::ostringstream& oss, const std::string& filename)
{
  auto tmp_dir = std::filesystem::temp_directory_path();

  auto file_path = tmp_dir / filename;
  std::ofstream file(file_path);
  if (!file)
  {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  file << oss.str();
  file.close();
}


int main()
{
  Kokkos::ScopeGuard _guad;
  const std::size_t np = 100;
  MC::ParticlesContainer<SerdeModel> container(np);
  MC::KPRNG::pool_type rng;
  Kokkos::parallel_for(np,KOKKOS_LAMBDA(const int i){
    SerdeModel::init(rng,i,container.model);
  });

  std::ostringstream buff(std::ios::binary);
  cereal::BinaryOutputArchive oarchive(buff);
  oarchive(container);
  std::istringstream iss(buff.str(), std::ios::binary);
  cereal::BinaryInputArchive iarchive(iss);
  MC::ParticlesContainer<SerdeModel> container2;
  iarchive(container2);


  assert(container.n_particles()==container2.n_particles());

   Kokkos::parallel_for(np,KOKKOS_LAMBDA(const int i){
    for(std::size_t j=0;j<SerdeModel::n_var;++j)
    {
      KOKKOS_ASSERT(container.model(i,j)==container2.model(i,j));
    }
    
  });

  
  
}