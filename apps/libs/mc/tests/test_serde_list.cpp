#include "cereal/details/helpers.hpp"
#include "common/kokkos_vector.hpp"
#include "mc/particles/mcparticles.hpp"
#include "mc/particles/particle_model.hpp"
#include "mc/particles/particles_container.hpp"
#include "mc/prng/prng.hpp"
#include "serdes_class.hpp"

#include "mc/particles/particle_list.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <cereal/types/vector.hpp>
#include <cstddef>
#include <fstream>
#include <sstream>


static constexpr std::size_t p_size = 100;


void write_to_file(const std::ostringstream &oss, const std::string &filename)
{
  std::ofstream file(filename);
  if (!file)
  {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }
  file << oss.str();
  file.close();
}

void test_vec_serde()
{
  std::vector<MC::Particle<SerdeModel>> de(p_size);

  de[0].properties.id = 1;
  de[80].data.property1 = 10;

  std::ostringstream buf(std::ios::binary);
  {
    cereal::BinaryOutputArchive ar(buf);
    ar(de);
  }

  std::vector<MC::Particle<SerdeModel>> se;
  std::istringstream iss(buf.str(), std::ios::binary);
  {
    cereal::BinaryInputArchive ar(iss);
    ar(se);
  }
  assert(se.size() == de.size());
  assert(se[0].properties.id == de[0].properties.id);
  assert(de[80].data.property1 == se[80].data.property1);
}

void test_vec_non_serde()
{
  std::vector<MC::Particle<NonSerdeModel>> de(p_size);

  de[0].properties.id = 1;
  de[80].data.property1 = 10;

  std::ostringstream buf(std::ios::binary);
  {
    cereal::BinaryOutputArchive ar(buf);
    ar(de);
  }

  std::vector<MC::Particle<NonSerdeModel>> se;
  std::istringstream iss(buf.str(), std::ios::binary);
  {
    cereal::BinaryInputArchive ar(iss);
    ar(se);
  }
  assert(se.size() == de.size());
  assert(se[0].properties.id == de[0].properties.id);
  assert(se[80].data.property1 == -1);
}

void test_list_serde()
{
  auto list = MC::ParticleList<ComputeSpace, SerdeModel>(p_size);

  MC::KPRNG rng;
  Kokkos::parallel_for(list.size(), KOKKOS_LAMBDA(auto &i) { list._owned_data(i).init(rng); });

  std::ostringstream buf(std::ios::binary);
  {
    cereal::BinaryOutputArchive ar(buf);
    ar(list);
  }
  //   auto buf = wrap_ser(list);
  write_to_file(buf, "tmp_serde");

  MC::ParticleList<ComputeSpace, SerdeModel> slist;
  std::istringstream iss(buf.str(), std::ios::binary);
  cereal::BinaryInputArchive ar(iss);
  ar(slist);

  assert(list.size() == slist.size());
  assert(list.capacity() == slist.capacity());
  assert(list._owned_data[9].data.property1 == slist._owned_data[9].data.property1);

}

void test_container_serde()
{
  auto container = MC::ParticlesContainer<SerdeModel>(p_size);

  std::ostringstream buf(std::ios::binary);
  {

    auto l = container.get_compute();

    Kokkos::parallel_for(p_size, KOKKOS_LAMBDA(auto i) { l._owned_data(i).data.property1 = i; });

    cereal::BinaryOutputArchive ar(buf);
    ar(container);
  }

//   write_to_file(buf, "tmp_serde");

  MC::ParticlesContainer<SerdeModel> secontainer;
  std::istringstream iss(buf.str(), std::ios::binary);
  cereal::BinaryInputArchive ar(iss);
  ar(secontainer);

  auto l = container.get_host();
  auto sl = secontainer.get_host();

  assert(l.size() == sl.size() && l.size()==p_size);
  assert(l.capacity() == sl.capacity());
  std::cerr<<l._owned_data[9].data.property1<<" "<<sl._owned_data[9].data.property1<<std::endl;
  assert(l._owned_data[9].data.property1 == sl._owned_data[9].data.property1 && l._owned_data[9].data.property1 == 9);
}

int main()
{
  test_vec_non_serde();
  
  test_vec_serde();

  Kokkos::initialize();
  test_list_serde();
  test_container_serde();
  Kokkos::finalize();
}