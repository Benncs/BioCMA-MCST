#include "mc/particles/data_holder.hpp"
#include "mc/particles/particle_model.hpp"
#include <cassert>
#include <mc/particles/mcparticles.hpp>

#include <cereal/archives/binary.hpp>
#include <sstream>

// NOLINTBEGIN
class NonSerdeModel
{
public:
  int property1 = 0;

  void init(MC::ParticleDataHolder &p, MC::KPRNG _rng)
  {
    p.status = MC::CellStatus::IDLE;
  }

  void update(double d_t,
              MC::ParticleDataHolder &p,
              const LocalConcentrationView &concentration,
              MC::KPRNG _rng)
  {
  }

  NonSerdeModel division(MC::ParticleDataHolder & /*p*/,MC::KPRNG)
  {
    return {};
  }

  void contribution(MC::ParticleDataHolder &p, ContributionView contrib)
  {
  }

  inline model_properties_detail_t get_properties()
  {
    return {};
  }
  double mass()const{return 0;}
};

class SerdeModel
{
public:
  int property1 = 0;

  void init(MC::ParticleDataHolder &p, MC::KPRNG _rng)
  {
    p.status = MC::CellStatus::IDLE;
  }

  template <class Archive> void serde(Archive &ar)
  {
    ar(property1);
  }

  void update(double d_t,
              MC::ParticleDataHolder &p,
              const LocalConcentrationView &concentration,
              MC::KPRNG _rng)
  {
  }

  SerdeModel division(MC::ParticleDataHolder & /*p*/,MC::KPRNG)
  {
    return {};
  }

  void contribution(MC::ParticleDataHolder &p, ContributionView contrib)
  {
  }

  inline model_properties_detail_t get_properties()
  {
    return {};
  }
   double mass()const{return 0;}
};
// NOLINTEND
template <typename T> std::ostringstream wrap_ser(T &t)
{
  std::ostringstream buff(std::ios::binary);
  cereal::BinaryOutputArchive ar(buff);
  t.serde(ar);
  return buff;
}

template <typename T, class StreamType> void wrap_de(T &t, StreamType &buff)
{
  std::istringstream iss(buff.str(), std::ios::binary);
  cereal::BinaryInputArchive ar(iss);
  t.serde(ar);
}

void test_data_holder()
{
  constexpr double init_weight = 0.81;
  constexpr size_t init_container = 50;
  MC::ParticleDataHolder raw_data(init_weight);
  raw_data.current_container = init_container;

  auto buff = wrap_ser(raw_data);

  MC::ParticleDataHolder loaded_data(0);

  wrap_de(loaded_data, buff);

  assert(loaded_data.weight == init_weight);
  assert(loaded_data.current_container == init_container);
  assert(loaded_data.current_domain == raw_data.current_domain);
  assert(loaded_data.id == raw_data.id);
  assert(loaded_data.random_seed == raw_data.random_seed);
}

void test_particle_with_non_model()
{
  // NonSerdeModel doesn't expose serde function so this should
  // compile but so serializion of model should be done
  MC::Particle<NonSerdeModel> part;
  part.data.property1 = 55;
  part.properties.id = 80;
  auto buff = wrap_ser(part);

  MC::Particle<NonSerdeModel> loaded_part;

  wrap_de(loaded_part, buff);

  assert(part.properties.id == loaded_part.properties.id);
  assert(loaded_part.data.property1 == 0); // It should not have changed
}

void test_particle_with_model()
{
  // NonSerdeModel expose serde function so this should
  // compile and serializion should happen
  MC::Particle<SerdeModel> part;
  part.data.property1 = 55;
  part.properties.id = 80;
  auto buff = wrap_ser(part);

  MC::Particle<SerdeModel> loaded_part;

  wrap_de(loaded_part, buff);

  assert(part.properties.id == loaded_part.properties.id);
  assert(loaded_part.data.property1 ==
         part.data.property1); // We serialized and desearilzed  model
}

int main()
{
  test_data_holder();
  test_particle_with_non_model();
  test_particle_with_model();

  //   MC::Particle<DefaultModel> part;

  //   part.properties.current_container = 50;

  //   std::ostringstream buff;

  //   {
  //     cereal::BinaryOutputArchive ar(buff);

  //     part.serde(ar);
  //   }

  //   std::string binaryData = buff.str();

  //   std::cout << "Dumb value: ";
  //   for (unsigned char c : binaryData)
  //   {
  //     std::cout << std::hex << std::setw(2) << std::setfill('0') << (int)c <<
  //     " ";
  //   }
  //   std::cout << std::endl;
}