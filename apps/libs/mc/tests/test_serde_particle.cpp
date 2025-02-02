#include "serdes_class.hpp"

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
  assert(loaded_part.data.property1 == -1); // It should not have changed
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
  assert(loaded_part.data.property1 == part.data.property1); // We serialized and desearilzed  model
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