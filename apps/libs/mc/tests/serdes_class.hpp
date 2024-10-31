#ifndef __UTILS_SERDES_HPP__
#define __UTILS_SERDES_HPP__

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
  int property1 = -1;

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
    property1 = _rng.uniform_u(0, 10000);
  }

  template <class Archive> void serialize(Archive &ar)
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
  t.serialize(ar);
  return buff;
}

template <typename T, class StreamType> void wrap_de(T &t, StreamType &buff)
{
  std::istringstream iss(buff.str(), std::ios::binary);
  cereal::BinaryInputArchive ar(iss);
  t.serialize(ar);
}

#endif 