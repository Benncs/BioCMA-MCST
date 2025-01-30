#include <iostream>
#include <simulation/feed_descriptor.hpp>
#include <stdexcept>
#include <utility>
#include <variant>

#define CHECK_TYPE_VARIANT(__variant_arg__, __ref__type)                                           \
  std::is_same_v<std::decay_t<decltype(__variant_arg__)>, __ref__type>

namespace Simulation::Feed
{
  FeedType get_type(const FeedTypeVariant& v)
  {
    return std::visit(
        [](auto&& arg) -> FeedType
        {
          if constexpr (CHECK_TYPE_VARIANT(arg, Constant))
          {
            return FeedType::Constant;
          }
          else if constexpr (CHECK_TYPE_VARIANT(arg, Step))
          {
            return FeedType::Step;
          }
          else if constexpr (CHECK_TYPE_VARIANT(arg, Pulse))
          {
            return FeedType::Pulse;
          }
          else if constexpr (CHECK_TYPE_VARIANT(arg, Custom))
          {
            return FeedType::Custom;
          }
          else
          {
            throw std::runtime_error("Unsupported type");
          }
        },
        v);
  }

  void FeedDescritor::update(double t, double d_t) noexcept
  {
  }

  FeedDescritor::FeedDescritor(double _f,
                               feed_value_t&& _target,
                               feed_position_t&& _position,
                               feed_species_t _species,
                               FeedTypeVariant _props)
      : flow_value(_f), position(std::move(_position)), species(std::move(_species)), props(_props),
        n_v(_target.size()), type(get_type(props)), target(std::move(_target))
  {

    value = target;

    // TODO check
    if (value.size() != position.size())
    {
      std::cout << value.size() << " " << position.size() << std::endl;
      throw std::invalid_argument("Bas feed size");
    }
  }

  FeedDescritor FeedFactory::constant(double _f,
                                      feed_value_t&& _target,
                                      feed_position_t&& _position,
                                      feed_species_t _species)
  {
    return {_f, std::move(_target), std::move(_position), std::move(_species), Constant{}};
  }

} // namespace Simulation::Feed