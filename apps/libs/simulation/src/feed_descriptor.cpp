#include <iostream>
#include <simulation/feed_descriptor.hpp>
#include <stdexcept>
#include <utility>
#include <variant>

#define CHECK_TYPE_VARIANT(__variant_arg__, __ref__type)                                           \
  std::is_same_v<std::decay_t<decltype(__variant_arg__)>, __ref__type>

// namespace {
//   struct FeedFunctor
//   {
//     void operator()(Simulation::Feed::DelayedConstant& c)
//     {

//     }
//   }

// } //namespace

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

    switch (this->type)
    {

    case FeedType::Constant:
    {
      return;
    }
    case FeedType::DelayedConstant:
    {
      auto& inner = std::get<DelayedConstant>(this->props);
      if (t >= inner.t_init && t <= inner.t_end)
      {
        flow_value = inner.stored_value;
      }
      else
      {
        flow_value = 0;
      }
      break;
    }
    case FeedType::Step:
    {
      __builtin_unreachable();
      break;
    }
    case FeedType::Pulse:
    {
      __builtin_unreachable();
      break;
    }
    case FeedType::Custom:
    {
      __builtin_unreachable();
    }
    break;
    }
  }

  FeedDescritor::FeedDescritor(double _f,
                               feed_value_t&& _target,
                               feed_position_t&& _position,
                               feed_species_t _species,
                               FeedTypeVariant _props,
                               bool _set_exit)
      : flow_value(_f), position(std::move(_position)), species(std::move(_species)), props(_props),
        set_exit(_set_exit), n_v(_target.size()), type(get_type(props)), target(std::move(_target))
  {

    value = target;

    // TODO check
    if (value.size() != position.size())
    {
      std::cout << value.size() << " " << position.size() << std::endl;
      throw std::invalid_argument("Bad feed size");
    }
  }

  FeedDescritor FeedFactory::constant(double _f,
                                      feed_value_t&& _target,
                                      feed_position_t&& _position,
                                      feed_species_t _species,
                                      bool set_exit)
  {
    return {
        _f, std::move(_target), std::move(_position), std::move(_species), Constant{}, set_exit};
  }

  FeedDescritor delayedconstant(double _f,
                                feed_value_t&& _target,
                                feed_position_t&& _position,
                                feed_species_t _species,
                                double t_init,
                                double t_end,
                                bool set_exit)
  {
    return {0.,
            std::move(_target),
            std::move(_position),
            std::move(_species),
            DelayedConstant{t_init, t_end, _f},
            set_exit};
  }

} // namespace Simulation::Feed