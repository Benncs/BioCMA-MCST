#include <cmath>
#include <simulation/feed_descriptor.hpp>
#include <stdexcept>
#include <variant>

#define CHECK_TYPE_VARIANT(__variant_arg__, __ref__type)                       \
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

  void FeedDescriptor::update(double t, double d_t) noexcept
  {

#pragma message("Update feed not implemented")
    (void)t;
    (void)d_t;
  }

  // FeedDescritor::FeedDescritor(double _f,
  //                              feed_value_t&& _target,
  //                              feed_position_t&& _position,
  //                              feed_species_t _species,
  //                              FeedTypeVariant _props,
  //                              bool _set_exit)
  //     : flow_value(_f), position(std::move(_position)),
  //     species(std::move(_species)), props(_props),
  //       set_exit(_set_exit), n_v(_target.size()), type(get_type(props)),
  //       target(std::move(_target))
  // {

  //   value = target;

  //   // TODO check
  //   if (value.size() != position.size())
  //   {
  //     std::cout << value.size() << " " << position.size() << std::endl;
  //     throw std::invalid_argument(
  //         "Feed descriptor: Number of position should be the same as value");
  //   }

  //   if (_set_exit)
  //   {
  //     ouput_position = position;
  //   }
  // }

  // FeedDescritor::FeedDescritor(double _f,
  //                              feed_value_t&& _target,
  //                              feed_position_t&& _position,
  //                              std::optional<feed_position_t>&&
  //                              _ouput_position, feed_species_t _species,
  //                              FeedTypeVariant _props,
  //                              bool _set_exit)
  //     : flow_value(_f), position(std::move(_position)),
  //     species(std::move(_species)), props(_props),
  //       set_exit(_set_exit), n_v(_target.size()), type(get_type(props)),
  //       target(std::move(_target))
  // {

  //   value = target;

  //   // TODO check
  //   if (value.size() != position.size())
  //   {
  //     std::cout << value.size() << " " << position.size() << std::endl;
  //     throw std::invalid_argument(
  //         "Feed descriptor: Number of position should be the same as value");
  //   }

  //   if (_set_exit)
  //   {
  //     if (_ouput_position)
  //     {

  //       ouput_position = std::move(*_ouput_position);
  //     }
  //     else
  //     {
  //       ouput_position = _position;
  //     }

  //     if (ouput_position.size() != position.size())
  //     {
  //       throw std::invalid_argument(
  //           "Feed descriptor: Number of position should be the same as
  //           ouput_position");
  //     }
  //   }
  // }

  FeedDescriptor
  FeedFactory::constant(double flow,
                        double concentration,
                        std::size_t species_index,
                        std::size_t input_position,
                        std::optional<std::size_t> _ouput_position,
                        bool set_output)
  {

    if (set_output && (!_ouput_position))
    {
      _ouput_position = input_position;
    }

    return {flow,
            concentration,
            species_index,
            input_position,
            _ouput_position,
            Constant{}};
  }

  // FeedDescriptor delayedconstant(double _f,
  //                                feed_value_t&& _target,
  //                                feed_position_t&& _position,
  //                                feed_species_t _species,
  //                                double t_init,
  //                                double t_end,
  //                                bool set_output)
  // {
  //   // return {0.,
  //   //         std::move(_target),
  //   //         std::move(_position),
  //   //         std::move(_species),
  //   //         DelayedConstant{t_init, t_end, _f},
  //   //         set_output};
  // }

  // FeedDescriptor FeedFactory::pulse(double _f,
  //                                   feed_value_t&& _target,
  //                                   feed_position_t&& _position,
  //                                   feed_species_t _species,
  //                                   double t_init,
  //                                   double t_end,
  //                                   double frequency,
  //                                   bool set_output)
  // {
  //   // return {_f,
  //   //         std::move(_target),
  //   //         std::move(_position),
  //   //         std::move(_species),
  //   //         Pulse{t_init, t_end, frequency, _f},
  //   //         set_output};
  // }

  void SimulationFeed::add_feed(FeedDescriptor&& fd, Phase phase)
  {
    auto& vec = phase == Phase::Liquid ? liquid : gas;
    if (!vec)
    {
      vec = std::vector<FeedDescriptor>();
    }
    vec->emplace_back(
        move_allow_trivial(fd)); // Use move_allow_trivial in case
                                 // FeedDescriptor become non trivial
  }

  void SimulationFeed::add_liquid(FeedDescriptor&& fd)
  {
    add_feed(move_allow_trivial(fd), Phase::Liquid);
  }

  void SimulationFeed::add_gas(FeedDescriptor&& fd)
  {
    add_feed(move_allow_trivial(fd), Phase::Gas);
  }

  std::size_t SimulationFeed::n_liquid_flow() const
  {
    return (liquid) ? liquid->size() : 0;
  }
  std::size_t SimulationFeed::n_gas_flow() const
  {
    return (gas) ? liquid->size() : 0;
  }

} // namespace Simulation::Feed