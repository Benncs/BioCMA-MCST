#include <Kokkos_Assert.hpp>
#include <cmath>
#include <iostream>
#include <optional>
#include <simulation/feed_descriptor.hpp>
#include <stdexcept>
#include <variant>
#define CHECK_TYPE_VARIANT(__variant_arg__, __ref__type)                       \
  std::is_same_v<std::decay_t<decltype(__variant_arg__)>, __ref__type>

namespace
{
  struct updateDispacther
  {
    Simulation::Feed::FeedDescriptor& self;
    double t;
    double d_t;

    inline void
    operator()(Simulation::Feed::Exponential& x) const noexcept
    {
      self.flow = x.f0 + std::exp(x.alpha * t);
    }

    inline void
    operator()(Simulation::Feed::Linear& x) const noexcept
    {
      self.flow = x.f0 + t * x.df;
    }

    template <typename T>
    inline void
    operator()(T& x) const noexcept
    {
      (void)x;
    }
  };

  template <typename variant_type>
  Simulation::Feed::FeedDescriptor
  gen_factory(variant_type variant,
              double flow,
              double concentration,
              std::size_t species_index,
              std::size_t input_position,
              std::optional<std::size_t> _ouput_position,
              bool set_output) noexcept
  {
    if (set_output && (!_ouput_position))
    {
      _ouput_position = input_position;
    }

    const auto value
        = Simulation::Feed::FeedValue{ concentration, species_index };
    const bool use_relative_time = true;
    return { flow,    { value },        input_position, _ouput_position,
             variant, use_relative_time };
  }

  template <typename variant> struct typeDispatcher
  {
  };

  template <> struct typeDispatcher<Simulation::Feed::Constant>
  {
    static constexpr FeedType type = FeedType::Constant;
  };
  template <> struct typeDispatcher<Simulation::Feed::Linear>
  {
    static constexpr FeedType type = FeedType::Linear;
  };
  template <> struct typeDispatcher<Simulation::Feed::Exponential>
  {
    static constexpr FeedType type = FeedType::Exponential;
  };
} // namespace

namespace Simulation::Feed
{
  FeedType
  get_type(const FeedTypeVariant& v)
  {
    return std::visit(
        [](auto&& arg) -> FeedType
        { return typeDispatcher<std::decay_t<decltype(arg)>>::type; },
        v);
  }

  void
  FeedDescriptor::update(double t, double d_t) noexcept
  {
    std::visit(updateDispacther{ *this, t, d_t }, extra);
  }

  FeedDescriptor
  FeedFactory::constant(double flow,
                        double concentration,
                        std::size_t species_index,
                        std::size_t input_position,
                        std::optional<std::size_t> _ouput_position,
                        bool set_output) noexcept
  {

    return gen_factory(Constant{},
                       flow,
                       concentration,
                       species_index,
                       input_position,
                       _ouput_position,
                       set_output);
  }

  FeedDescriptor
  FeedFactory::linear(double flow,
                      double df,
                      double concentration,
                      std::size_t species_index,
                      std::size_t input_position,
                      std::optional<std::size_t> _ouput_position,
                      bool set_output) noexcept
  {
    return gen_factory(Linear{ .f0 = flow, .df = df },
                       flow,
                       concentration,
                       species_index,
                       input_position,
                       _ouput_position,
                       set_output);
  }

  void
  SimulationFeed::add_feed(FeedDescriptor&& fd, Phase phase)
  {
    auto& vec = phase == Phase::Liquid ? liquid : gas;
    if (!vec)
    {
      vec = std::vector<FeedDescriptor>();
    }
    KOKKOS_ASSERT(fd.flow >= 0.);

    vec->emplace_back(
        move_allow_trivial(fd)); // Use move_allow_trivial in case
                                 // FeedDescriptor become non trivial
  }

  void
  SimulationFeed::add_liquid(FeedDescriptor&& fd)
  {
    add_feed(move_allow_trivial(fd), Phase::Liquid);
  }

  void
  SimulationFeed::add_gas(FeedDescriptor&& fd)
  {
    add_feed(move_allow_trivial(fd), Phase::Gas);
  }

  std::size_t
  SimulationFeed::n_liquid_flow() const noexcept
  {
    return (liquid) ? liquid->size() : 0;
  }
  std::size_t
  SimulationFeed::n_gas_flow() const noexcept
  {
    return (gas) ? liquid->size() : 0;
  }

  SimulationFeed
  SimulationFeed::empty() noexcept
  {
    return { .liquid = std::nullopt, .gas = std::nullopt };
  }

} // namespace Simulation::Feed
