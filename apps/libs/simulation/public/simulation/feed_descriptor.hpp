#ifndef __SIMULATION_FEED_DESCRIPTOR_HPP__
#define __SIMULATION_FEED_DESCRIPTOR_HPP__

#include <cstddef>
#include <cstdint>
#include <optional>
#include <ranges>
#include <variant>
#include <vector>

enum class Phase
{
  Liquid,
  Gas
};

enum class FeedType : std::uint8_t
{
  Constant,
  DelayedConstant,
  Step,
  Pulse,
  Custom
};

template<typename T>
constexpr decltype(auto) move_allow_trivial(T&& t) noexcept
{
    return std::move(t); // NOLINT
}

namespace Simulation::Feed
{
  using feed_value_t = std::vector<double>;
  using feed_position_t = std::vector<std::size_t>;
  using feed_species_t = std::vector<std::size_t>;

  struct Constant
  {
  };

  struct DelayedConstant
  {
    double t_init;
    double t_end;
    double stored_value;
  };

  struct Step
  {
    double t_init = 0.;
  };

  struct Pulse
  {
    double t_init;
    double t_end;
    double frequency;
    double stored_value;
  };

  struct Custom
  {
  };

  using FeedTypeVariant = std::variant<Constant, Step, Pulse, Custom, DelayedConstant>;

  FeedType get_type(const FeedTypeVariant& v);

  struct FeedDescriptor
  {
    double flow{};
    double concentration{};
    std::size_t species_index;
    std::size_t input_position{};
    std::optional<std::size_t> output_position;
    FeedTypeVariant extra;

    void update(double t, double d_t) noexcept;
  };

  struct FeedFactory
  {
    static FeedDescriptor constant(double flow,
                                  double concentration,
                                  std::size_t species_index,
                                  std::size_t input_position,
                                  std::optional<std::size_t > _ouput_position = std::nullopt,
                                  bool set_output = true);

    // static FeedDescriptor delayedconstant(double _f,
    //                                      feed_value_t&& _target,
    //                                      feed_position_t&& _position,
    //                                      feed_species_t _species,
    //                                      double t_init,
    //                                      double t_end,
    //                                      bool set_output = true);

    // static FeedDescriptor pulse(double _f,
    //                            feed_value_t&& _target,
    //                            feed_position_t&& _position,
    //                            feed_species_t _species,
    //                            double t_init,
    //                            double t_end,
    //                            double frequency,
    //                            bool set_output = true);
  };

  class SimulationFeed
  {
  public:
    std::optional<std::vector<FeedDescriptor>> liquid;
    std::optional<std::vector<FeedDescriptor>> gas;

    void add_liquid(FeedDescriptor&& fd);

    void add_gas(FeedDescriptor&& fd);

    void add_feed(FeedDescriptor&& fd, Phase phase);

    std::size_t n_liquid_flow()const;
    std::size_t n_gas_flow()const;

    auto liquid_feeds()
    {
      if (liquid)
      {
        return std::ranges::subrange(liquid->begin(), liquid->end());
      }
      return std::ranges::subrange(std::vector<FeedDescriptor>::iterator(),
                                   std::vector<FeedDescriptor>::iterator());
    }

    auto gas_feeds()
    {
      if (gas)
      {
        return std::ranges::subrange(gas->begin(), gas->end());
      }
      return std::ranges::subrange(std::vector<FeedDescriptor>::iterator(),
                                   std::vector<FeedDescriptor>::iterator());
    }
  };

} // namespace Simulation::Feed

#endif