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
  Exponential,
  Linear,
  // DelayedConstant,
  // Step,
  // Pulse,
  // Custom
};

template <typename T>
constexpr decltype(auto)
move_allow_trivial(T&& t) noexcept
{
  return std::move(t); // NOLINT
}

namespace Simulation::Feed
{

  struct Constant
  {
  };

  struct Exponential
  {
    double f0;
    double alpha;
  };

  struct Linear
  {
    double f0;
    double df;
  };

  // struct DelayedConstant
  // {
  //   double t_init;
  //   double t_end;
  //   double stored_value;
  // };

  // struct Step
  // {
  //   double t_init = 0.;
  // };

  // struct Pulse
  // {
  //   double t_init;
  //   double t_end;
  //   double frequency;
  //   double stored_value;
  // };

  // struct Custom
  // {
  // };

  // using FeedTypeVariant
  //     = std::variant<Constant, Step, Pulse, Custom, DelayedConstant>;

  using FeedTypeVariant = std::variant<Constant, Linear, Exponential>;

  FeedType get_type(const FeedTypeVariant& v);

  struct FeedValue
  {
    double concentration;
    std::size_t species_index;
  };

  struct FeedDescriptor
  {
    double flow{};
    std::vector<FeedValue> values;
    std::size_t input_position{};
    std::optional<std::size_t> output_position;
    FeedTypeVariant extra;
    bool use_relative_time = true;
    void update(double t, double d_t) noexcept;
  };

  struct FeedFactory
  {
    static FeedDescriptor constant(double flow,
                                   double concentration,
                                   std::size_t species_index,
                                   std::size_t input_position,
                                   std::optional<std::size_t> _ouput_position
                                   = std::nullopt,
                                   bool set_output = true) noexcept;
  };

  class SimulationFeed
  {
  public:
    std::optional<std::vector<FeedDescriptor>> liquid;
    std::optional<std::vector<FeedDescriptor>> gas;

    void add_liquid(FeedDescriptor&& fd);

    void add_gas(FeedDescriptor&& fd);

    void add_feed(FeedDescriptor&& fd, Phase phase);

    [[nodiscard]] std::size_t n_liquid_flow() const noexcept;
    [[nodiscard]] std::size_t n_gas_flow() const noexcept;

    auto
    liquid_feeds()
    {
      if (liquid)
      {
        return std::ranges::subrange(liquid->begin(), liquid->end());
      }
      return std::ranges::subrange(std::vector<FeedDescriptor>::iterator(),
                                   std::vector<FeedDescriptor>::iterator());
    }

    auto
    gas_feeds()
    {
      if (gas)
      {
        return std::ranges::subrange(gas->begin(), gas->end());
      }
      return std::ranges::subrange(std::vector<FeedDescriptor>::iterator(),
                                   std::vector<FeedDescriptor>::iterator());
    }

    static SimulationFeed empty() noexcept;
  };

} // namespace Simulation::Feed

#endif
