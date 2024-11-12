#ifndef __SIMULATION_FEED_DESCRIPTOR_HPP__
#define __SIMULATION_FEED_DESCRIPTOR_HPP__

#include <optional>
#include <variant>
#include <vector>

enum class FeedType
{
  Constant,
  Step,
  Pulse,
  Custom
};

namespace Simulation::Feed
{
  using feed_value_t = std::vector<double>;
  using feed_position_t = std::vector<std::size_t>;
  using feed_species_t = std::vector<std::size_t>;
  struct Constant
  {
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
  };

  struct Custom
  {
  };

  using FeedTypeVariant = std::variant<Constant, Step, Pulse, Custom>;

  FeedType get_type(const FeedTypeVariant &v);

  class FeedDescritor
  {
  public:
    FeedDescritor() = default;
    FeedDescritor(double _f, feed_value_t &&_target, feed_position_t &&_position, feed_species_t _species, FeedTypeVariant _props);

    double flow_value{};
    feed_value_t value;
    feed_position_t position;
    feed_species_t species;
    FeedTypeVariant props;
    size_t n_v{};
    void update(double t, double d_t) noexcept;

  private:
    FeedType type{};
    feed_value_t target;
  };

  struct FeedFactory
  {
      static FeedDescritor constant(double _f, feed_value_t &&_target, feed_position_t &&_position, feed_species_t _species);
  };

  struct SimulationFeed
  {
    std::optional<std::vector<FeedDescritor>> liquid;
    std::optional<std::vector<FeedDescritor>> gas;
  };

  // struct Visitor
  // {
  //   void operator()(Uniform args);
  //   void operator()(Local args);
  //   void operator()(File filepath);
  //   void operator()(CustomScript path);
  // };
} // namespace Simulation::Feed

#endif