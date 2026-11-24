#include <Kokkos_Assert.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <optional>
#include <simulation/feed_descriptor.hpp>
#include <simulation/simulation_exception.hpp>
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
              bool set_output)
  {
    if (set_output && (!_ouput_position))
    {
      _ouput_position = input_position;
    }

    if (flow < 0.)
    {
      throw Simulation::FeedException(
          Simulation::FeedExceptionError::NegativeFlow);
    }
    if (concentration <= 0.)
    {
      throw Simulation::FeedException(
          Simulation::FeedExceptionError::NegativeConcentration);
    }

    const Simulation::Feed::FeedValue value{ concentration, species_index };
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
                        bool set_output)
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
  FeedFactory::sink(double flow, std::size_t _ouput_position)
  {
    if (flow < 0.)
    {
      throw Simulation::FeedException(
          Simulation::FeedExceptionError::NegativeFlow);
    }
    FeedDescriptor fd;
    fd.flow = flow;
    fd.extra = Constant{};
    fd.input_position = 0;
    fd.values = {};
    fd.output_position = _ouput_position;
    return fd;
  }

  FeedDescriptor
  FeedFactory::linear(double flow,
                      double df,
                      double concentration,
                      std::size_t species_index,
                      std::size_t input_position,
                      std::optional<std::size_t> _ouput_position,
                      bool set_output)
  {
    return gen_factory(Linear{ .f0 = flow, .df = df },
                       flow,
                       concentration,
                       species_index,
                       input_position,
                       _ouput_position,
                       set_output);
  }

  std::vector<FeedDescriptor>
  FeedFactory::gas_feed(double concentration,
                        std::size_t species_index,
                        std::span<std::size_t> inlet_i,
                        std::span<double> inlet_f,
                        std::span<std::size_t> outlet_i,
                        std::span<double> outlet_f,
                        double tol)
  {
    if (inlet_f.size() != inlet_i.size())
    {
      throw Simulation::FeedException(FeedExceptionError::MismatchSize);
    }

    if (outlet_i.size() != outlet_f.size())
    {
      throw Simulation::FeedException(FeedExceptionError::MismatchSize);
    }

    const double q_in = std::accumulate(inlet_f.begin(), inlet_f.end(), 0.);
    const double q_out = std::accumulate(outlet_f.begin(), outlet_f.end(), 0.);
    if (std::abs(q_in - q_out) > tol * std::max(q_in, q_out))
    {
      throw Simulation::FeedException(FeedExceptionError::UnbalancedFlow);
    }

    std::vector<FeedDescriptor> vi;
    vi.reserve(inlet_i.size() + outlet_i.size());
    for (std::size_t i = 0; i < inlet_i.size(); ++i)
    {
      FeedDescriptor fd = FeedFactory::constant(inlet_f[i],
                                                concentration,
                                                species_index,
                                                inlet_i[i],
                                                std::nullopt,
                                                false);
      // this->add_gas(std::move(fd));
      vi.emplace_back(fd);
    }
    for (std::size_t i = 0; i < outlet_i.size(); ++i)
    {
      FeedDescriptor fd = FeedFactory::sink(outlet_f[i], outlet_i[i]);
      // this->add_gas(std::move(fd));
      vi.emplace_back(fd);
    }

    return vi;
  }

  void
  SimulationFeed::add_feed(FeedDescriptor&& fd, Phase phase)
  {
    auto& vec = phase == Phase::Liquid ? m_liquid : m_gas;
    if (!vec)
    {
      vec = std::vector<FeedDescriptor>();
    }
    KOKKOS_ASSERT(fd.flow >= 0.);

    // Use move_allow_trivial in case
    // FeedDescriptor become non trivial
    vec->emplace_back(move_allow_trivial(fd));
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
    return (m_liquid) ? m_liquid->size() : 0;
  }
  std::size_t
  SimulationFeed::n_gas_flow() const noexcept
  {
    return (m_gas) ? m_liquid->size() : 0;
  }

  SimulationFeed
  SimulationFeed::empty() noexcept
  {
    SimulationFeed s;
    s.m_liquid = std::nullopt;
    s.m_gas = std::nullopt;
    return s;
  }

} // namespace Simulation::Feed
