#ifndef __SIMULATIONS_UNIT_GETTER_HPP__
#define __SIMULATIONS_UNIT_GETTER_HPP__
#include <mc/unit.hpp>
#include <optional>
#include <simulation/feed_descriptor.hpp>
#include <simulation/probe.hpp>
#include <span>

namespace Simulation
{
  class SimulationUnit;
  struct Dimensions;

  class Getter
  {
  private:
    SimulationUnit const* a_;

  public:
    explicit Getter(SimulationUnit const* a);

    Getter(const Getter& m);

    Getter& operator=(const Getter& m) = delete;
    Getter& operator=(Getter&& m) = delete;

    ~Getter();

    [[nodiscard]] double start_time() const noexcept;
    [[nodiscard]] double endtime() const noexcept;
    [[nodiscard]] bool two_phase_flow() const noexcept;
    [[nodiscard]] Dimensions getDimensions() const noexcept;

    [[nodiscard]] const std::unique_ptr<MC::MonteCarloUnit>&
    mc_unit() const noexcept;
    [[nodiscard]] const Simulation::Feed::SimulationFeed&
    get_feed() const noexcept;

    [[nodiscard]] std::span<const double> getCliqData() const;
    [[nodiscard]] std::optional<std::span<const double>> getCgasData() const;
    [[nodiscard]] std::span<const double> getContributionData() const;
    [[nodiscard]] std::optional<std::span<const double>> getMTRData() const;

    [[nodiscard]] std::span<double> getCliqData_mut() const;
    [[nodiscard]] std::span<double> getContributionData_mut() const;
    [[nodiscard]] std::ranges::subrange<ProbeIterator, ProbeIterator>
    it_probes() const;
  };

}; // namespace Simulation

#endif
