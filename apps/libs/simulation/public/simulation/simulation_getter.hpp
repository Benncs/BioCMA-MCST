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
    SimulationUnit* const a_;

  public:
    explicit Getter(SimulationUnit* a);
    // explicit Getter(const std::unique_ptr<SimulationUnit>& a);

    Getter(const Getter& m);

    ~Getter();

    [[nodiscard]] const Simulation::Feed::SimulationFeed& get_feed() const;
    [[nodiscard]] double& get_start_time_mut() const noexcept;
    [[nodiscard]] double& get_end_time_mut() noexcept;
    [[nodiscard]] Dimensions getDimensions() const noexcept;
    [[nodiscard]] std::span<double> getCliqData() const;
    [[nodiscard]] std::optional<std::span<const double>> getCgasData() const;
    [[nodiscard]] std::span<const double> getContributionData() const;
    [[nodiscard]] bool two_phase_flow() const noexcept;
    [[nodiscard]] std::optional<std::span<const double>> getMTRData() const;
    [[nodiscard]] std::span<double> getContributionData_mut() const;
    [[nodiscard]] std::ranges::subrange<ProbeIterator, ProbeIterator>
    it_probes() const;

    [[nodiscard]] const std::unique_ptr<MC::MonteCarloUnit>& mc_unit() const;
  };

}; // namespace Simulation

#endif
