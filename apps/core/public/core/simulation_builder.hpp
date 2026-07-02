#ifndef __SIMULATION_BUILDER_HPP__
#define __SIMULATION_BUILDER_HPP__

#include "simulation/mass_transfer.hpp"
#include <common/results.hpp>
#include <mc/unit.hpp>
#include <memory>
#include <mixture/species_descriptor.hpp>
#include <optional>
#include <simulation/feed_descriptor.hpp>
#include <simulation/scalar_initializer.hpp>
#include <simulation/simulation_getter.hpp>

namespace Simulation
{
  class SimulationUnit;
} // namespace Simulation

namespace Core
{

  enum class BuildErrorCode
  {
    MissingUnit,
    MissingScalarInit,
    MissingSpecieTable,
    ErrorSimulation,
  };

  struct BuildError
  {
    BuildErrorCode code;
    std::string message;
  };

  class SimulationUnitBuilder
  {
  public:
    SimulationUnitBuilder() = default;
    ~SimulationUnitBuilder() = default;

    SimulationUnitBuilder(const SimulationUnitBuilder&) = delete;
    SimulationUnitBuilder& operator=(const SimulationUnitBuilder&) = delete;

    SimulationUnitBuilder(SimulationUnitBuilder&&) = default;
    SimulationUnitBuilder& operator=(SimulationUnitBuilder&&) = default;

    SimulationUnitBuilder with_unit(std::unique_ptr<MC::MonteCarloUnit>&& unit);

    SimulationUnitBuilder
    with_scalar(Simulation::ScalarInitializer&& scalar_init);

    SimulationUnitBuilder
    with_specie_table(std::shared_ptr<Mixture::SpecieTable> specie_table);

    SimulationUnitBuilder with_params(Simulation::Sparam params);

    SimulationUnitBuilder with_mt_model(
        std::optional<Simulation::MassTransfer::Type::MtrTypeVariant>&&
            mt_model_variant);

    // Truly optional
    SimulationUnitBuilder with_feed(Simulation::Feed::SimulationFeed feed);

    [[nodiscard]] Result<std::unique_ptr<Simulation::SimulationUnit>,
                         BuildError>
    build();

  private:
    std::unique_ptr<MC::MonteCarloUnit> m_unit;
    std::optional<Simulation::ScalarInitializer> m_scalar_init;
    std::shared_ptr<Mixture::SpecieTable> m_specie_table;

    std::optional<Simulation::MassTransfer::Type::MtrTypeVariant>
        m_mt_model_variant;

    // Default value
    Simulation::Sparam m_params = Simulation::Sparam::fdefault();

    // genuinely optional
    std::optional<Simulation::Feed::SimulationFeed> m_feed;
  };

} // namespace Core

#endif
