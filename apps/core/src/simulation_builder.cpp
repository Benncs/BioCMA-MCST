#include "common/logger.hpp"
#include "simulation/simulation_exception.hpp"
#include <core/simulation_builder.hpp>
#include <iostream>
#include <mc/unit.hpp>
#include <mixture/species_descriptor.hpp>
#include <simulation/scalar_initializer.hpp>
#include <simulation/simulation.hpp>
namespace Core
{

  SimulationUnitBuilder
  SimulationUnitBuilder::with_unit(std::unique_ptr<MC::MonteCarloUnit>&& unit)
  {
    m_unit = std::move(unit);
    return std::move(*this);
  }

  SimulationUnitBuilder
  SimulationUnitBuilder::with_scalar(
      Simulation::ScalarInitializer&& scalar_init)
  {

    m_scalar_init = std::move(scalar_init);
    return std::move(*this);
  }

  SimulationUnitBuilder
  SimulationUnitBuilder::with_specie_table(
      std::shared_ptr<Mixture::SpecieTable> specie_table)
  {
    m_specie_table = std::move(specie_table);
    return std::move(*this);
  }

  SimulationUnitBuilder
  SimulationUnitBuilder::with_params(Simulation::Sparam params)
  {
    m_params = params;
    return std::move(*this);
  }

  SimulationUnitBuilder
  SimulationUnitBuilder::with_mt_model(
      std::optional<Simulation::MassTransfer::Type::MtrTypeVariant>&&
          mt_model_variant)
  {
    m_mt_model_variant = std::move(mt_model_variant);
    return std::move(*this);
  }

  // Truly optional
  SimulationUnitBuilder
  SimulationUnitBuilder::with_feed(Simulation::Feed::SimulationFeed feed)
  {
    m_feed = std::move(feed);
    return std::move(*this);
  }

  [[nodiscard]] Result<std::unique_ptr<Simulation::SimulationUnit>, BuildError>
  SimulationUnitBuilder::build()
  {
    using ret_type
        = Result<std::unique_ptr<Simulation::SimulationUnit>, BuildError>;
    if (m_unit == nullptr)
    {
      BuildError error
          = BuildError{ BuildErrorCode::MissingUnit,
                        "SimulationUnitBuilder: MonteCarloUnit not set (call "
                        "with_unit())" };
      return ret_type(std::move(error));
    }
    if (!m_scalar_init.has_value())
    {
      auto error
          = BuildError{ BuildErrorCode::MissingScalarInit,
                        "SimulationUnitBuilder: ScalarInitializer not set "
                        "(call with_scalar_init())" };
      return ret_type(std::move(error));
    }
    if (m_specie_table == nullptr)
    {
      auto error
          = BuildError{ BuildErrorCode::MissingSpecieTable,
                        "SimulationUnitBuilder: SpecieTable not set (call "
                        "with_specie_table())" };
      return ret_type(std::move(error));
    }

    try
    {
      return ret_type(std::make_unique<Simulation::SimulationUnit>(
          std::move(m_unit),
          std::move(*m_scalar_init),
          std::move(m_feed),
          std::move(m_specie_table),
          std::move(m_mt_model_variant),
          m_params));
    }
    catch (Simulation::BaseSimulationException& e)
    {
      auto error = BuildError{ BuildErrorCode::ErrorSimulation,
                               IO::format("SimulationUnitBuilder: From ",
                                          e.where().function_name(),
                                          ": ",
                                          e.what()) };
      return ret_type(std::move(error));
    }
  }

} // namespace Core
