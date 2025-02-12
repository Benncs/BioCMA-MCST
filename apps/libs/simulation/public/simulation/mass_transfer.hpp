#ifndef __SIMULATION_HYDRO_MASS_TRANSFER_HPP__
#define __SIMULATION_HYDRO_MASS_TRANSFER_HPP__

#include <cma_read/reactorstate.hpp>
#include <memory>

namespace Simulation
{
  class ScalarSimulation;
};

namespace Simulation::MassTransfer
{
  enum class MTRType
  {
    Flowmap,
    FixedKla,
  };

  struct MassTransferProxy;

  class MassTransferModel
  {
  public:
    explicit MassTransferModel(MassTransfer::MTRType _type,
                               std::shared_ptr<Simulation::ScalarSimulation> _liquid_scalar,
                               std::shared_ptr<Simulation::ScalarSimulation> _gas_scalar);

    void gas_liquid_mass_transfer(const CmaRead::ReactorState& state) const;

    [[nodiscard]] MassTransferProxy* proxy() const;

    MassTransferModel(MassTransferModel&& rhs) noexcept
        : type(rhs.type), _proxy(rhs._proxy), liquid_scalar(std::move(rhs.liquid_scalar)),
          gas_scalar(std::move(rhs.gas_scalar))
    {
      rhs._proxy = nullptr; 
    }

    MassTransferModel& operator=(MassTransferModel&& rhs) noexcept
    {
      if (this != &rhs)
      { 
        type = rhs.type;
        _proxy = rhs._proxy;
        liquid_scalar = std::move(rhs.liquid_scalar);
        gas_scalar = std::move(rhs.gas_scalar);

        rhs._proxy = nullptr; 
      }
      return *this;
    }

    // MassTransferModel(const MassTransferModel& rhs) = default;
    // MassTransferModel& operator=(const MassTransferModel& rhs) = default;
    ~MassTransferModel();

    MassTransferModel();

  private:
    MassTransfer::MTRType type;
    MassTransferProxy* _proxy;
    std::shared_ptr<Simulation::ScalarSimulation> liquid_scalar;
    std::shared_ptr<Simulation::ScalarSimulation> gas_scalar;
  };

} // namespace Simulation::MassTransfer

#endif