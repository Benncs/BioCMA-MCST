#ifndef __SIMULATION_HYDRO_MASS_TRANSFER_HPP__
#define __SIMULATION_HYDRO_MASS_TRANSFER_HPP__

#include <cma_utils/iteration_state.hpp>
#include <memory>
#include <variant>
#include <vector>

namespace Simulation
{
  class ScalarSimulation;
};

namespace Simulation::MassTransfer
{

  namespace Type
  {
    enum class Names
    {
      FlowmapTurb,
      FixedKla,
      FlowmapKla
    };

    struct FlowmapTurbulence
    {
    };

    struct FlowmapKla
    {
    };

    struct FixedKla
    {
      std::vector<double> value;
    };

    using MtrTypeVariant =
        std::variant<FlowmapTurbulence, FixedKla, FlowmapKla>;
  } // namespace Type

  enum class Sign : int
  {
    GasToLiquid = -1,
    LiquidToGas = 1,
  };

  static_assert(static_cast<float>(Sign::GasToLiquid) == -1., "Sign mtr");

  struct MassTransferProxy;

  class MassTransferModel
  {
  public:
    explicit MassTransferModel(
        MassTransfer::Type::MtrTypeVariant _type,
        std::shared_ptr<Simulation::ScalarSimulation> _liquid_scalar,
        std::shared_ptr<Simulation::ScalarSimulation> _gas_scalar);

    void gas_liquid_mass_transfer(const CmaUtils::IterationState& state) const;

    [[nodiscard]] const std::shared_ptr<MassTransferProxy>& proxy() const;

    [[nodiscard]] std::optional<std::span<const double>> mtr_data() const;

    MassTransferModel(MassTransferModel&& rhs) noexcept;
    MassTransferModel& operator=(MassTransferModel&& rhs) noexcept;

    MassTransferModel(const MassTransferModel& rhs) = delete;
    MassTransferModel& operator=(const MassTransferModel& rhs) = delete;

    ~MassTransferModel();

    MassTransferModel();

  private:
    MassTransfer::Type::MtrTypeVariant type;
    std::shared_ptr<MassTransferProxy> _proxy;
    std::shared_ptr<Simulation::ScalarSimulation> liquid_scalar;
    std::shared_ptr<Simulation::ScalarSimulation> gas_scalar;
  };

} // namespace Simulation::MassTransfer

#endif