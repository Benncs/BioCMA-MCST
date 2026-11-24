#ifndef __SIMULATION_HYDRO_MASS_TRANSFER_HPP__
#define __SIMULATION_HYDRO_MASS_TRANSFER_HPP__

#include <cma_utils/alias.hpp>
#include <memory>
#include <mixture/species_descriptor.hpp>
#include <variant>
#include <vector>

namespace Simulation
{
  class ScalarSimulation;
}; // namespace Simulation

namespace Simulation::MassTransfer::Type
{
  enum class Names
  {
    FlowmapTurb,
    FixedKla,
    FlowmapKla
  };

  struct FlowmapTurbulence
  {
    static constexpr double default_bubble_diameter = 5e-3;
    double db = default_bubble_diameter;
  };

  struct FlowmapKla
  {
  };

  struct FixedKla
  {
    std::vector<double> value;
  };

  struct Auto
  {
  };

  using MtrTypeVariant
      = std::variant<FlowmapTurbulence, FixedKla, FlowmapKla, Auto>;

  std::string mtr_variant_info(MtrTypeVariant& v) noexcept;

} // namespace Simulation::MassTransfer::Type

namespace Simulation::MassTransfer
{

  enum class Sign : int
  {
    GasToLiquid = -1,
    LiquidToGas = 1,
  };

  struct MassTransferProxy;

  class MassTransferModel
  {
  public:
    explicit MassTransferModel(
        const Mixture::SpecieTable& species,
        MassTransfer::Type::MtrTypeVariant _type,
        std::shared_ptr<Simulation::ScalarSimulation> _liquid_scalar,
        std::shared_ptr<Simulation::ScalarSimulation> _gas_scalar);

    void gas_liquid_mass_transfer() const;

    void update(const CmaUtils::IterationStatePtrType& state);

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

  static_assert(static_cast<float>(Sign::GasToLiquid) == -1., "Sign mtr");

} // namespace Simulation::MassTransfer

#endif
