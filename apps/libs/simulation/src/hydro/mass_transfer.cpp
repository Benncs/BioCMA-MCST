#include <cassert>
#include <hydro/impl_mass_transfer.hpp>
#include <memory>
#include <optional>
#include <scalar_simulation.hpp>
#include <simulation/mass_transfer.hpp>
#include <stdexcept>
#include <utility>

namespace Simulation::MassTransfer
{

  MassTransferModel::MassTransferModel(MTRType _type,
                                       std::shared_ptr<Simulation::ScalarSimulation> _liquid_scalar,
                                       std::shared_ptr<Simulation::ScalarSimulation> _gas_scalar)
      : type(_type), liquid_scalar(std::move(_liquid_scalar)), gas_scalar(std::move(_gas_scalar))
  {

    const auto nrow = liquid_scalar->n_row();
    const auto ncol = liquid_scalar->n_col();

    // _proxy = // NOLINT
    //     new MassTransferProxy{MatrixType(nrow, ncol), Eigen::ArrayXXd(nrow, ncol)};
    _proxy = std::make_shared<MassTransferProxy>();
    _proxy->mtr = MatrixType(nrow, ncol);
    _proxy->kla = Eigen::ArrayXXd(nrow, ncol);

    _proxy->kla.setZero();
  }

  void MassTransferModel::gas_liquid_mass_transfer(const CmaRead::ReactorState& state) const
  {
    PROFILE_SECTION("gas_liquid_mass_transfer")
    if (gas_scalar == nullptr || _proxy == nullptr)
    {
      throw std::invalid_argument(
          "gas_liquid_mass_transfer should not be called if gas not intialized");
    }

    switch (type)
    {
    case Simulation::MassTransfer::MTRType::FixedKla:
    {
      break;
    };
    case Simulation::MassTransfer::MTRType::Flowmap:
    {
      Impl::flowmap_gas_liquid_mass_transfer(*_proxy, liquid_scalar, gas_scalar, state);
      break;
    }
    default:
    {
      assert(0 && "gas_liquid_mass_transfer switch");
      __builtin_unreachable(); // TODO use c++23 cross plateform unreachable
    };
    }
  }

  std::optional<std::span<const double>> MassTransferModel::mtr_data() const
  {
    if (_proxy != nullptr)
    {
      return std::make_optional<std::span<const double>>(
          {_proxy->mtr.data(), static_cast<size_t>(_proxy->mtr.size())});
    }
    return std::nullopt;
  }

  MassTransferModel::MassTransferModel()
      : type(MTRType::Flowmap), _proxy(nullptr), liquid_scalar(nullptr), gas_scalar(nullptr)
  {
  }

  [[nodiscard]] const std::shared_ptr<MassTransferProxy>& MassTransferModel::proxy() const
  {
    return _proxy;
  }

  MassTransferModel::~MassTransferModel()
  {
  }

  MassTransferModel::MassTransferModel(MassTransferModel&& rhs) noexcept
      : type(rhs.type), _proxy(std::move(rhs._proxy)), liquid_scalar(std::move(rhs.liquid_scalar)),
        gas_scalar(std::move(rhs.gas_scalar))
  {
  }

  MassTransferModel& MassTransferModel::operator=(MassTransferModel&& rhs) noexcept
  {
    if (this != &rhs)
    {
      type = rhs.type;
      _proxy = std::move(rhs._proxy);
      liquid_scalar = std::move(rhs.liquid_scalar);
      gas_scalar = std::move(rhs.gas_scalar);

      // rhs._proxy = nullptr;
    }
    return *this;
  }

}; // namespace Simulation::MassTransfer
