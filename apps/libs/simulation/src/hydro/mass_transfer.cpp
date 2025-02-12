#include <cassert>
#include <hydro/impl_mass_transfer.hpp>
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
    _proxy = new MassTransferProxy;
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

  MassTransferModel::MassTransferModel()
      : type(MTRType::Flowmap), _proxy(nullptr),liquid_scalar(nullptr), gas_scalar(nullptr)
  {
  }

  [[nodiscard]] MassTransferProxy* MassTransferModel::proxy() const
  {
    return _proxy;
  }

  MassTransferModel::~MassTransferModel()
  {
    if (_proxy != nullptr)
    {
      delete _proxy;
      _proxy = nullptr;
    }
  }

}; // namespace Simulation::MassTransfer
