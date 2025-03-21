#include "common/common.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>
#include <hydro/impl_mass_transfer.hpp>
#include <memory>
#include <optional>
#include <scalar_simulation.hpp>
#include <simulation/mass_transfer.hpp>
#include <stdexcept>
#include <utility>

namespace
{
  struct FunctorKla
  {
    const std::shared_ptr<Simulation::MassTransfer::MassTransferProxy>& proxy;
    std::size_t nrow;

    void operator()(Simulation::MassTransfer::Type::FixedKla& kla) const
    {
      if (kla.value.size() != nrow)
      {
        throw std::invalid_argument("Given kla dimension doesn’t match with CM dimensions");
      }
      for (std::size_t i = 0; i < nrow; ++i)
      {
        proxy->kla.row(EIGEN_INDEX(i)).setConstant(kla.value[i]);
      }
    }

    void operator()(Simulation::MassTransfer::Type::Flowmap&) const
    {
      proxy->kla.setZero();
    }
  };

  struct MtrVisitor
  {
    const std::shared_ptr<Simulation::MassTransfer::MassTransferProxy>& proxy;
    const std::shared_ptr<Simulation::ScalarSimulation>& liquid_scalar;
    const std::shared_ptr<Simulation::ScalarSimulation>& gas_scalar;
    const CmaUtils::IterationState& state;

    void operator()(const Simulation::MassTransfer::Type::FixedKla& _) const
    {
      (void)_;
      Simulation::MassTransfer::Impl::fixed_kla_gas_liquid_mass_transfer(
          *proxy,
          liquid_scalar->getConcentrationArray(),
          gas_scalar->getConcentrationArray(),
          liquid_scalar->getVolume(),
          state);
    }

    void operator()(const Simulation::MassTransfer::Type::Flowmap& _) const
    {
      (void)_;
      Simulation::MassTransfer::Impl::flowmap_gas_liquid_mass_transfer(
          *proxy,
          liquid_scalar->getConcentrationArray(),
          gas_scalar->getConcentrationArray(),
          liquid_scalar->getVolume(),
          state);
    }
  };

} // namespace

namespace Simulation::MassTransfer
{

  MassTransferModel::MassTransferModel(MassTransfer::Type::MtrTypeVariant _type,
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
    _proxy->Henry = Eigen::ArrayXXd(liquid_scalar->n_row(), 1);
    _proxy->Henry.setZero();
    _proxy->Henry(1) = 3.181e-2;

    std::visit(FunctorKla{_proxy, nrow}, _type);

    _proxy->db = 5e-3;
  }

  void MassTransferModel::gas_liquid_mass_transfer(const CmaUtils::IterationState& state) const
  {
    PROFILE_SECTION("gas_liquid_mass_transfer")
    if (gas_scalar == nullptr || _proxy == nullptr)
    {
      throw std::invalid_argument(
          "gas_liquid_mass_transfer should not be called if gas not intialized");
    }
    std::visit(MtrVisitor{_proxy, liquid_scalar, gas_scalar, state}, type);
    // switch (type)
    // {
    // case Simulation::MassTransfer::MTRType::FixedKla:
    // {
    //   Impl::fixed_kla_gas_liquid_mass_transfer(*_proxy,
    //                                            liquid_scalar->getConcentrationArray(),
    //                                            gas_scalar->getConcentrationArray(),
    //                                            liquid_scalar->getVolume(),
    //                                            state);
    //   break;
    // };
    // case Simulation::MassTransfer::MTRType::Flowmap:
    // {
    //   Impl::flowmap_gas_liquid_mass_transfer(*_proxy,
    //                                          liquid_scalar->getConcentrationArray(),
    //                                          gas_scalar->getConcentrationArray(),
    //                                          liquid_scalar->getVolume(),
    //                                          state);
    //   break;
    // }
    // default:
    // {
    //   assert(0 && "gas_liquid_mass_transfer switch");
    //   __builtin_unreachable(); // TODO use c++23 cross plateform unreachable
    // };
    // }
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
      : type(Type::Flowmap{}), _proxy(nullptr), liquid_scalar(nullptr), gas_scalar(nullptr)
  {
  }

  [[nodiscard]] const std::shared_ptr<MassTransferProxy>& MassTransferModel::proxy() const
  {
    return _proxy;
  }

  MassTransferModel::~MassTransferModel() = default;

  MassTransferModel::MassTransferModel(MassTransferModel&& rhs) noexcept = default;

  MassTransferModel& MassTransferModel::operator=(MassTransferModel&& rhs) noexcept = default;

}; // namespace Simulation::MassTransfer
