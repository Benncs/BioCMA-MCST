#include <common/eigen_diag.hpp>
EIGEN_DIAG_PUSH
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
EIGEN_DIAG_POP

#include <cassert>
#include <common/common.hpp>
#include <hydro/impl_mass_transfer.hpp>
#include <memory>
#include <mixture/species_descriptor.hpp>
#include <optional>
#include <scalar.hpp>
#include <simulation/mass_transfer.hpp>
#include <stdexcept>
#include <utility>

namespace
{
  struct FunctorKla
  {
    const std::shared_ptr<Simulation::MassTransfer::MassTransferProxy>& proxy;
    std::size_t nrow;

    void
    operator()(Simulation::MassTransfer::Type::FixedKla& kla) const
    {
      if (kla.value.size() != nrow)
      {
        throw std::invalid_argument(
            "Given kla dimension doesn’t match with CM dimensions");
      }
      for (std::size_t i = 0; i < nrow; ++i)
      {
        proxy->kla.row(EIGEN_INDEX(i)).setConstant(kla.value[i]);
      }

      // proxy->kla(1,0)=0;
    }

    void
    operator()(Simulation::MassTransfer::Type::FlowmapTurbulence&) const
    {
      proxy->kla.setZero();
    }
    void
    operator()(Simulation::MassTransfer::Type::FlowmapKla&) const
    {
      proxy->kla.setZero();
    }
    void
    operator()(Simulation::MassTransfer::Type::Auto&) const
    {
      for (std::size_t i = 0; i < nrow; ++i)
      {
        // proxy->kla.row(EIGEN_INDEX(i)).setConstant(0.2); // 700h/1

        proxy->kla.row(EIGEN_INDEX(i)).setConstant(0.011); // 700h/1
      }
    }
  };

  struct MtrVisitor
  {
    const std::shared_ptr<Simulation::MassTransfer::MassTransferProxy>& proxy;
    const std::shared_ptr<Simulation::impl::ScalarSimulation>& scalar;
    const CmaUtils::IterationStatePtrType& state;

    void
    operator()(const Simulation::MassTransfer::Type::Auto& _) const
    {
      (void)_;
      Simulation::MassTransfer::Impl::fixed_kla_gas_liquid_mass_transfer(
          *proxy,
          scalar->phase(Phase::Liquid).concentration.array(),
          scalar->phase(Phase::Gas).concentration.array(),
          scalar->phase(Phase::Liquid).volume.asDiagonal(),
          state);
    }

    void
    operator()(const Simulation::MassTransfer::Type::FixedKla& _) const
    {
      (void)_;
      Simulation::MassTransfer::Impl::fixed_kla_gas_liquid_mass_transfer(
          *proxy,
          scalar->phase(Phase::Liquid).concentration.array(),
          scalar->phase(Phase::Gas).concentration.array(),
          scalar->phase(Phase::Liquid).volume.asDiagonal(),
          state);
    }

    void
    operator()(
        const Simulation::MassTransfer::Type::FlowmapTurbulence& variant) const
    {
      const double db = variant.db;

      Simulation::MassTransfer::Impl::flowmap_gas_liquid_mass_transfer(
          *proxy,
          db,
          scalar->phase(Phase::Liquid).concentration.array(),
          scalar->phase(Phase::Gas).concentration.array(),
          scalar->phase(Phase::Liquid).volume.asDiagonal(),
          state);
    }

    void
    operator()(const Simulation::MassTransfer::Type::FlowmapKla& _) const
    {
      (void)_;
      // Simulation::MassTransfer::Impl::flowmap_gas_liquid_mass_transfer(
      //     *proxy,
      //     liquid_scalar->getConcentrationArray(),
      //     gas_scalar->getConcentrationArray(),
      //     liquid_scalar->getVolume(),
      //     state);
    }
  };

} // namespace

namespace Simulation::MassTransfer
{

  namespace Type
  {
    struct MtrVariantName
    {
      std::string
      operator()(const Auto&) const
      {
        return "auto mtr";
      }

      std::string
      operator()(const FixedKla&) const
      {
        return "Fixed kla mtr";
      }

      std::string
      operator()(const FlowmapKla&) const
      {
        return "Flowmap kla mtr";
      }

      std::string
      operator()(const FlowmapTurbulence& v) const
      {
        return "Flowmap turb with db=" + std::to_string(v.db);
      }
    };

    std::string
    mtr_variant_info(Simulation::MassTransfer::Type::MtrTypeVariant& v) noexcept
    {
      return std::visit(MtrVariantName{}, v);
    }
  } // namespace Type

  MassTransferModel::MassTransferModel(
      const Mixture::SpecieTable& species,
      MassTransfer::Type::MtrTypeVariant _type,
      std::shared_ptr<Simulation::impl::ScalarSimulation> _scalar)
      : type(_type), scalar(std::move(_scalar))
  {

    const auto nrow = scalar->n_species();      // nspecies
    const auto ncol = scalar->n_compartments(); // n compartment

    _proxy = std::make_shared<MassTransferProxy>();
    _proxy->mtr = KokkosEigen::Alias::ColMajorMatrixtype<double>(nrow, ncol);
    _proxy->kla = Eigen::ArrayXXd(nrow, ncol);
    _proxy->flag_transfer = Eigen::ArrayXXd(nrow, 1);

    const auto henry = species.henry();

    _proxy->Henry.resize(EIGEN_INDEX(henry.size()), 1);

    std::visit(FunctorKla{ _proxy, nrow }, _type);
    // Fixme, how to desactivate transfer
    // Set kla to 0
    int i = 0;
    for (const auto& h : henry)
    {
      _proxy->flag_transfer.coeffRef(i, 0) = (h == 0.) ? 0. : 1.;

      // if (h == 0.)
      // {
      //   _proxy->kla.row(EIGEN_INDEX(i)).setConstant(0.);
      // }
      _proxy->Henry.coeffRef(i, 0) = h;
      i++;
    }

    // _proxy->db = 5e-3; // FIXME
  }

  void
  MassTransferModel::update(const CmaUtils::IterationStatePtrType& state)
  {
    PROFILE_SECTION("gas_liquid_mass_transfer")
    if (scalar == nullptr || !scalar->has_gas() || _proxy == nullptr)
    {
      throw std::invalid_argument("gas_liquid_mass_transfer should not be "
                                  "called if gas not intialized");
    }
    std::visit(MtrVisitor{ _proxy, scalar, state }, type);
  }

  void
  MassTransferModel::gas_liquid_mass_transfer() const
  {
    PROFILE_SECTION("gas_liquid_mass_transfer")
    if (scalar == nullptr || !scalar->has_gas() || _proxy == nullptr)
    {
      throw std::invalid_argument("gas_liquid_mass_transfer should not be "
                                  "called if gas not intialized");
    }
    // std::visit(MtrVisitor{_proxy, liquid_scalar, gas_scalar, state}, type);

    const auto& liquid = this->scalar->phase(Phase::Liquid);
    const auto& gas = this->scalar->phase(Phase::Gas);

    auto gas_concentration = gas.concentration.array();
    auto liquid_concentration = liquid.concentration.array();
    auto liquid_volume = liquid.volume.asDiagonal();

    auto& flag = this->_proxy->flag_transfer;

    // _proxy->mtr = (_proxy->kla
    //                * (gas_concentration.colwise() * _proxy->Henry
    //                   - liquid_concentration))
    //                   .matrix()
    //               * liquid_volume ;

    // _proxy->mtr = ((_proxy->kla
    //                 * (gas_concentration.colwise() * _proxy->Henry
    //                    - liquid_concentration))
    //                    .matrix()
    //                * liquid_volume)
    //                   .array()
    //                   .colwise()
    //               * flag.array();

    // Using temp variable shouln't theoretically introduce overhead because of
    // Eigen laziness

    auto concentration_diff
        = gas_concentration.colwise() * _proxy->Henry - liquid_concentration;

    auto flux = (_proxy->kla * concentration_diff).matrix() * liquid_volume;

    _proxy->mtr = flux.array().colwise() * flag.array();
  }

  std::optional<std::span<const double>>
  MassTransferModel::mtr_data() const
  {
    if (_proxy != nullptr)
    {
      return std::make_optional<std::span<const double>>(
          { _proxy->mtr.data(), static_cast<size_t>(_proxy->mtr.size()) });
    }
    return std::nullopt;
  }

  MassTransferModel::MassTransferModel()
      : type(Type::FlowmapTurbulence{}), _proxy(nullptr), scalar(nullptr)
  {
  }

  [[nodiscard]] const std::shared_ptr<MassTransferProxy>&
  MassTransferModel::proxy() const
  {
    return _proxy;
  }

  MassTransferModel::~MassTransferModel() = default;

  MassTransferModel::MassTransferModel(MassTransferModel&& rhs) noexcept
      = default;

  MassTransferModel&
  MassTransferModel::operator=(MassTransferModel&& rhs) noexcept
      = default;

}; // namespace Simulation::MassTransfer
