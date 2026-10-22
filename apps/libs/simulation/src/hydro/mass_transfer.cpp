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
    const CmaUtils::IterationStatePtrType& state;

    void
    operator()(const Simulation::MassTransfer::Type::Auto& _) const
    {
      (void)_;
      Simulation::MassTransfer::Impl::fixed_kla_gas_liquid_mass_transfer(*proxy,
                                                                         state);
    }

    void
    operator()(const Simulation::MassTransfer::Type::FixedKla& _) const
    {
      (void)_;
      Simulation::MassTransfer::Impl::fixed_kla_gas_liquid_mass_transfer(*proxy,
                                                                         state);
    }

    void
    operator()(
        const Simulation::MassTransfer::Type::FlowmapTurbulence& variant) const
    {
      const double db = variant.db;

      Simulation::MassTransfer::Impl::flowmap_gas_liquid_mass_transfer(
          *proxy, db, state);
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

  // Warning with lifetime
  // mtrInput's lifetime < this
  // But the views inside of input are located in the scalar object and have
  // enough lifetime then this is correct to map the data only once in this
  // constructor
  MassTransferModel::MassTransferModel(const Mixture::SpecieTable& species,
                                       MassTransfer::Type::MtrTypeVariant _type,
                                       const MtrInput& input)
      : type(_type)
  {

    const auto nrow = input.liquid_concentration.extent(0); // nspecies
    const auto ncol = input.liquid_concentration.extent(1); // n compartment

    // Views on the scalar
    const auto liquid_concentration
        = Eigen::Map<const Eigen::ArrayXXd>(input.liquid_concentration.data(),
                                            EIGEN_INDEX(nrow),
                                            EIGEN_INDEX(ncol));
    const auto gas_concentration = Eigen::Map<const Eigen::ArrayXXd>(
        input.gas_concentration.data(), EIGEN_INDEX(nrow), EIGEN_INDEX(ncol));
    const auto liquid_volume = Eigen::Map<const Eigen::VectorXd>(
        input.liquid_volume.data(), EIGEN_INDEX(ncol));

    _proxy = std::make_shared<MassTransferProxy>(
        liquid_concentration, gas_concentration, liquid_volume);
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
  }

  void
  MassTransferModel::update(const CmaUtils::IterationStatePtrType& state)
  {
    PROFILE_SECTION("gas_liquid_mass_transfer")
    if (_proxy == nullptr)
    {
      throw std::invalid_argument("gas_liquid_mass_transfer should not be "
                                  "called if gas not intialized");
    }
    std::visit(MtrVisitor{ _proxy, state }, type);
  }

  void
  MassTransferModel::gas_liquid_mass_transfer() const
  {
    PROFILE_SECTION("gas_liquid_mass_transfer")
    if (_proxy == nullptr)
    {
      throw std::invalid_argument("gas_liquid_mass_transfer should not be "
                                  "called if gas not intialized");
    }
    // std::visit(MtrVisitor{_proxy, liquid_scalar, gas_scalar, state}, type);

    const auto& gas_concentration = _proxy->gas_concentration;
    const auto& liquid_concentration = _proxy->liquid_concentration;
    auto liquid_volume = _proxy->liquid_volume.asDiagonal();

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
      : type(Type::FlowmapTurbulence{}), _proxy(nullptr)
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
