#include <Kokkos_Assert.hpp>
#include <common/eigen_diag.hpp>
#include <mc/alias.hpp>
#include <type_traits>
EIGEN_DIAG_PUSH
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <kokkos_eigen.hpp>
EIGEN_DIAG_POP
#include <Kokkos_Core.hpp>
#include <common/common.hpp>
#include <scalar_simulation.hpp>
#include <simulation/simulation_exception.hpp>
#include <stdexcept>

namespace
{

  template <typename ViewType1, typename ViewType2>
    requires(ViewType1::rank() == ViewType2::rank()
             && std::is_convertible_v<typename ViewType2::non_const_value_type,
                                      typename ViewType1::non_const_value_type>)
  void
  smrt_deep_copy2d(ViewType1 dst, const ViewType2& src)
  {

    if (dst.extent(0) != src.extent(0) || dst.extent(1) != src.extent(1))
    {
      throw std::runtime_error("Dimension mismatch");
    }

    if constexpr (std::is_same_v<typename ViewType1::array_layout,
                                 typename ViewType2::array_layout>
                  && std::is_same_v<typename ViewType1::non_const_value_type,
                                    typename ViewType2::non_const_value_type>)
    {
      Kokkos::deep_copy(dst, src);
    }
    else
    {
      using exec_space = typename ViewType1::execution_space;
      auto src_view = Kokkos::create_mirror_view_and_copy(exec_space(), src);
      Kokkos::parallel_for(
          "smrt_deep_copy",
          Kokkos::MDRangePolicy<
              typename ViewType1::execution_space,
              Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>(
              { 0, 0 }, { src.extent(0), src.extent(1) }),
          KOKKOS_LAMBDA(int i, int j) { dst(i, j) = src_view(i, j); });
    }
  }

  // void
  // sparse_from_coo(Eigen::SparseMatrix<double>& res,
  //                 const std::size_t n_c,
  //                 std::span<const size_t>&& rows,
  //                 std::span<const size_t>&& cols,
  //                 std::span<const double>&& vals)
  // {
  //   KOKKOS_ASSERT(res.cols() == EIGEN_INDEX(n_c));
  //   using T = Eigen::Triplet<double>;
  //   const auto n = vals.size();
  //   KOKKOS_ASSERT(rows.size() == n);
  //   KOKKOS_ASSERT(cols.size() == n);
  //   KOKKOS_ASSERT(vals.size() == n);
  //   std::vector<T> triplets(n);

  //   for (std::size_t i = 0; i < n; i++)
  //   {
  //     triplets[i]
  //         = T(static_cast<int>(rows[i]), static_cast<int>(cols[i]), vals[i]);
  //   }
  //   res.setFromTriplets(triplets.begin(), triplets.end());
  //   res.makeCompressed();
  // }

#define UNROLL_CHECK(i, s)                                                     \
  if ((i) < (s))                                                               \
  {                                                                            \
    triplets[i] = T(rows[i], cols[i], vals[i]);                                \
  }

#define UNROLL(idx) triplets[idx] = T(rows[idx], cols[idx], vals[idx]);
  void
  sparse_from_coo(Eigen::SparseMatrix<double>& res,
                  std::size_t n_c,
                  std::span<const size_t>&& rows,
                  std::span<const size_t>&& cols,
                  std::span<const double>&& vals)
  {

    KOKKOS_ASSERT(res.cols() == EIGEN_INDEX(n_c));

    using T = Eigen::Triplet<double>;
    const auto n = vals.size();
    std::vector<T> triplets(n);

    // #pragma omp simd
    for (std::size_t i = 0; i < n; i += 8)
    {
      const std::size_t idx0 = i;
      const std::size_t idx1 = i + 1;
      const std::size_t idx2 = i + 2;
      const std::size_t idx3 = i + 3;
      const std::size_t idx4 = i + 4;
      const std::size_t idx5 = i + 5;
      const std::size_t idx6 = i + 6;
      const std::size_t idx7 = i + 7;

      if (idx7 < n)
      {
        UNROLL(idx0);
        UNROLL(idx1);
        UNROLL(idx2);
        UNROLL(idx3);
        UNROLL(idx4);
        UNROLL(idx5);
        UNROLL(idx6);
        UNROLL(idx7);
      }
      else
      {
        UNROLL_CHECK(idx0, n)
        UNROLL_CHECK(idx1, n)
        UNROLL_CHECK(idx2, n)
        UNROLL_CHECK(idx3, n)
        UNROLL_CHECK(idx4, n)
        UNROLL_CHECK(idx5, n)
        UNROLL_CHECK(idx6, n)
      }
    }
    res.setFromTriplets(triplets.begin(), triplets.end());
    res.makeCompressed();
  }

} // namespace

namespace Simulation
{

  ScalarSimulation::ScalarSimulation(size_t n_compartments,
                                     size_t n_species,
                                     std::span<double> volumes)
      : n_r(n_species), n_c(n_compartments), total_mass(n_r, n_c),
        concentrations("concentrations", n_r, n_c),
        sources("sources", n_r, n_c), contribs("contribs", n_r, n_c),
        m_transition(
            FlowMatrixType<double>(EIGEN_INDEX(n_c), EIGEN_INDEX(n_c))),
        m_volumes(KokkosEigen::Alias::DiagonalType<mass_balance_float_type>(
            EIGEN_INDEX(n_c))),
        volumes_inverse(
            KokkosEigen::Alias::DiagonalType<mass_balance_float_type>(
                EIGEN_INDEX(n_c))),
        sink(KokkosEigen::Alias::DiagonalType<mass_balance_float_type>(
            EIGEN_INDEX(n_c)))
  {
    if (volumes.size() != n_compartments)
    {
      throw SimulationException(ErrorCodes::MismatchSizeVolume);
    }

    this->m_volumes.diagonal() = Eigen::Map<const Eigen::VectorXd>(
        volumes.data(), static_cast<int>(volumes.size()));

    this->volumes_inverse.setIdentity();

    this->total_mass.setZero();

    this->sink.setZero();
  }

  // simple getters

  [[nodiscard]] ScalarSimulation::concentration_view_t
  ScalarSimulation::get_concentration() noexcept
  {
    return concentrations.host_view();
  }

  std::size_t
  ScalarSimulation::n_col() const noexcept
  {
    return n_c;
  }

  std::size_t
  ScalarSimulation::n_row() const noexcept
  {
    return n_r;
  }

  void
  ScalarSimulation::synchro_sources()
  {

    static_assert(Kokkos::SpaceAccessibility<
                      decltype(contribs)::execution_space::memory_space,
                      decltype(sources)::device_view_type::execution_space::
                          memory_space>::accessible
                  != 0U);

    Kokkos::deep_copy(sources.device_view(), contribs);
    sources.device_to_host_sync();
  }

  [[nodiscard]] MC::KernelConcentrationType
  ScalarSimulation::get_device_concentration() const
  {
    return concentrations.device_view_cst();
  }

  // void
  // ScalarSimulation::reduce_contribs(std::span<const double> data)
  // {
  //   assert(data.size() == (n_c * n_r));
  //   using eigen_map_type = decltype(sources)::eigen_map_type;

  //   sources.eigen().noalias() += eigen_map_type(
  //       const_cast<double*>(data.data()), EIGEN_INDEX(n_r),
  //       EIGEN_INDEX(n_c));

  //   // sources.eigen_data.noalias() += Eigen::Map<eigen_type>(
  //   //     const_cast<double*>(data.data()), EIGEN_INDEX(n_r),
  //   //     EIGEN_INDEX(n_c));
  // }

  void
  ScalarSimulation::performStepGL(
      double d_t,
      const KokkosEigen::Alias::ColMajorMatrixtype<double>& mtr,
      MassTransfer::Sign sign)
  {
    PROFILE_SECTION("performStep_gl")

    auto& c = concentrations.eigen();
    const auto& _sources = sources.eigen();

    auto dmdt = c * m_transition - c * sink + _sources
                + static_cast<float>(sign) * mtr;

    total_mass.noalias() += d_t * dmdt;

    c.noalias() = total_mass * volumes_inverse;

    // Make accessible new computed concentration to ComputeSpace
    concentrations.host_to_device_sync();
  }

  void
  ScalarSimulation::performStep(double d_t)
  {
    PROFILE_SECTION("performStep_l")

    auto& c = concentrations.eigen();
    const auto& _sources = sources.cst_eigen();
    auto dmdt = c * m_transition - c * sink + _sources;

    total_mass.noalias() += d_t * dmdt;
    c.noalias() = total_mass * volumes_inverse;

    // Make accessible new computed concentration to ComputeSpace
    concentrations.host_to_device_sync();
  }

  void
  ScalarSimulation::clearNegs()
  {

    using float_type = decltype(concentrations)::float_type;
    // constexpr float_type TOL = -1e-8;
    // order of magniture species to clip
    constexpr float_type max_species_value = 5e-3;
    // O(d_t)=1e-4
    constexpr float_type scheme_relative_error = 1e-4;

    constexpr float_type TOL = scheme_relative_error * max_species_value;

    using space = decltype(concentrations)::host_view_type::execution_space;
    auto hv = concentrations.host_view();
    Kokkos::parallel_for(
        "clear_negs",
        Kokkos::MDRangePolicy<
            space,
            Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>(
            { 0, 0 }, { n_r, n_c }),
        KOKKOS_LAMBDA(int i, int j) {
          const auto val = hv(i, j);
          if (val < static_cast<float_type>(0) && Kokkos::abs(val) < TOL)
          {
            hv(i, j) = static_cast<float_type>(0);
          }
        });

    // concentrations.update_host_to_compute();
    concentrations.host_to_device_sync();
  }

  bool
  ScalarSimulation::deep_copy_concentration(const std::vector<double>& data)
  {
    if (data.size() != n_c * n_r)
    {
      return false;
    }
    using view_layout = decltype(this->concentrations)::view_layout;
    Kokkos::View<const double**, view_layout> unmanaged_host_view(
        data.data(), n_r, n_c);

    Kokkos::deep_copy(concentrations.host_view(), unmanaged_host_view);
    concentrations.host_to_device_sync();

    return true;
  }

  void
  ScalarSimulation::set_mass()
  {
    total_mass = this->concentrations.eigen() * m_volumes;
  }

  void
  ScalarSimulation::set_transition(CmaUtils::StateCooMatrixType&& transition)
  {
    sparse_from_coo(this->m_transition,
                    transition->nrows(),
                    transition->row_indices(),
                    transition->col_indices(),
                    transition->values());
  }

} // namespace Simulation
