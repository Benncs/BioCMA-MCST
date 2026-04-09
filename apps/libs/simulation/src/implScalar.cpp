#include "Kokkos_CheckUsage.hpp"
#include "eigen_kokkos.hpp"
#include <Kokkos_Assert.hpp>
#include <algorithm>
#include <mc/alias.hpp>
#include <type_traits>

#ifndef NDEBUG
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#  pragma GCC diagnostic ignored "-Wnan-infinity-disabled"
#endif
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#ifndef NDEBUG
#  pragma GCC diagnostic pop
#endif
#include <Kokkos_Core.hpp>
#include <common/common.hpp>
#include <scalar_simulation.hpp>
#include <simulation/simulation_exception.hpp>
#include <stdexcept>
namespace
{

  template <typename ViewType1, typename ViewType2>
  void
  smrt_deep_copy(ViewType1 dst, ViewType2 src)
  {
    if (dst.extent(0) != src.extent(0) || dst.extent(1) != dst.extent(1))
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
      Kokkos::parallel_for(
          "synchro_source",
          Kokkos::MDRangePolicy<
              typename ViewType1::execution_space,
              Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>(
              { 0, 0 }, { src.extent(0), src.extent(1) }),
          KOKKOS_LAMBDA(int i, int j) { dst(i, j) = src(i, j); });
    }
  }

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

    // KOKKOS_ASSERT(res.cols() == EIGEN_INDEX(n_c));
    // using T = Eigen::Triplet<double>;
    // std::vector<T> triplets(vals.size());

    // triplets.reserve(vals.size());

    // for (std::size_t i = 0; i < vals.size(); ++i)
    // {
    //   triplets[i] = T(rows[i], cols[i], vals[i]);
    // }

    // // res.
    // res.setFromTriplets(triplets.begin(), triplets.end());
    // res.makeCompressed();

    KOKKOS_ASSERT(res.cols() == EIGEN_INDEX(n_c));

    using T = Eigen::Triplet<double>;
    const auto n = vals.size();
    std::vector<T> triplets(n);
    triplets.reserve(n);

    // #pragma omp simd
    for (std::size_t i = 0; i < n; i += 8)
    {
      std::size_t idx0 = i;
      std::size_t idx1 = i + 1;
      std::size_t idx2 = i + 2;
      std::size_t idx3 = i + 3;
      std::size_t idx4 = i + 4;
      std::size_t idx5 = i + 5;
      std::size_t idx6 = i + 6;
      std::size_t idx7 = i + 7;

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
      : n_r(n_species), n_c(n_compartments), concentrations(n_r, n_c),
        sources(n_r, n_c)
  {
    if (volumes.size() != n_compartments)
    {
      throw SimulationException(ErrorCodes::MismatchSizeVolume);
    }

    const int n_row = EIGEN_INDEX(n_r);
    const int n_col = EIGEN_INDEX(n_c);

    m_volumes = DiagonalType(n_col);

    this->m_volumes.diagonal() = Eigen::Map<const Eigen::VectorXd>(
        volumes.data(), static_cast<int>(volumes.size()));

    contribs = kernelContribution("contribs", n_r, n_c);

    m_transition = FlowMatrixType(n_col, n_col);

    volumes_inverse = DiagonalType(n_col);
    volumes_inverse.setIdentity();

    this->total_mass = ColMajorMatrixtype<double>(n_row, n_col);
    this->total_mass.setZero();

    this->sink = DiagonalType(n_col);
    this->sink.setZero();
  }

  // simple getters

  [[nodiscard]] ColMajorMatrixtype<double>&
  ScalarSimulation::get_concentration() noexcept
  {
    return concentrations.eigen_data;
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
    Kokkos::deep_copy(sources.compute, contribs);
  }

  [[nodiscard]] ColMajorKokkosScalarMatrix<double>
  ScalarSimulation::get_device_concentration() const
  {
    // smrt_deep_copy(kc, concentrations.compute);
    // return kc;
    return concentrations.compute;
  }

  void
  ScalarSimulation::reduce_contribs(std::span<const double> data)
  {
    assert(data.size() == (n_c * n_r));
    using eigen_type = decltype(sources)::EigenMatrix;
    sources.eigen_data.noalias() += Eigen::Map<eigen_type>(
        const_cast<double*>(data.data()), EIGEN_INDEX(n_r), EIGEN_INDEX(n_c));
  }
  void
  ScalarSimulation::performStepGL(double d_t,
                                  const ColMajorMatrixtype<double>& mtr,
                                  MassTransfer::Sign sign)
  {
    PROFILE_SECTION("performStep_gl")
#define c concentrations.eigen_data

    total_mass.noalias() += d_t
                            * (c * m_transition - c * sink + sources.eigen_data
                               + static_cast<float>(sign) * mtr);
    c.noalias() = total_mass * volumes_inverse;

    // Make accessible new computed concentration to ComputeSpace
    concentrations.update_host_to_compute();
  }

  void
  ScalarSimulation::performStep(double d_t)
  {
    PROFILE_SECTION("performStep_l")
#define c concentrations.eigen_data

    total_mass.noalias()
        += d_t * (c * m_transition - c * sink + sources.eigen_data);
    c.noalias() = total_mass * volumes_inverse;

    // Make accessible new computed concentration to ComputeSpace
    concentrations.update_host_to_compute();
  }

  void
  ScalarSimulation::clearNegs()
  {
    auto s = concentrations.get_span();
    using float_t = decltype(concentrations)::float_t;
    constexpr float_t TOL = -1e-8;
    using space = decltype(concentrations)::HostView::execution_space;
    auto hv = concentrations.host;
    Kokkos::parallel_for(
        "clear_negs",
        Kokkos::MDRangePolicy<
            space,
            Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>(
            { 0, 0 }, { n_r, n_c }),
        KOKKOS_LAMBDA(int i, int j) {
          auto val = hv(i, j);
          if (val < static_cast<float_t>(0) && val >= TOL)
          {
            hv(i, j) = static_cast<float_t>(0);
          }
        });

    // std::replace_if(
    //     s.begin(), s.end(), [](auto&& f) { return f > TOL && f < 0.; }, 0.);

    concentrations.update_host_to_compute();
  }

  bool
  ScalarSimulation::deep_copy_concentration(const std::vector<double>& data)
  {
    if (data.size() != n_c * n_r)
    {
      return false;
    }
    using eigen_type = decltype(this->concentrations)::EigenMatrix;
    Eigen::Map<const eigen_type> temp_map(
        data.data(), EIGEN_INDEX(n_r), EIGEN_INDEX(n_c));
    this->concentrations.eigen_data = temp_map; // Performs deep copy
    return true;
  }

  void
  ScalarSimulation::set_mass()
  {
    total_mass = this->concentrations.eigen_data * m_volumes;
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
