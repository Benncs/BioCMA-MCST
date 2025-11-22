#include "eigen_kokkos.hpp"
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
#include <simulation/alias.hpp>
#include <stdexcept>

namespace
{

#define UNROLL_CHECK(i, s)                                                     \
  if (i < s)                                                                   \
  {                                                                            \
    triplets[i] = T(rows[i], cols[i], vals[i]);                                \
  }

#define UNROLL(idx) triplets[idx] = T(rows[idx], cols[idx], vals[idx]);

  void sparse_from_coo(Eigen::SparseMatrix<double>& res,
                       std::size_t n_c,
                       std::span<const size_t> rows,
                       std::span<const size_t> cols,
                       std::span<const double> vals)
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

    // After filling the triplets, create the sparse matrix
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
      throw std::invalid_argument("Volumes size mismatch");
    }

    const int n_row = EIGEN_INDEX(n_r);
    const int n_col = EIGEN_INDEX(n_c);

    m_volumes = Eigen::DiagonalMatrix<double, -1>(n_col);

    this->m_volumes.diagonal() = Eigen::Map<const Eigen::VectorXd>(
        volumes.data(), static_cast<int>(volumes.size()));

    contribs = kernelContribution("contribs", n_r, n_c);

    m_transition = FlowMatrixType(n_c, n_c);

    volumes_inverse = Eigen::DiagonalMatrix<double, -1>(n_col);
    volumes_inverse.setIdentity();

    this->total_mass = ColMajorMatrixtype<double>(n_row, n_col);
    this->total_mass.setZero();

    this->sink = DiagonalType(n_col);
    this->sink.setZero();
  }

  void ScalarSimulation::synchro_sources()
  {
    Kokkos::deep_copy(sources.compute, contribs);
  }

  [[nodiscard]] ColMajorMatrixtype<double>&
  ScalarSimulation::get_concentration()
  {
    return concentrations.eigen_data;
  }

  [[nodiscard]] ColMajorKokkosScalarMatrix<double>
  ScalarSimulation::get_device_concentration() const
  {
    return concentrations.compute;
  }

  std::size_t ScalarSimulation::n_col() const noexcept
  {
    return n_c;
  }

  std::size_t ScalarSimulation::n_row() const noexcept
  {
    return n_r;
  }

  void ScalarSimulation::reduce_contribs(std::span<const double> data)
  {
    assert(data.size() == (n_c * n_r));
    using eigen_type = decltype(sources)::EigenMatrix;
    sources.eigen_data.noalias() += Eigen::Map<eigen_type>(
        const_cast<double*>(data.data()), EIGEN_INDEX(n_r), EIGEN_INDEX(n_c));
  }

  void ScalarSimulation::performStepGL(double d_t,
                                       const ColMajorMatrixtype<double>& mtr,
                                       MassTransfer::Sign sign)
  {
    PROFILE_SECTION("performStep_gl")
#define c concentrations.eigen_data

    total_mass =
        total_mass + d_t * (c * m_transition - c * sink + sources.eigen_data +
                            static_cast<float>(sign) * mtr);
    c = total_mass * volumes_inverse;

    // Make accessible new computed concentration to ComputeSpace
    concentrations.update_host_to_compute();
  }

  void ScalarSimulation::performStep(double d_t)
  {
    PROFILE_SECTION("performStep_l")
#define c concentrations.eigen_data

    total_mass =
        total_mass + d_t * (c * m_transition - c * sink + sources.eigen_data);
    c = total_mass * volumes_inverse;

    // Make accessible new computed concentration to ComputeSpace
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

  void ScalarSimulation::set_mass()
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
