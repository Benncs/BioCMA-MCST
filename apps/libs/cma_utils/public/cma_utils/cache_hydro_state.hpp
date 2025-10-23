#ifndef __CACHE_HYDRO_STATE__
#define __CACHE_HYDRO_STATE__

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

#include <common/common.hpp>
#include <vector>

template <typename ExecSpace>
using DiagonalView = Kokkos::View<double*,
                                  Kokkos::LayoutLeft,
                                  ExecSpace,
                                  Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

template <typename Space>
using CumulativeProbabilityView =
    Kokkos::View<double**,
                 Kokkos::LayoutRight,
                 Space,
                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

using FlowMatrixType = Eigen::SparseMatrix<double>;

/**
  @brief Namespace to handle algorithms and structures related to reading
  compartment mesh
  @see Simulation::KernelInline::MoveFunctor

  @note See documentation page for full explanation
*/
namespace CmaUtils
{
  class ProxyPreCalculatedHydroState;

  /**
   * @brief Structure to store hydrodynamic properties from a compartment mesh.
   *
   * This class is populated during the reading and processing of a
   * Compartment-Mesh case and is later accessed during transient simulations.
   * The primary advantage is that it allows reusing precomputed states,
   * eliminating the need to recalculate them from the mesh. Access and
   * modification are performed via the ProxyCalculatedHydroState class, which
   * ensures controlled interaction with the private members of this class.
   */
  class PreCalculatedHydroState
  {
  public:
    friend class ProxyCalculatedHydroState;

    /**
     * @brief Default constructor.
     */
    PreCalculatedHydroState();

    /**
     * @brief Default destructor.
     */
    ~PreCalculatedHydroState() = default;

    /**
     * @brief Constructor initializing with a given flow matrix.
     * @param _tm Reference to a FlowMatrixType object.
     */
    explicit PreCalculatedHydroState(const FlowMatrixType& _tm);

    /**
     * @brief Default move constructor.
     */
    PreCalculatedHydroState& operator=(PreCalculatedHydroState&& rhs) = default;
    /**
     * @brief Deleted copy constructor.
     */
    PreCalculatedHydroState&
    operator=(const PreCalculatedHydroState& rhs) = delete;
    /**
     * @brief Default move assigment.
     */
    PreCalculatedHydroState(PreCalculatedHydroState&& other) = default;
    /**
     * @brief Deleted copy assigment.
     */
    PreCalculatedHydroState(const PreCalculatedHydroState& other) = delete;
    /**
     * @brief Get the transition matrix.
     * @return Constant reference to FlowMatrixType.
     */
    [[nodiscard]] const FlowMatrixType& get_transition() const;

    /**
     * @brief Get the kernel diagonal view.
     * @return DiagonalView object for kernel computation.
     */
    [[nodiscard]] DiagonalView<ComputeSpace> get_kernel_diagonal() const;

    FlowMatrixType transition_matrix; ///< Matrix representing flow transitions.

    Eigen::Matrix<double, -1, -1, Eigen::RowMajor>
        cumulative_probability; ///< Cumulative probability matrix.

    std::vector<double> inverse_volume; ///< Inverse of compartment volumes.
    std::vector<double> volume;         ///< Volumes of compartments.
    DiagonalView<ComputeSpace>
        diagonal_compute; ///< Diagonal view for compute operations.

    // CumulativeProbabilityView<ComputeSpace>
    //     compute_cumulative_probability; ///< View for cumulative probability
    //     computation.

  private:
  };

  inline DiagonalView<ComputeSpace>
  PreCalculatedHydroState::get_kernel_diagonal() const
  {
    return diagonal_compute;
  }

} // namespace CmaUtils

#endif //__CACHE_HYDRO_STATE__