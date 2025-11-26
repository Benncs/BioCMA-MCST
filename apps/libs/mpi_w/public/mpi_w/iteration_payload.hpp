#ifndef __ITERATION_PAYLOAD_HPP__
#define __ITERATION_PAYLOAD_HPP__

#include <cma_utils/alias.hpp>
#include <mpi.h>
#include <span>
#include <vector>

namespace WrapMPI
{

  /**
   * @class IterationPayload
   * @brief Represents the payload of data exchanged during an iteration.
   *
   * This class encapsulates the data related to liquid flows, liquid volumes,
   * gas volumes, and neighbor information for a given iteration.
   */
  class IterationPayload
  {
  public:
    std::vector<double> liquid_volumes;             // n_c
    std::vector<std::size_t> liquid_neighbors_flat; // m*n_c
    std::vector<double> proba_leaving_flat;         // m*n_c
    std::vector<double> liquid_out_flows;           // n_c

    /*
      Estimated payload size with
      n_c = 500
      m =6
      size = 56000 bytes ≈ 54.69 KiB ≈ 0.055 MB > Eager limit MPI should
      triggers SM
    */

    /**
     * @brief Constructs an IterationPayload with specified sizes for flows and
     * volumes.
     *
     * @param size_flows The number of elements in the liquid flows vector.
     * @param volumes The number of elements in the liquid and gas volumes
     * vectors.
     * @note: Size are needed to alloc vector, this allow to use preallocated
     * chunk when transfer
     */
    explicit IterationPayload(size_t volumes);

    /**
     * @brief Receives data for this payload from a specified source.
     *
     * This function receives liquid flows, liquid volumes, and gas volumes
     * from a given MPI source rank.
     *
     * @param source The MPI rank of the source process sending the data.
     * @param status Pointer to an MPI_Status object to store information about
     * the receive operation.
     *
     * @note This method uses MPI to perform the receive operation and assumes
     * the MPI environment is initialized.
     */
    bool recv(size_t source, MPI_Status* status) noexcept;
  };

  /**
   * @class HostIterationPayload
   * @brief Represents the payload of data on the host side for an iteration.
   *
   * This class encapsulates the data to be sent
   */
  class HostIterationPayload
  {
  public:
    void fill(const CmaUtils::IterationStatePtrType& current_reactor_state);

    [[nodiscard]] bool sendAll(std::size_t n_rank) noexcept;

    void wait() noexcept;
#ifdef NDEBUG // FIXME Needed for unite test
  private:
#endif
    std::span<const double> liquid_volumes;
    std::span<const std::size_t> liquid_neighbors_flat;
    std::span<const double> proba_leaving_flat;
    std::span<const double> liquid_out_flows;

    bool to_wait;

    /**
     * @brief Sends this payload to a specified MPI rank.
     *
     * This function sends the liquid flows, liquid volumes, and gas volumes
     * to a specified destination rank using MPI.
     *
     * @param rank The MPI rank of the destination process.
     *
     * @note This method uses MPI to perform the send operation and assumes the
     * MPI environment is initialized.
     */
    [[nodiscard]] bool send(size_t rank) noexcept;

    static constexpr std::size_t n_vector_send = 4;
    std::array<MPI_Request, n_vector_send> requests;
  };

} // namespace WrapMPI

#endif //__ITERATION_PAYLOAD_HPP__
