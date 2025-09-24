#ifndef __ITERATION_PAYLOAD_HPP__
#define __ITERATION_PAYLOAD_HPP__

#include <cma_read/reactorstate.hpp>
#include <cma_read/flow_iterator.hpp>
#include <cma_read/neighbors.hpp>
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
    /// Vector containing liquid flow values for the current iteration.
    std::vector<double> liquid_flows;

    /// Vector containing liquid volume values for the current iteration.
    std::vector<double> liquid_volumes;

    /// Vector containing gas volume values for the current iteration.
    std::vector<double> gas_volumes;

    /// Vector containing raw neighbor indices for the current iteration.
    std::vector<size_t> raw_neighbors;

    /// Constant view of neighbor data for the current iteration.
    CmaRead::Neighbors::Neighbors_const_view_t neighbors;

    /**
     * @brief Constructs an IterationPayload with specified sizes for flows and volumes.
     *
     * @param size_flows The number of elements in the liquid flows vector.
     * @param volumes The number of elements in the liquid and gas volumes vectors.
     * @note: Size are needed to alloc vector, this allow to use preallocated chunk when transfer
     */
    explicit IterationPayload(size_t size_flows, size_t volumes);

    /**
     * @brief Receives data for this payload from a specified source.
     *
     * This function receives liquid flows, liquid volumes, and gas volumes
     * from a given MPI source rank.
     *
     * @param source The MPI rank of the source process sending the data.
     * @param status Pointer to an MPI_Status object to store information about the receive
     * operation.
     *
     * @note This method uses MPI to perform the receive operation and assumes the MPI environment
     * is initialized.
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
    /// View liquid flow values to be sent for the current iteration.
    std::span<const double> liquid_flows;
    /// View liquid volmes values to be sent for the current iteration.
    std::span<const double> liquid_volumes;
    /// View gas flow values to be sent for the current iteration.
    std::span<const double> gas_volumes;
    /// View neighbors values to be sent for the current iteration.
    CmaRead::Neighbors::Neighbors_const_view_t neighbors;

    void fill(const CmaRead::ReactorState& current_reactor_state);

    [[nodiscard]] bool sendAll(std::size_t n_rank)  noexcept;
    private:
    /**
     * @brief Sends this payload to a specified MPI rank.
     *
     * This function sends the liquid flows, liquid volumes, and gas volumes
     * to a specified destination rank using MPI.
     *
     * @param rank The MPI rank of the destination process.
     *
     * @note This method uses MPI to perform the send operation and assumes the MPI environment is
     * initialized.
     */
    [[nodiscard]] bool send(size_t rank)  noexcept;
    
    static constexpr std::size_t n_vector_send = 4; 
    std::array<MPI_Request,n_vector_send> requests;
  };

} // namespace WrapMPI

#endif //__ITERATION_PAYLOAD_HPP__
