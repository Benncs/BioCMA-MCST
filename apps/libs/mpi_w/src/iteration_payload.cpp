#include <cma_read/neighbors.hpp>
#include <cstddef>
#include <mpi.h>
#include <mpi_w/impl_op.hpp>
#include <mpi_w/iteration_payload.hpp>
#include <mpi_w/message_t.hpp>

namespace WrapMPI
{

  IterationPayload::IterationPayload(const size_t size_flows, const size_t volumes)
  {
    this->liquid_flows.resize(size_flows);
    this->liquid_volumes.resize(volumes);
    this->gas_volumes.resize(volumes);
  }

  bool HostIterationPayload::send(const size_t rank) const noexcept
  {
    int rc1 = WrapMPI::send_v<double>(liquid_flows, rank, 0, false);

    int rc2 = WrapMPI::send_v<double>(liquid_volumes, rank, 1, false);

    int rc3 = WrapMPI::send_v<double>(gas_volumes, rank, 2, false);

    int rc4 = WrapMPI::send_v(neighbors.data(), rank, 3, true);

    return rc1 == MPI_SUCCESS && rc2 == MPI_SUCCESS && rc3 == MPI_SUCCESS && rc4 == MPI_SUCCESS;
  }

  bool IterationPayload::recv(const size_t source, MPI_Status* status) noexcept
  {

    int rc1 = WrapMPI::recv_span<double>(liquid_flows, source, status, 0);
    int rc2 = WrapMPI::recv_span<double>(liquid_volumes, source, status, 1);
    int rc3 = WrapMPI::recv_span<double>(gas_volumes, source, status, 2);

    auto opt = WrapMPI::recv_v<size_t>(source, status, 3);
    if (!opt.has_value())
    {
      return false;
    }
    raw_neighbors = opt.value();
    auto n_col = raw_neighbors.size() / gas_volumes.size();
    neighbors =
        CmaRead::Neighbors::Neighbors_const_view_t(raw_neighbors, gas_volumes.size(), n_col, true);
    return rc1 == MPI_SUCCESS && rc2 == MPI_SUCCESS && rc3 == MPI_SUCCESS;
  }

} // namespace WrapMPI