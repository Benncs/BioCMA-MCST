#include "messages/impl_op.hpp"
#include "messages/message_t.hpp"
#include <messages/iteration_payload.hpp>
#include <mpi.h>

namespace MPI_W
{

  IterationPayload::IterationPayload(size_t size_flows, size_t volumes)
  {
    this->liquid_flows.resize(size_flows);
    this->liquid_volumes.resize(volumes);
    this->gas_volumes.resize(volumes);
  }

  void HostIterationPayload::send(size_t rank) const
  {
    int rc1 = MPI_W::send_v<double>(liquid_flows, rank, 0, false);

    int rc2 = MPI_W::send_v<double>(liquid_volumes, rank, 1, false);

    int rc3 = MPI_W::send_v<double>(gas_volumes, rank, 2, false);

    if (rc1 != MPI_SUCCESS || rc2 != MPI_SUCCESS || rc3 != MPI_SUCCESS)
    {
      MPI_W::critical_error();
    }
  }

  void IterationPayload::recv(size_t source, MPI_Status *status)
  {

    int rc1 = MPI_W::recv_span<double>(liquid_flows, source, status, 0);
    int rc2 = MPI_W::recv_span<double>(liquid_volumes, source, status, 1);
    int rc3 = MPI_W::recv_span<double>(gas_volumes, source, status, 2);

    if (rc1 != MPI_SUCCESS || rc2 != MPI_SUCCESS || rc3 != MPI_SUCCESS)
    {
      MPI_W::critical_error();
    }
  }

} // namespace MPI_W