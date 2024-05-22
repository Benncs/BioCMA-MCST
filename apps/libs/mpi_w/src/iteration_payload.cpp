
#include "cma_read/light_2d_view.hpp"
#include "cma_read/neighbors.hpp"
#include <cstddef>
#include <mpi_w/iteration_payload.hpp>
#include <mpi_w/message_t.hpp>
#include <mpi_w/impl_op.hpp>
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

    
    int rc4 = MPI_W::send_v(neigbors.data(), rank,3,true);


    if (rc1 != MPI_SUCCESS || rc2 != MPI_SUCCESS || rc3 != MPI_SUCCESS|| rc4 != MPI_SUCCESS)
    {
      MPI_W::critical_error();
    }
  }

  void IterationPayload::recv(size_t source, MPI_Status *status)
  {

    int rc1 = MPI_W::recv_span<double>(liquid_flows, source, status, 0);
    int rc2 = MPI_W::recv_span<double>(liquid_volumes, source, status, 1);
    int rc3 = MPI_W::recv_span<double>(gas_volumes, source, status, 2);

    auto opt  = MPI_W::recv_v<size_t>(source,status,3);
    if(!opt.has_value())
    {
      MPI_W::critical_error();
    }
    raw_neigbors = opt.value();
    auto n_col = raw_neigbors.size() / liquid_flows.size();
    neigbors = Neighbors::Neighbors_const_view_t(
        raw_neigbors, liquid_flows.size(), n_col, true);


    if (rc1 != MPI_SUCCESS || rc2 != MPI_SUCCESS || rc3 != MPI_SUCCESS )
    {
      MPI_W::critical_error();
    }
  }

} // namespace MPI_W