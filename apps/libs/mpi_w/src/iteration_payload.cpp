#include <cma_read/neighbors.hpp>
#include <cstddef>
#include <mpi.h>
#include <mpi_w/iteration_payload.hpp>
#include <mpi_w/wrap_mpi.hpp>

namespace WrapMPI
{

  IterationPayload::IterationPayload(const size_t size_flows, const size_t volumes)
  {
    this->liquid_flows.resize(size_flows);
    this->liquid_volumes.resize(volumes);
    this->gas_volumes.resize(volumes);
  }

  [[nodiscard]] bool HostIterationPayload::sendAll(std::size_t n_rank) noexcept
  {
    PROFILE_SECTION("host:All")
    bool flag = false;
    if(n_rank==1){return true;}
    // MPI_Request req_flag;
    // for (size_t __macro_j = 1; __macro_j < n_rank; ++__macro_j)
    // {
    //   auto _ = WrapMPI::send(WrapMPI::SIGNALS::RUN, __macro_j,tag);
    // }
    // MPI_Wait(&req_flag, MPI_STATUS_IGNORE); // Wait for each send to finish

    for (size_t __macro_j = 1; __macro_j < n_rank; ++__macro_j)
    {
      auto _ = WrapMPI::send(WrapMPI::SIGNALS::RUN, __macro_j);
      flag = this->send(__macro_j);
    }

    for (auto& req : requests)
    {
      MPI_Wait(&req, MPI_STATUS_IGNORE); // Wait for each send to finish
    }

    return flag;
  }

  bool HostIterationPayload::send(const size_t rank) noexcept
  {
    PROFILE_SECTION("host:to_node")
    int rc1 = WrapMPI::Async::send_v<double>(requests[0], liquid_flows, rank, 0, false);

    int rc2 = WrapMPI::Async::send_v<double>(requests[1], liquid_volumes, rank, 1, false);

    int rc3 = WrapMPI::Async::send_v<double>(requests[2], gas_volumes, rank, 2, false);

    int rc4 = WrapMPI::Async::send_v(requests[3], neighbors.data(), rank, 3, true);

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

  void HostIterationPayload::fill(const CmaRead::ReactorState& current_reactor_state)
  {
    liquid_flows = current_reactor_state.liquid_flow.getViewFlows().data();
    liquid_volumes = current_reactor_state.liquidVolume;
    gas_volumes = current_reactor_state.gasVolume;
    neighbors = current_reactor_state.liquid_flow.getViewNeighors().to_const();
  }

} // namespace WrapMPI