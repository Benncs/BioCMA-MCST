#include "mpi_w/impl_op.hpp"
#include <cstddef>
#include <mpi.h>
#include <mpi_w/iteration_payload.hpp>
#include <mpi_w/wrap_mpi.hpp>

namespace WrapMPI
{

  enum TagExchange
  {
    volumes = 0,
    flows,
    proba,
    neighbors
  };

  IterationPayload::IterationPayload(const size_t volumes)
      : liquid_volumes(volumes), liquid_out_flows(volumes)
  {
  }

  [[nodiscard]] bool HostIterationPayload::sendAll(std::size_t n_rank) noexcept
  {
    PROFILE_SECTION("host:All")
    bool flag = false;
    if (n_rank == 1)
    {
      return true;
    }

    for (size_t __macro_j = 1; __macro_j < n_rank; ++__macro_j)
    {
      auto _ = WrapMPI::send(WrapMPI::SIGNALS::HydroUpdate,
                             __macro_j); // Need to be synchro
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

    int rc1 = WrapMPI::Async::send_v(
        requests[0], liquid_volumes, rank, TagExchange::volumes, false);

    int rc2 = WrapMPI::Async::send_v(
        requests[1], liquid_out_flows, rank, TagExchange::flows, false);

    int rc3 = WrapMPI::Async::send_v(
        requests[2], proba_leaving_flat, rank, TagExchange::proba, true);

    int rc4 = WrapMPI::Async::send_v(
        requests[3], liquid_neighbors_flat, rank, TagExchange::neighbors, true);

    return rc1 == MPI_SUCCESS && rc2 == MPI_SUCCESS && rc3 == MPI_SUCCESS &&
           rc4 == MPI_SUCCESS;
  }

  bool IterationPayload::recv(const size_t source, MPI_Status* status) noexcept
  {

    int rc1 = WrapMPI::recv_span<double>(
        liquid_volumes, source, status, TagExchange::volumes);
    int rc2 = WrapMPI::recv_span<double>(
        liquid_out_flows, source, status, TagExchange::flows);

    auto opt_proba =
        WrapMPI::recv_v<double>(source, status, TagExchange::proba);
    if (!opt_proba.has_value())
    {
      return false;
    }
    proba_leaving_flat = opt_proba.value();

    auto opt = WrapMPI::recv_v<size_t>(source, status, TagExchange::neighbors);
    if (!opt.has_value())
    {
      return false;
    }

    liquid_neighbors_flat = opt.value();

    return rc1 == MPI_SUCCESS && rc2 == MPI_SUCCESS;
  }

  void HostIterationPayload::fill(const CmaUtils::IterationStatePtrType& state)
  {
    auto liq = state->get_liquid();
    liquid_out_flows = liq->out_flows();
    liquid_volumes = liq->volume();
    liquid_neighbors_flat = state->flat_neighobrs();
    proba_leaving_flat = state->flat_probability_leaving();
  }

} // namespace WrapMPI