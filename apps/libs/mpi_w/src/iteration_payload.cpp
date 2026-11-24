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

  [[nodiscard]] bool
  HostIterationPayload::sendAll(std::size_t n_rank) noexcept
  {
    PROFILE_SECTION("host:sendAll")
    bool flag = false;
    if (n_rank == 1)
    {
      to_wait = false;
      return true;
    }
    // Hotfix: introduce pending vector to store all requests done by all rank
    // size of 'requests' array is not sufficient and all request object are
    // lost (then maybe unclosed) except the last one
    // TODO: make requests a vector resized at construction with n_rank and use
    // the correct index in HostIterationPayload::send
    for (size_t __macro_j = 1; __macro_j < n_rank; ++__macro_j)
    {
      auto _ = WrapMPI::send(WrapMPI::SIGNALS::HydroUpdate,
                             __macro_j); // Need to be synchro
      flag = this->send(__macro_j);
      for (auto& req : requests)
      {
        pending.emplace_back(req);
      }
    }
    to_wait = true;

    return flag;
  }

  void
  HostIterationPayload::wait() noexcept
  {
    if (to_wait)
    {
      for (auto& req : pending)
      {
        MPI_Wait(&req, MPI_STATUS_IGNORE); // Wait for each send to finish
      }
      pending.clear();
      to_wait = false;
    }
  }

  bool
  HostIterationPayload::send(const size_t rank) noexcept
  {
    PROFILE_SECTION("host:to_node")

    // Hotfix: sending local variable size via async crash in some condition
    // in send_v size is a local vriable that live out of scope then with async
    // this is UB (ARM 72 ranks)
    // TODO: stabilize the send_v function to correctly send the size of the
    // array.
    // Sending the size as a sync operation prevent crash but this is a bit
    // non-sens to have sync operation amoung asyncs
    int rc1 = WrapMPI::Async::send_v(
        requests[0], liquid_volumes, rank, TagExchange::volumes, false);

    int rc2 = WrapMPI::Async::send_v(
        requests[1], liquid_out_flows, rank, TagExchange::flows, false);

    auto _s3 = WrapMPI::send<size_t>(
        proba_leaving_flat.size(), rank, TagExchange::proba);
    int rc3 = WrapMPI::Async::send_v(
        requests[2], proba_leaving_flat, rank, TagExchange::proba, false);

    auto _s4 = WrapMPI::send<size_t>(
        liquid_neighbors_flat.size(), rank, TagExchange::neighbors);
    int rc4 = WrapMPI::Async::send_v(requests[3],
                                     liquid_neighbors_flat,
                                     rank,
                                     TagExchange::neighbors,
                                     false);

    return rc1 == MPI_SUCCESS && rc2 == MPI_SUCCESS && rc3 == MPI_SUCCESS
           && rc4 == MPI_SUCCESS;
  }

  bool
  IterationPayload::recv(const size_t source, MPI_Status* status) noexcept
  {

    int rc1 = WrapMPI::recv_span<double>(
        liquid_volumes, source, status, TagExchange::volumes);
    int rc2 = WrapMPI::recv_span<double>(
        liquid_out_flows, source, status, TagExchange::flows);

    auto opt_proba
        = WrapMPI::recv_v<double>(source, status, TagExchange::proba);
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

  void
  HostIterationPayload::fill(const CmaUtils::IterationStatePtrType& state)
  {
    auto liq = state->get_liquid();
    liquid_out_flows = liq->out_flows();
    liquid_volumes = liq->volume();
    liquid_neighbors_flat = state->flat_neighobrs();
    proba_leaving_flat = state->flat_probability_leaving();
  }

} // namespace WrapMPI
