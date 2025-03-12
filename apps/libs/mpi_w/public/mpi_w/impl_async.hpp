#ifndef __IMPL_ASYNC_MPI_HPP__
#define __IMPL_ASYNC_MPI_HPP__

#include "common/common.hpp"
#include <mpi_w/message_t.hpp>
#include <mpi_w/mpi_types.hpp>

#include <common/execinfo.hpp>
#include <cstddef>
#include <limits>
#include <math.h>
#include <mpi.h>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

namespace WrapMPI::Async
{

  inline MPI_Status wait(MPI_Request& request)
  {
    MPI_Status status;
    MPI_Wait(&request, &status);
    return status;
  }

  template <POD_t DataType>
  static int _send_unsafe(
      MPI_Request& request, DataType* buf, size_t buf_size, size_t dest, size_t tag) noexcept
  {
    return MPI_Isend(buf, buf_size, get_type<DataType>(), dest, tag, MPI_COMM_WORLD, &request);
  }
  template <POD_t DataType> int send(MPI_Request& request, DataType data, size_t dest, size_t tag)
  {
    return _send_unsafe<DataType>(request, &data, 1, dest, tag);
  }

  template <POD_t DataType>
  int send_v(MPI_Request& request,
             std::span<const DataType> data,
             size_t dest,
             size_t tag,
             bool send_size) noexcept
  {
    int send_status = MPI_SUCCESS;

    if (send_size)
    {
      send_status = send<size_t>(request, data.size(), dest, tag);
    }

    if (send_status == MPI_SUCCESS)
    {
      send_status = _send_unsafe(request, data.data(), data.size(), dest, tag);
    }

    return send_status;
  }

  template <POD_t DataType>
  int recv_span(MPI_Request& request, std::span<DataType> buf, size_t src, size_t tag) noexcept
  {
    return MPI_Irecv(
        buf.data(), buf.size(), get_type<DataType>(), src, tag, MPI_COMM_WORLD, &request);
  }

} // namespace WrapMPI::Async

#endif