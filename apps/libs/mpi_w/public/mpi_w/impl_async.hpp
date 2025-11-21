#ifndef __IMPL_ASYNC_MPI_HPP__
#define __IMPL_ASYNC_MPI_HPP__

#include <common/common.hpp>
#include <common/execinfo.hpp>
#include <cstddef>
#include <math.h>
#include <mpi.h>
#include <mpi_w/message_t.hpp>
#include <mpi_w/mpi_types.hpp>
#include <span>

namespace WrapMPI::Async
{
  namespace
  {
    template <POD_t DataType>
    int _send_unsafe(MPI_Request& request,
                     DataType* buf,
                     size_t buf_size,
                     size_t dest,
                     size_t tag) noexcept
    {
      return MPI_Isend(buf,
                       buf_size,
                       get_type<DataType>(),
                       dest,
                       tag,
                       MPI_COMM_WORLD,
                       &request);
    }

  } // namespace

  inline MPI_Status wait(MPI_Request& request)
  {
    PROFILE_SECTION("WrapMPI::wait");
    MPI_Status status;
    MPI_Wait(&request, &status);
    return status;
  }

  inline void wait(MPI_Request& request, MPI_Status* status)
  {
    MPI_Wait(&request, status);
  }

  template <POD_t DataType>
  int send(MPI_Request& request, DataType data, size_t dest, size_t tag)
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
  int recv_span(MPI_Request& request,
                std::span<DataType> buf,
                size_t src,
                size_t tag) noexcept
  {
    return MPI_Irecv(buf.data(),
                     buf.size(),
                     get_type<DataType>(),
                     src,
                     tag,
                     MPI_COMM_WORLD,
                     &request);
  }

  template <POD_t DataType>
  std::optional<DataType>
  recv(size_t src, MPI_Request& request, size_t tag) noexcept
  {
    DataType buf;

    int recv_status = MPI_Irecv(
        &buf, sizeof(DataType), MPI_BYTE, src, tag, MPI_COMM_WORLD, &request);
    if (recv_status != MPI_SUCCESS)
    {
      return std::nullopt;
    }
    return buf;
  }

  template <POD_t T>
  int _broadcast_unsafe(T* data,
                        size_t _size,
                        size_t root,
                        MPI_Request& request)
  {
    if (data == nullptr)
    {
      throw std::invalid_argument("Data pointer is null");
    }
    if (_size == 0 || _size > std::numeric_limits<size_t>::max())
    {
      throw std::invalid_argument("Error size");
    }

    int comm_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    if (root >= static_cast<size_t>(comm_size))
    {
      throw std::invalid_argument("Root process rank is out of range");
    }

    return MPI_Ibcast(data,
                      _size,
                      get_type<T>(),
                      static_cast<int>(root),
                      MPI_COMM_WORLD,
                      &request);
  }

  template <POD_t T>
  int broadcast_span(std::span<T> data, size_t root, MPI_Request& request)
  {
    return _broadcast_unsafe(data.data(), data.size(), root, request);
  }

} // namespace WrapMPI::Async

#endif