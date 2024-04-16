#ifndef __IMPL_MPI_OP_HPP__
#define __IMPL_MPI_OP_HPP__

#include "messages/mpi_types.hpp"
#include "mpi.h"
#include <common/execinfo.hpp>
#include <cstddef>
#include <messages/message_t.hpp>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>
#include <limits>

namespace MPI_W
{

  template <typename T>
  concept POD = std::is_standard_layout<T>::value &&
                std::is_trivially_copyable<T>::value &&
                std::is_trivially_destructible<T>::value &&
                std::is_trivially_default_constructible<T>::value;

  template <POD DataType>
  int send(DataType &&data, size_t dest, size_t tag = MPI_ANY_TAG)
  {

    int status =
        MPI_Send(&data, 1, get_type<DataType>(), dest, tag, MPI_COMM_WORLD);
    return status;
  }

  template <POD DataType>
  int send_v(std::span<const DataType> data,
             size_t dest,
             size_t tag = MPI_ANY_TAG)
  {
    int status{};
    status = send<size_t>(data.size(), dest, tag);
    if (status == MPI_SUCCESS)
    {
      status = MPI_Send(data.data(),
                        data.size(),
                        get_type<DataType>(),
                        dest,
                        tag,
                        MPI_COMM_WORLD);
    }

    return status;
  }

  template <POD DataType>
  std::optional<DataType>
  recv(size_t src, MPI_Status *status = nullptr, size_t tag = MPI_ANY_TAG)
  {
    DataType buf;

    int recv_status = MPI_Recv(
        &buf, sizeof(DataType), MPI_BYTE, src, tag, MPI_COMM_WORLD, status);
    if (recv_status != MPI_SUCCESS)
    {
      return std::nullopt;
    }
    return buf;
  }

  template <POD DataType>
  DataType
  try_recv(size_t src, MPI_Status *status = nullptr, size_t tag = MPI_ANY_TAG)
  {
    auto opt_data = MPI_W::recv<DataType>(src, status, tag);
    if (!opt_data.has_value())
    {
      MPI_W::critical_error();
    }
    else
    {
      return opt_data.value();
    }
  }

  template <typename T>
  std::optional<std::vector<T>> recv_v(size_t source,
                                       MPI_Status *status = nullptr,
                                       size_t tag = MPI_ANY_TAG) noexcept
  {
    std::vector<T> buf;

    // Receive the size of the vector
    auto opt_size = recv<size_t>(source, status, tag);
    if (!opt_size.has_value())
    {
      return std::nullopt; // Return early if size reception fails
    }
    size_t buf_size = opt_size.value();

    // Resize the buffer
    buf.resize(buf_size);

    MPI_Datatype datatype = get_type<T>();

    // Receive the vector data
    int recv_status = MPI_Recv(buf.data(),
                               static_cast<int>(buf_size),
                               datatype,
                               source,
                               tag,
                               MPI_COMM_WORLD,
                               status);

    if (recv_status != MPI_SUCCESS)
    {
      return std::nullopt; // Return early if MPI_Recv fails
    }

    return buf; // Return the received vector
  }

  template <typename T>
  std::vector<T>
  try_recv_v(size_t src, MPI_Status *status = nullptr, size_t tag = MPI_ANY_TAG)
  {
    auto opt_data = MPI_W::recv_v<T>(src, status, tag);
    if (!opt_data.has_value())
    {
      MPI_W::critical_error();
    }
    else
    {
      return opt_data.value();
    }
  }

  template <POD DataType> int broadcast(DataType &data, size_t root)
  {
    return MPI_Bcast(&data,
                     sizeof(DataType),
                     MPI_BYTE,
                     static_cast<int>(root),
                     MPI_COMM_WORLD);
  }

  // template <> int broadcast(size_t &data, size_t root)
  // {
  //   return MPI_Bcast(&data,
  //                    sizeof(size_t),
  //                    MPI_UNSIGNED_LONG,
  //                    static_cast<int>(root),
  //                    MPI_COMM_WORLD);
  // }

  template <typename T>
  int broadcast(std::vector<T> &data, size_t root, size_t current_rank)
  {

    size_t data_size = 0;
    if (current_rank == root)
    {
      data_size = data.size();
    }
    broadcast(data_size, root);
    if (current_rank != root)
    {
      data.resize(data_size);
    }

    return MPI_Bcast(data.data(),
                     data_size,
                     MPI_TYPES<T>::value,
                     static_cast<int>(root),
                     MPI_COMM_WORLD);
  }

  template <typename T>
  int _broadcast_unsafe(T *data, size_t _size, size_t root)
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

    // Broadcast operation
    return MPI_Bcast(data,
                     _size,
                     MPI_TYPES<T>::value,
                     static_cast<int>(root),
                     MPI_COMM_WORLD);
  }

  template <typename T> int broadcast_span(std::span<T> data, size_t root)
  {
    return _broadcast_unsafe(data.data(), data.size(), root);
  }

  template <typename T>
  std::vector<T>
  gather(std::span<T> &&local_data, size_t n_rank, size_t root = 0)
  {
    int src_size = static_cast<int>(local_data.size());
    std::vector<T> total_data(src_size * n_rank);
    T *src_data = local_data.data();
    T *dest_data = total_data.data();
    auto mpi_type = get_type<T>();

    int gather_result = MPI_Gather(src_data,
                                   src_size,
                                   mpi_type,
                                   dest_data,
                                   src_size,
                                   mpi_type,
                                   root,
                                   MPI_COMM_WORLD);
    if (gather_result != MPI_SUCCESS)
    {
      throw std::runtime_error("MPI_Gather failed");
    }

    return total_data;
  }

  template <typename... Args>
  void host_dispatch(const ExecInfo &info, SIGNALS &&sign, Args &&...args)
  {

    for (int j = 1; j < static_cast<int>(info.n_rank); ++j)
    {
      if (sign != MPI_W::SIGNALS::NOP)
      {
        MPI_Send(&sign, sizeof(sign), MPI_CHAR, j, 0, MPI_COMM_WORLD);
      }

      (
          [&]<typename T>(T &&arg)
          {
            size_t s = 0;
            void *buf = nullptr;
            if constexpr (std::is_same_v<std::decay_t<T>, std::span<double>>)
            {
              s = arg.size();
              buf = arg.data();
              MPI_Send(&s, 1, MPI_UNSIGNED_LONG, j, 0, MPI_COMM_WORLD);
            }
            else
            {
              s = sizeof(T);
              buf = &arg;
            }

            MPI_Send(buf, static_cast<int>(s), MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
          }(std::forward<Args>(args)),
          ...);
    }
  }

} // namespace MPI_W

#endif //__IMPL_MPI_OP_HPP__