#ifndef __IMPL_MPI_OP_HPP__
#define __IMPL_MPI_OP_HPP__

#include "common/common.hpp"
#include <mpi_w/message_t.hpp>
#include <mpi_w/mpi_types.hpp>

#include <common/execinfo.hpp>
#include <cstddef>
#include <limits>
#include <mpi.h>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

/**
  @brief Namespace to correclty wrap MPI C API for modern C++
*/
namespace WrapMPI
{

  template <typename T>
  concept POD_t =
      std::is_standard_layout_v<T> && std::is_trivially_copyable_v<T> &&
      std::is_trivially_destructible_v<T> && std::is_trivially_default_constructible_v<T>;

  // SENDING

  /**
   * @brief Sends raw data to a destination in an unsafe manner.
   *
   * This function sends raw data of type `DataType` to the specified destination.
   * It assumes the provided buffer and its size are valid
   *
   * @tparam DataType A type satisfying the `POD` concept.
   * @param buf Pointer to the buffer containing data to send.
   * @param buf_size The size of the buffer in bytes.
   * @param dest The destination identifier for the data.
   * @param tag Optional tag to identify the message (default is 0).
   * @return An integer indicating success or failure of the operation.
   *
   * @note Use this function with caution as it performs no validation on the input.
   */
  template <POD_t DataType>
  [[nodiscard]] static int
  _send_unsafe(DataType* buf, size_t buf_size, size_t dest, size_t tag = 0) noexcept;

  /**
   * @brief Sends a single instance of data to a destination.
   *
   * This function sends a single object of type `DataType` to the specified destination.
   *
   * @tparam DataType A type satisfying the `POD` concept.
   * @param data The data object to send.
   * @param dest The destination identifier for the data.
   * @param tag Optional tag to identify the message (default is 0).
   * @return An integer indicating success or failure of the operation.
   */
  template <POD_t DataType> int send(DataType data, size_t dest, size_t tag = 0);

  /**
   * @brief Sends a vector of data to a destination.
   *
   * This function sends a collection of objects of type `DataType` to the specified destination.
   * Optionally, the size of the collection can also be transmitted.
   *
   * @tparam DataType A type satisfying the `POD` concept.
   * @param data A `std::span` containing the data to send.
   * @param dest The destination identifier for the data.
   * @param tag Optional tag to identify the message (default is 0).
   * @param send_size If `true`, the size of the data will be sent before the data itself (default is `true`).
   * @return An integer indicating success or failure of the operation.
   */
  template <POD_t DataType>
  int send_v(std::span<const DataType> data,
             size_t dest,
             size_t tag = 0,
             bool send_size = true) noexcept;

  // RECEIVE

  /**
   * @brief Receives a single data item from a source.
   *
   * This function attempts to receive a single object of type `DataType` from the specified source.
   * If the receive operation fails, it returns an empty `std::optional`.
   *
   * @tparam DataType A type satisfying the `POD` concept.
   * @param src The source identifier from which to receive data.
   * @param status Optional pointer to an `MPI_Status` structure to retrieve message details
   * (default is `nullptr`).
   * @param tag Optional tag to filter messages (default is 0).
   * @return A `std::optional` containing the received data if successful, or empty if the operation
   * fails.
   */
  template <POD_t DataType>
  std::optional<DataType> recv(size_t src, MPI_Status* status = nullptr, size_t tag = 0) noexcept;

  /**
   * @brief Receives data into a span buffer.
   *
   * This function attempts to receive data into the provided span buffer. The buffer must
   * be pre-allocated and of sufficient size to hold the incoming data.
   *
   * @tparam DataType A type satisfying the `POD` concept.
   * @param buf A `std::span` representing the destination buffer for the data.
   * @param src The source identifier from which to receive data.
   * @param status Optional pointer to an `MPI_Status` structure to retrieve message details
   * (default is `nullptr`).
   * @param tag Optional tag to filter messages (default is 0).
   * @return An integer indicating the success or failure of the operation.
   */
  template <POD_t DataType>
  int recv_span(std::span<DataType> buf,
                size_t src,
                MPI_Status* status = nullptr,
                size_t tag = 0) noexcept;

  /**
   * @brief Attempts to receive a single data item from a source.
   *
   * This function attempts to receive a single object of type `DataType` from the specified source.
   * If the operation fails, it throws an exception.
   *
   * @tparam DataType A type satisfying the `POD` concept.
   * @param src The source identifier from which to receive data.
   * @param status Optional pointer to an `MPI_Status` structure to retrieve message details
   * (default is `nullptr`).
   * @param tag Optional tag to filter messages (default is 0).
   * @return The received data item.
   * @throws An exception if the receive operation fails.
   */
  template <POD_t DataType>
  DataType try_recv(size_t src, MPI_Status* status = nullptr, size_t tag = 0);

  /**
   * @brief Receives a vector of data items from a source.
   *
   * This function attempts to receive a vector of objects of type `T` from the specified source.
   * If the receive operation fails, it returns an empty `std::optional`.
   *
   * @tparam T The type of data to receive.
   * @param source The source identifier from which to receive data.
   * @param status Optional pointer to an `MPI_Status` structure to retrieve message details
   * (default is `nullptr`).
   * @param tag Optional tag to filter messages (default is 0).
   * @return A `std::optional` containing the received vector if successful, or empty if the
   * operation fails.
   */
  template <POD_t T>
  std::optional<std::vector<T>>
  recv_v(size_t source, MPI_Status* status = nullptr, size_t tag = 0) noexcept;

  /**
   * @brief Attempts to receive a vector of data items from a source.
   *
   * This function attempts to receive a vector of objects of type `T` from the specified source.
   * If the operation fails, it throws an exception.
   *
   * @tparam T The type of data to receive.
   * @param src The source identifier from which to receive data.
   * @param status Optional pointer to an `MPI_Status` structure to retrieve message details
   * (default is `nullptr`).
   * @param tag Optional tag to filter messages (default is 0).
   * @return A vector containing the received data items.
   * @throws An exception if the receive operation fails.
   */
  template <POD_t T>
  std::vector<T> try_recv_v(size_t src, MPI_Status* status = nullptr, size_t tag = 0);

  // BROADCASTING

  /**
   * @brief Broadcasts raw data to all processes in an unsafe manner.
   *
   * This function sends raw data of type `T` from the root process to all other processes.
   * It performs only basic safety checks
   *
   * @tparam T The type of the data to broadcast.
   * @param data Pointer to the data to be broadcasted.
   * @param _size The size of the data in bytes.
   * @param root The identifier of the root process that initiates the broadcast.
   * @return An integer indicating the success or failure of the operation.
   *
   * @note Use this function with caution
   */
  template <POD_t T> [[nodiscard]] int _broadcast_unsafe(T* data, size_t _size, size_t root);

  /**
   * @brief Broadcasts a single data item to all processes.
   *
   * This function sends a single object of type `DataType` from the root process to all other
   * processes.
   *
   * @tparam DataType A type satisfying the `POD` concept.
   * @param data Reference to the data to be broadcasted.
   * @param root The identifier of the root process that initiates the broadcast.
   * @return An integer indicating the success or failure of the operation.
   */
  template <POD_t DataType> int broadcast(DataType& data, size_t root) noexcept;

  /**
   * @brief Broadcasts data stored in a span to all processes.
   *
   * This function sends data stored in a `std::span` from the root process to all other processes.
   *
   * @tparam T The type of the data to broadcast.
   * @param data A `std::span` containing the data to be broadcasted.
   * @param root The identifier of the root process that initiates the broadcast.
   * @return An integer indicating the success or failure of the operation.
   */
  template <POD_t T> int broadcast_span(std::span<T> data, size_t root);

  /**
   * @brief Dispatches a task to the host with signal and arguments.
   *
   * This function sends execution information and arguments to the host for processing, triggered
   * by a specific signal.
   *
   * @tparam Args Variadic template types for the arguments.
   * @param info An `ExecInfo` object containing execution metadata.
   * @param sign A `SIGNALS` enumeration indicating the signal to trigger the task.
   * @param args Variadic arguments to pass to the host.
   */
  template <POD_t... Args> void host_dispatch(const ExecInfo& info, SIGNALS sign, Args&&... args);

  /**
   * @brief Gathers raw data from all processes in an unsafe manner.
   *
   * This function collects raw data of type `T` from all processes into a single vector on the root
   * process. It performs no safety checks and assumes the input pointers and sizes are valid.
   *
   * @tparam T The type of the data to gather.
   * @param src_data Pointer to the local data to be gathered.
   * @param size The size of the local data in bytes.
   * @param n_rank The number of processes in the communicator.
   * @param root The identifier of the root process that gathers the data (default is 0).
   * @return A vector containing the gathered data on the root process.
   *
   * @note Use this function with caution as it performs no validation on the input.
   */
  template <POD_t T>
  std::vector<T> _gather_unsafe(T* src_data, size_t size, size_t n_rank, size_t root = 0);

  template <POD_t T>
  int _gather_unsafe_to_buffer(T* dest, T* src_data, size_t size, size_t root = 0) noexcept;

  /**
   * @brief Gathers data from all processes.
   *
   * This function collects data stored in a `std::span` from all processes into a single vector on
   * the root process.
   *
   * @tparam T The type of the data to gather.
   * @param local_data A `std::span` containing the local data to be gathered.
   * @param n_rank The number of processes in the communicator.
   * @param root The identifier of the root process that gathers the data (default is 0).
   * @return A vector containing the gathered data on the root process.
   */
  template <POD_t T> std::vector<T> gather(std::span<T> local_data, size_t n_rank, size_t root = 0);

  template <POD_t T> void gather_span(std::span<T> dest, std::span<const T> local_data, size_t root = 0);

  template <POD_t T>
  std::vector<T> gather(std::span<const T> local_data, size_t n_rank, size_t root = 0);

  /**
   * @brief Gathers and reduces data to a single value.
   *
   * This function collects data of type `T` from all processes, performs a reduction operation
   * (e.g., sum), and returns the result to the root process.
   *
   * @tparam T A type satisfying the `NumberType` concept.
   * @param data The local data to be reduced and gathered.
   * @param root The identifier of the root process that gathers and reduces the data (default is
   * 0).
   * @return The reduced value on the root process.
   */
  template <NumberType T> T gather_reduce(T data, size_t root = 0);

  /**
   * @brief Gathers a vector of data from all processes.
   *
   * This function collects vectors of data from all processes into a single vector on the root
   * process.
   *
   * @tparam T The type of the data to gather.
   * @param local_data A vector containing the local data to be gathered.
   * @param n_rank The number of processes in the communicator.
   * @param root The identifier of the root process that gathers the data (default is 0).
   * @return A vector containing the gathered data on the root process.
   */
  template <POD_t T>
  std::vector<T> gather_v(const std::vector<T>& local_data, size_t n_rank, size_t root = 0);

  //**
  // IMPL
  //**
  template <POD_t DataType>
  static int _send_unsafe(DataType* buf, size_t buf_size, size_t dest, size_t tag) noexcept
  {
    return MPI_Send(buf, buf_size, get_type<DataType>(), dest, tag, MPI_COMM_WORLD);
  }

  template <POD_t DataType> DataType try_recv(size_t src, MPI_Status* status, size_t tag)
  {
    auto opt_data = WrapMPI::recv<DataType>(src, status, tag);
    if (!opt_data.has_value())
    {
      WrapMPI::critical_error();
      exit(-1); // critical_error should exit before reaching this statement
    }
    else
    {
      return opt_data.value();
    }
  }

  template <POD_t T>
  std::optional<std::vector<T>> recv_v(size_t source, MPI_Status* status, size_t tag) noexcept
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
    int recv_status = MPI_Recv(
        buf.data(), static_cast<int>(buf_size), datatype, source, tag, MPI_COMM_WORLD, status);

    if (recv_status != MPI_SUCCESS)
    {
      return std::nullopt; // Return early if MPI_Recv fails
    }

    return buf; // Return the received vector
  }

  template <POD_t DataType>
  int recv_span(std::span<DataType> buf, size_t src, MPI_Status* status, size_t tag) noexcept
  {
    return MPI_Recv(buf.data(), buf.size(), get_type<DataType>(), src, tag, MPI_COMM_WORLD, status);
  }

  template <POD_t T> std::vector<T> try_recv_v(size_t src, MPI_Status* status, size_t tag)
  {
    auto opt_data = WrapMPI::recv_v<T>(src, status, tag);
    if (!opt_data.has_value())
    {
      WrapMPI::critical_error();
    }
    else
    {
      return opt_data.value();
    }
  }

  template <POD_t DataType> int broadcast(DataType& data, size_t root) noexcept
  {
    return MPI_Bcast(&data, 1, get_type<DataType>(), static_cast<int>(root), MPI_COMM_WORLD);
  }

  // template <> int broadcast(size_t &data, size_t root)
  // {
  //   return MPI_Bcast(&data,
  //                    sizeof(size_t),
  //                    MPI_UNSIGNED_LONG,
  //                    static_cast<int>(root),
  //                    MPI_COMM_WORLD);
  // }

  template <POD_t T> int broadcast(std::vector<T>& data, size_t root, size_t current_rank)
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

    return MPI_Bcast(data.data(), data_size, get_type<T>(), static_cast<int>(root), MPI_COMM_WORLD);
  }

  template <POD_t T> int _broadcast_unsafe(T* data, size_t _size, size_t root)
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
    return MPI_Bcast(data, _size, get_type<T>(), static_cast<int>(root), MPI_COMM_WORLD);
  }

  template <POD_t T> int broadcast_span(std::span<T> data, size_t root)
  {
    return _broadcast_unsafe(data.data(), data.size(), root);
  }

  template <POD_t T>
  std::vector<T> _gather_unsafe(T* src_data, size_t size, size_t n_rank, size_t root)
  {
    int src_size = static_cast<int>(size);
    std::vector<T> total_data(src_size * n_rank);
    T* dest_data = total_data.data();
    // auto mpi_type = get_type<T>();

    // int gather_result = MPI_Gather(
    //     src_data, src_size, mpi_type, dest_data, src_size, mpi_type, root, MPI_COMM_WORLD);
    // if (gather_result != MPI_SUCCESS)
    // {
    //   throw std::runtime_error("MPI_Gather failed");
    // }
    if (_gather_unsafe_to_buffer(dest_data, src_data, size, root) != MPI_SUCCESS)
    {
      throw std::runtime_error("MPI_Gather failed");
    }

    return total_data;
  }

  template <POD_t T>
  int _gather_unsafe_to_buffer(T* const dest, T* src_data, size_t size, size_t root) noexcept
  {
    int src_size = static_cast<int>(size);
    auto mpi_type = get_type<T>();

    return MPI_Gather(src_data, src_size, mpi_type, dest, src_size, mpi_type, root, MPI_COMM_WORLD);
  }

  template <POD_t T> void gather_span(std::span<T> dest, std::span<const T> local_data, size_t root)
  {
    T* dest_data = dest.data();
#ifndef NDEBUG
    int size{};
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    assert(dest.size() == local_data.size() * size);
#endif
    _gather_unsafe_to_buffer(dest_data, const_cast<T*>(local_data.data()), local_data.size(), root);
  }

  template <POD_t T> std::vector<T> gather(std::span<T> local_data, size_t n_rank, size_t root)
  {
    return _gather_unsafe(local_data.data(), local_data.size(), n_rank, root);
  }
  // FIXME
  template <POD_t T>
  std::vector<T> gather(std::span<const T> local_data, size_t n_rank, size_t root)
  {
    return _gather_unsafe(const_cast<T*>(local_data.data()), local_data.size(), n_rank, root);
  }

  template <NumberType T> T gather_reduce(T data, size_t root)
  {

    T result{};

    MPI_Reduce(&data, &result, 1, WrapMPI::get_type<T>(), MPI_SUM, root, MPI_COMM_WORLD);

    return result;
  }

  template <POD_t T>
  std::vector<T> gather_v(const std::vector<T>& local_data, size_t n_rank, size_t root)
  {
    return _gather_unsafe(local_data.data(), local_data.size(), n_rank, root);
  }

  template <POD_t... Args> void host_dispatch(const ExecInfo& info, SIGNALS sign, Args&&... args)
  {

    for (int j = 1; j < static_cast<int>(info.n_rank); ++j)
    {
      if (sign != WrapMPI::SIGNALS::NOP)
      {
        MPI_Send(&sign, sizeof(sign), MPI_CHAR, j, 0, MPI_COMM_WORLD);
      }

      (
          [&]<POD_t T>(T& arg)
          {
            size_t s = 0;
            void* buf = nullptr;
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

  template <POD_t DataType> int send(DataType data, size_t dest, size_t tag)
  {
    return _send_unsafe<DataType>(&data, 1, dest, tag);
  }

  template <POD_t DataType>
  int send_v(std::span<const DataType> data, size_t dest, size_t tag, bool send_size) noexcept
  {
    int send_status = MPI_SUCCESS;

    if (send_size)
    {
      send_status = send<size_t>(data.size(), dest, tag);
    }

    if (send_status == MPI_SUCCESS)
    {
      send_status = _send_unsafe(data.data(), data.size(), dest, tag);
    }

    return send_status;
  }

  template <POD_t DataType>
  std::optional<DataType> recv(size_t src, MPI_Status* status, size_t tag) noexcept
  {
    DataType buf;

    int recv_status = MPI_Recv(&buf, sizeof(DataType), MPI_BYTE, src, tag, MPI_COMM_WORLD, status);
    if (recv_status != MPI_SUCCESS)
    {
      return std::nullopt;
    }
    return buf;
  }

} // namespace WrapMPI

#endif //__IMPL_MPI_OP_HPP__