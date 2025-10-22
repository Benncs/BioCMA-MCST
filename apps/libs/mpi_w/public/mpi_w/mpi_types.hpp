#ifndef __MPI_TYPES_HPP__
#define __MPI_TYPES_HPP__

#include <mpi.h>
#include <mpi_w/message_t.hpp>
#include <type_traits>

namespace WrapMPI
{

  template <typename T> constexpr MPI_Datatype get_type() noexcept
  {
    MPI_Datatype datatype{};

    using _type = std::remove_const_t<std::remove_reference_t<T>>;

    if constexpr (std::is_same_v<_type, size_t>)
    {
      datatype = MPI_UNSIGNED_LONG;
    }
    else if constexpr (std::is_same_v<_type, double>)
    {
      datatype = MPI_DOUBLE;
    }
    else if constexpr (std::is_same_v<_type, int>)
    {
      datatype = MPI_INT;
    }
    else if constexpr (std::is_same_v<_type, bool>)
    {
      datatype = MPI_CXX_BOOL;
    }
    else if constexpr (std::is_same_v<_type, char> ||
                       std::is_same_v<_type, WrapMPI::SIGNALS>)
    {
      datatype = MPI_BYTE;
    }
    else
    {
      []<bool flag = false>() { static_assert(flag, "no match"); }();
    }

    return datatype;
  }
} // namespace WrapMPI

#endif //__MPI_TYPES_HPP__