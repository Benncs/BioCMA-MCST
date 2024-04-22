#ifndef __MPI_TYPES_HPP__
#define __MPI_TYPES_HPP__

#include "messages/message_t.hpp"
#include <mpi.h>
#include <stdexcept>
#include <type_traits>

namespace MPI_W
{

  template <typename T> consteval MPI_Datatype get_type()
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
    else if constexpr (std::is_same_v<_type, char> ||
                       std::is_same_v<_type, MPI_W::SIGNALS>)
    {
      datatype = MPI_BYTE;
    }
    else
    {
      []<bool flag = false>()
      {
        static_assert(flag, "no match");
      }
      ();
    }

    return datatype;
  }
} // namespace MPI_W

#endif //__MPI_TYPES_HPP__