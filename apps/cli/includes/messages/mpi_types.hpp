#ifndef __MPI_TYPES_HPP__
#define __MPI_TYPES_HPP__

#include <mpi.h>

namespace MPI_W
{
  template <typename T> struct MPI_TYPES
  {
    static_assert(std::is_same_v<T, int> || std::is_same_v<T, double> ||
                      std::is_same_v<T, size_t>,
                  "Unsupported type ");
    static MPI_Datatype value;
  };

  template <typename T> constexpr MPI_Datatype get_type()
  {
    MPI_Datatype datatype = MPI_BYTE;
    if constexpr (std::is_same_v<T, size_t>)
    {
      datatype = MPI_UNSIGNED_LONG;
    }
    else if constexpr (std::is_same_v<T, double>)
    {
      datatype = MPI_DOUBLE;
    }
    else if constexpr (std::is_same_v<T, int>)
    {
      datatype = MPI_INT;
    }
    else if constexpr (std::is_same_v<T, char>)
    {
      datatype = MPI_BYTE;
    }
    else
    {

      throw std::runtime_error("Error");
    }

    return datatype;
  }
} // namespace MPI_W

#endif //__MPI_TYPES_HPP__