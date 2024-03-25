#ifndef __MPI_TYPES_HPP__
#define __MPI_TYPES_HPP__

#include <mpi.h>

namespace MPI_W
{
  template <typename T> struct MPI_TYPES
  {
    static_assert(std::is_same_v<T, int> || std::is_same_v<T, double>|| std::is_same_v<T, size_t>,
                  "Unsupported type ");
    static MPI_Datatype value;
  };

 

} // namespace MPI_W




#endif //__MPI_TYPES_HPP__