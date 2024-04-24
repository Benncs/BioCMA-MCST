#include <common/common.hpp>
#include <cstdlib>
#include <messages/message_t.hpp>
#include <messages/mpi_types.hpp>
#include <messages/wrap_mpi.hpp>
#include <mpi.h>
#include <omp.h>

namespace MPI_W
{

  void finalize() noexcept
  {
    MPI_Finalize();
  }

  int critical_error() noexcept
  {
    return MPI_Abort(MPI_COMM_WORLD, 0);
  }

  void barrier() noexcept
  {
    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
    {
      critical_error();
    }
  }
} // namespace MPI_W