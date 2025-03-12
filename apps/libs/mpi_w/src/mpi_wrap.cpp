#include <cstdlib>
#include <mpi.h>
#include <mpi_w/message_t.hpp>
#include <omp.h>
#include <mpi_w/impl_async.hpp>

namespace WrapMPI
{
  bool is_initialized() noexcept
  {
    int initialized{};

    MPI_Initialized(&initialized);

    return initialized != 0;
  }

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

  namespace Async
  {
  
  } // namespace Async
} // namespace WrapMPI
