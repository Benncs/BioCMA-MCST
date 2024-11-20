#include <common/common.hpp>
#include <cstdlib>
#include <mpi_w/message_t.hpp>
#include <mpi_w/mpi_types.hpp>
#include <mpi_w/wrap_mpi.hpp>
#include <mpi.h>
#include <omp.h>

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
} // namespace WrapMPI
