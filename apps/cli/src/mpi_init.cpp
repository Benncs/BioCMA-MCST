
#include <cstdlib>
#include <messages/wrap_mpi.hpp>
#include <mpi.h>
#include <omp.h>
#include <common/common.hpp>

namespace MPI_W
{
  

  ExecInfo init_mpi(int argc, char **argv)
  {
    ExecInfo info;

    int rank, size;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    info.current_rank = static_cast<size_t>(rank);
    info.n_rank = static_cast<size_t>(size);

    size_t omp_n_thread = omp_get_max_threads();
    if (omp_n_thread >= info.n_rank)
    {
      info.thread_per_process = omp_n_thread / info.n_rank;

      if (omp_n_thread % info.n_rank != 0 &&
          info.current_rank == info.n_rank - 1)
      {
        info.thread_per_process += omp_n_thread % info.n_rank;
      }
    }
    else
    {
      info.thread_per_process = 1;
    }
    omp_set_num_threads(static_cast<int>(info.thread_per_process));

    // std::cout << "NUM thread per process " << info.thread_per_process
    //           << std::endl;
    std::atexit(MPI_W::finalize);
    MPI_W::is_mpi_init = true;
    return info;
  }

} // namespace MPI_W
