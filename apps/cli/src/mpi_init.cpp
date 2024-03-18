#include <messages/init.hpp>
#include <mpi.h>

ExecInfo init_mpi(int argc, char **argv)
{
  ExecInfo info;

  int rank, size;

  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Get the rank of the current process
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Get the total number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  info.current_rank = static_cast<size_t>(rank);
  info.n_rank = static_cast<size_t>(size);

  return info;
}