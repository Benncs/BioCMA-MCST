#include <cassert>
#include <mpi_w/iteration_payload.hpp>
#include <iostream>






int main(int argc, char** argv)
{
  
  int rank = 0;
  int size = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<double> flows = {1.1, 2.2, 3.3};
  std::vector<double> volumes = {4.4, 5.5, 6.6};
  std::vector<double> gases = {7.7, 8.8, 9.9};

  std::vector<size_t> raw_neighbors = {1, 2, 3, 1, 2, 3, 1, 2, 3};

  CmaRead::L2DView<const size_t> n_view(raw_neighbors, 3);
  if (size != 1)
  {

    MPI_Status status;
    if (rank == 0)
    { // Sender

      WrapMPI::HostIterationPayload host_payload;
      host_payload.liquid_flows= flows;
      host_payload.liquid_volumes=volumes;
      host_payload.gas_volumes=gases;
      host_payload.neighbors = n_view;
      auto _ = host_payload.sendAll(size);
    }
    else if (rank == 1)
    {
      WrapMPI::IterationPayload payload(3, 3);
      payload.recv(0, &status);

      assert(payload.liquid_flows == flows);
      assert(payload.liquid_volumes == volumes);
      assert(payload.gas_volumes == gases);
      assert(raw_neighbors == payload.raw_neighbors);
    }
  }
  std::cout << "OK" << std::endl;
  MPI_Finalize();

  return 0;
}
