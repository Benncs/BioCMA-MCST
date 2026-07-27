#ifdef FIX_SEED
#  undef DFIX_SEED
#endif

#include <core/case_data.hpp>
#include <core/init_rng_seed.hpp>
#include <iostream>
#ifndef NO_MPI
#  include <mpi_w/wrap_mpi.hpp>
#endif
int
main()
{

  const auto exec = Core::runtime_init(0, nullptr, std::cout);
  const auto seed = Core::get_rng_seed(exec);

#ifndef NO_MPI
  if (exec.current_rank == 0)
  {
    for (int i = 1; i < exec.n_rank; ++i)
    {
      const auto worker_seed = WrapMPI::try_recv<uint64_t>(i, nullptr, i);
      if (exec.current_rank == 0)
      {
        std::cout << "Seed " << i << " " << worker_seed << std::endl;
      }
      assert(worker_seed != seed);
    }
  }
  else
  {
    WrapMPI::send(seed, 0, exec.current_rank);
  }
#endif

  if (exec.current_rank == 0)
  {
    std::cout << "Seeds OK" << std::endl;
  }
  return 0;
}
