
#include <rt_init.hpp>

#include <Eigen/Core>
#include <common/common.hpp>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mpi_w/wrap_mpi.hpp>
#include <sstream>

#ifndef USE_PYTHON_MODULE
#  include <omp.h>

#else
#  define omp_get_max_threads() 1 // 1 thread
#  define omp_get_num_procs() 1
#  define omp_set_num_threads(__arg__)
#endif

std::string env_file_path()
{
  return ".bmc_info";
}

size_t generate_run_id();

void set_openmp_threads(int rank,
                        int size,
                        ExecInfo &info,
                        SimulationParameters &params)
{
  // Casting rank and size to size_t
  info.current_rank = static_cast<size_t>(rank);
  info.n_rank = static_cast<size_t>(size);

  // Determining the number of OpenMP threads
  size_t omp_n_thread = (params.n_threads > 0)
                            ? static_cast<size_t>(params.n_threads)
                            : static_cast<size_t>(omp_get_max_threads());

  int num_core_per_node = omp_get_num_procs();

  size_t threads_per_process = 1;

  if (omp_n_thread >= info.n_rank)
  {
    threads_per_process = omp_n_thread / info.n_rank;
    if (omp_n_thread % info.n_rank != 0 && info.current_rank == info.n_rank - 1)
    {
      threads_per_process += omp_n_thread % info.n_rank;
    }
  }

  if (threads_per_process > static_cast<size_t>(num_core_per_node))
  {
    threads_per_process = static_cast<size_t>(num_core_per_node);
  }
  info.thread_per_process = threads_per_process;

  omp_set_num_threads(static_cast<int>(info.thread_per_process));
}

ExecInfo runtime_init(int argc, char **argv, SimulationParameters &params)
{
  ExecInfo info{};

  int rank = 0;
  int size = 0;
  int mpi_thread_level{};
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &mpi_thread_level);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm_size(MPI_COMM_WORLD, &size);

  set_openmp_threads(rank, size, info, params);

  Eigen::setNbThreads(1); //FIXME: 4 threads is clearly enough for 500*500 matrix 
  

#ifdef USE_PYTHON_MODULE
  info.thread_per_process = 1; // Set one thread because of PYthon GIL
#endif
      // TODO FIXME : disable OMP feature
  std::cout << "Numberof thread per process " << info.thread_per_process
            << std::endl;
  std::atexit(MPI_W::finalize);
  info.run_id =
      static_cast<size_t>(time(nullptr) * size * info.thread_per_process);
  // MPI_W::is_mpi_init = true;
  return info;
}

void init_environment()
{
  const auto env_path = env_file_path();

  static const std::string cma_data_folder_path = "cma_data";
  static const std::string results_folder_path = "results";
  static bool environment_initialized = false;

  if (!environment_initialized)
  {
    if (!std::filesystem::exists(cma_data_folder_path))
    {
      std::filesystem::create_directory(cma_data_folder_path);
    }

    if (!std::filesystem::exists(results_folder_path))
    {
      std::filesystem::create_directory(results_folder_path);
    }

    if (!std::filesystem::exists(env_path))
    {
      std::cout << "Setting up environment ...";

      std::ofstream bmc_info_file(env_path);
      if (bmc_info_file.is_open())
      {
        std::filesystem::permissions(env_path,
                                     std::filesystem::perms::owner_write |
                                         std::filesystem::perms::owner_read |
                                         std::filesystem::perms::others_read,
                                     std::filesystem::perm_options::replace);
        append_date_time(bmc_info_file);
        bmc_info_file << std::endl;
        bmc_info_file.close();
        std::cout << "Done" << '\n';
      }
      else
      {
        std::cout << "Failed to create .bmc_info file" << '\n';
      }
    }

    environment_initialized = true;
  }

}

// void append_date_time(std::ofstream &fd)
// {
//   auto now =
//       std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
//   fd << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S");
// }

std::string sappend_date_time(std::string_view string)
{
  std::stringstream fd;
  fd<<string;
  auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  fd << std::put_time(std::localtime(&now), "%Y-%m-%d-%H:%M:%S");
  return fd.str();

}

void register_run(ExecInfo &exec, SimulationParameters &params)
{
  // Open the file in append mode
  std::ofstream env(env_file_path(), std::ios_base::app);
  if (env.is_open())
  {
    append_date_time(env);
    env<< "\t";
    env << exec;
    env << params;
    env << std::endl;
  }
  else
  {
    std::cerr << "Error: Unable to open file for writing\n";
  }
}