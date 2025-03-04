#include <Kokkos_Core.hpp>
#include <biocma_cst_config.hpp>
#include <cassert>
#include <chrono>
#include <common/execinfo.hpp>
#include <core/simulation_parameters.hpp>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <rt_init.hpp>
#include <sstream>
#include <string>
#include <string_view>
#include <Eigen/Core>
#include <api/api.hpp>

#ifndef NO_MPI
#  include <mpi_w/wrap_mpi.hpp>
#endif

#ifndef USE_PYTHON_MODULE
#  include <omp.h>
#else
#  define omp_get_max_threads() 1 // 1 thread
#  define omp_get_num_procs() 1
#  define omp_set_num_threads(__arg__)
#endif

std::string env_file_path() noexcept;
void init_environment();
size_t generate_run_id();

void check_results_file_name(Core::UserControlParameters& params);

ExecInfo runtime_init(int argc, char** argv, Core::UserControlParameters& params) noexcept
{
  init_environment();
  ExecInfo info{};

  int rank = 0;
  int size = 0;

  if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
  {
#ifndef NO_MPI

    // int mpi_thread_level{};
    // MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &mpi_thread_level);
    // int mpi_thread_level{};
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  }
  else
  {
    size = 1;
  }

  set_n_thread_current_rank(rank, size, info, params);
  Kokkos::initialize(Kokkos::InitializationSettings() 
                         .set_disable_warnings(false)
                         .set_num_threads(static_cast<int32_t>(info.thread_per_process))
                         .set_map_device_id_by("random"));
#ifndef NO_MPI
  WrapMPI::barrier();
#endif
  if (rank == 0)
  {
    if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
    {
        
      std::cout << "USING MPI" << std::endl;
    }
    Kokkos::print_configuration(std::cout);
    
  }


  Eigen::setNbThreads(6);
  // Eigen::setNbThreads(std::min(omp_get_num_procs(), 4));

  #ifdef USE_PYTHON_MODULE
    info.thread_per_process = 1; // Set one thread because of PYthon GIL
  #endif
  // TODO FIXME : disable OMP feature


  if constexpr (AutoGenerated::FlagCompileTime::use_mpi)
  {
#ifndef NO_MPI
    std::atexit(WrapMPI::finalize);
#endif
  }
  std::atexit(Api::finalise);

  const auto id_seed = static_cast<size_t>(time(nullptr) * info.n_rank * info.thread_per_process);

  std::string s = std::to_string(id_seed);
  info.run_id = std::hash<std::string>{}(s);
  check_results_file_name(params);
  return info;
}

void init_environment()
{

  const auto env_path = env_file_path();

  static const std::string cma_data_folder_path = "cma_data";
  static bool environment_initialized = false;

  if (!environment_initialized)
  {
    if (!std::filesystem::exists(cma_data_folder_path))
    {
      std::filesystem::create_directory(cma_data_folder_path);
    }

    if (!std::filesystem::exists(AutoGenerated::result_dir_path))
    {
      std::filesystem::create_directory(AutoGenerated::result_dir_path);
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

std::string sappend_date_time(std::string_view string) noexcept
{
  std::stringstream fd;
  fd << string;
  auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  fd << std::put_time(std::localtime(&now), "%Y-%m-%d-%H:%M:%S");
  return fd.str();
}

std::string env_file_path() noexcept
{
  return ".bmc_info";
}

void register_run(const ExecInfo& exec, const Core::UserControlParameters& /*params*/) noexcept
{
  // Open the file in append mode
  std::ofstream env(env_file_path(), std::ios_base::app);
  if (env.is_open())
  {
    append_date_time(env);
    env << "\t";
    env << exec;
    // env << params;
    env << std::endl;
  }
  else
  {
    std::cerr << "Error: Unable to open file for writing\n";
  }
}

void set_n_thread_current_rank(const int rank,
                               const int size,
                               ExecInfo& info,
                               const Core::UserControlParameters& params) noexcept
{
  // Casting rank and size to size_t
  info.current_rank = static_cast<size_t>(rank);
  info.n_rank = static_cast<size_t>(size);

  // Determining the number of OpenMP threads
  size_t omp_n_thread = (params.n_thread > 0) ? static_cast<size_t>(params.n_thread)
                                              : static_cast<size_t>(omp_get_max_threads());
  // int num_core_per_node = omp_get_num_procs();

  size_t threads_per_process = 1;

  if (omp_n_thread >= info.n_rank)
  {
    threads_per_process = omp_n_thread / info.n_rank;
    if (omp_n_thread % info.n_rank != 0 && info.current_rank == info.n_rank - 1)
    {
      threads_per_process += omp_n_thread % info.n_rank;
    }
  }

  info.thread_per_process = threads_per_process;

  assert(info.thread_per_process > 0);

  omp_set_num_threads(static_cast<int>(info.thread_per_process));
}

void check_results_file_name(Core::UserControlParameters& params)
{
  // TODO

  std::filesystem::path path = AutoGenerated::result_dir_path;
  std::filesystem::path folder_name = AutoGenerated::result_dir_path;

  if (params.results_file_name.empty())
  {
    folder_name.append(sappend_date_time("result_"));
    path.append(sappend_date_time("result_"));

    path.append(sappend_date_time("result_"));
  }
  else
  {
    folder_name.append(params.results_file_name);
    path.append(params.results_file_name);
    path.append(params.results_file_name);
  }
  std::filesystem::create_directory(folder_name);
  params.results_file_name = path.string();
}