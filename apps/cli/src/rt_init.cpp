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


