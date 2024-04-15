
#include "messages/message_t.hpp"
#include <rt_init.hpp>

#include <common/common.hpp>
#include <cstdlib>
#include <filesystem>
#include <messages/wrap_mpi.hpp>
#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <Eigen/Core>

void init_bmc_info(std::ofstream& fd);

void set_openmp_threads(int rank,
    int size,
    ExecInfo& info,
    SimulationParameters& params)
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

    if (omp_n_thread >= info.n_rank) {
        threads_per_process = omp_n_thread / info.n_rank;
        if (omp_n_thread % info.n_rank != 0 && info.current_rank == info.n_rank - 1) {
            threads_per_process += omp_n_thread % info.n_rank;
        }
    }   
    
    if (threads_per_process > static_cast<size_t>(num_core_per_node)) {
        threads_per_process = static_cast<size_t>(num_core_per_node);
    }
    info.thread_per_process = threads_per_process;

    omp_set_num_threads(static_cast<int>(info.thread_per_process));
}

ExecInfo runtime_init(int argc, char** argv, SimulationParameters& params)
{
    ExecInfo info;

    int rank, size;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    set_openmp_threads(rank, size, info, params);

    std::cout << "NUM thread per process " << info.thread_per_process
              << std::endl;
    std::atexit(MPI_W::finalize);

    // MPI_W::is_mpi_init = true;
    return info;
}



void init_environement() {
    if (std::filesystem::exists(".bmc_info")) {
        return;
    }

    std::cout << "Setting up environment ...";

    std::filesystem::create_directory("cma_data");
    std::filesystem::create_directory("results");

    std::ofstream bmc_info_file(".bmc_info");
    if (bmc_info_file.is_open()) {
        std::filesystem::permissions(".bmc_info", std::filesystem::perms::owner_read);
        init_bmc_info(bmc_info_file);
        bmc_info_file.close();
        std::cout << "Done" << std::endl;
    } else {
        std::cout << "Failed to create .bmc_info file" << std::endl;
    }
}


void init_bmc_info(std::ofstream& fd) {
    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    fd  << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S") << std::endl;
}