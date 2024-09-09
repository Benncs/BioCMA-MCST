"""
Script to perform bench scaling and profiling of an external running program.

Author: CASALE Benjamin
Date: 04/26/2024
Version: 1.0
"""

import os
import csv
import subprocess
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import sys
# from wakepy import keep
# from ..exec import exec
import datetime 
BENCH_OMP_THREADS = [1, 6,12]  # List of thread numbers when running scaling
EXECUTABLE_PATH = "./builddir/release_gcc/apps/cli"  # Path to executable to run
EXECUTABLE_NAME = "biocma_mcst_cli_app_shared"  # Name of executable to run
BENCH_SCRIPT_PATH = (
    "./devutils/benchs/bench.sh"  # Intermediate script used to perform bench
)
date = datetime.datetime.today().strftime('%Y_%m_%d')
FILENAME = f"./devutils/benchs/bench_records_{date}.csv"  # Record filename
OUTPUT_PDF = f"./devutils/benchs/results_bench_{date}.pdf"  # Output path
MODEL_NAME = "default"
FINAL_TIME = 1  # Reference simulation time
DELTA_TIME = 1e-3  # Reference delta time fixed
# CMA_DATA_PATH = "./cma_data/bench/"
CMA_DATA_PATH = "./cma_data/0d/"
RECURSIVE = False

def format_cli(number_particle, final_time):

    rec  = [""]
    if RECURSIVE:
        rec = [    "-r",
        "1",]

    return [
        f"{EXECUTABLE_PATH}/{EXECUTABLE_NAME}",
        "-mn",
        MODEL_NAME,
        "-np",
        f"{number_particle}",
        "-d",
        f"{final_time}",
        "-dt",
        f"{DELTA_TIME}",
        *rec,
        "-f",
        CMA_DATA_PATH,
        "-er",
        "bench",
    ]


def execute(n_thread, script_path, command):
    env_var = os.environ.copy()
    env_var["OMP_NUM_THREADS"] = str(n_thread)
    env_var["OMP_PLACES"] = "threads"
    env_var["OMP_PROC_BIND"] = "spread"
    # commands = [
    #     "mpiexec",
    #     "--allow-run-as-root",
    #     "-n",
    #     "6",
    #     BENCH_SCRIPT_PATH,
    #     *command,
    # ]

    commands = [
        BENCH_SCRIPT_PATH,
        *command,
    ]

    result = subprocess.run(
        commands,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        env=env_var,
    )
    if result.returncode != 0:
        raise Exception("error", result.stderr)
    time_resutlt = result.stderr

    match = re.search(r"Time elapsed: (\d+\.\d+) seconds", time_resutlt)
    if match:
        time_elapsed = float(match.group(1))
        return time_elapsed
    else:
        print("Time elapsed not found in the output.")


def scalling(np, tf):
    records = []
    commanDELTA_TIMEo_run = format_cli(np, tf)
    print("Running: ", *commanDELTA_TIMEo_run, " ....")

    for n_thread in BENCH_OMP_THREADS:
        print("OMP_NUM_THREADS=", n_thread)
        run_time = execute(n_thread, BENCH_SCRIPT_PATH, commanDELTA_TIMEo_run)
        records.append(run_time)

    return records


def dump_to_csv(records, particles, n_iteration):
    # Create a CSV file with a header row
    is_empty = not os.path.exists(FILENAME) or os.path.getsize(FILENAME) == 0
    with open(FILENAME, mode="a") as file:
        writer = csv.writer(file)
        if is_empty:
            writer.writerow(["Thread", "particles", "iteration", "Record"])
        # Write each record to the CSV file
        for i, record in enumerate(records):
            writer.writerow([BENCH_OMP_THREADS[i], particles, n_iteration, record])


def read_csv_to_dict():
    # Read the CSV file into a dictionary
    data = []
    with open(FILENAME, mode="r") as file:
        reader = csv.DictReader(file)
        for row in reader:
            data.append(row)
    return data


def plot_thread_vs_time(threads, particles, iterations, records):
    unique_particles = np.unique(particles)
    unique_iterations = np.unique(iterations)
    figs = []
    # Plot time vs. number of threads for each particle
    for particle in unique_particles:
        for iteration in unique_iterations:
            f = plt.figure()
            plt.title(
                f"Time=f(num_thread) for n_particle={particle}, iteration={iteration}"
            )
            mask = (particles == particle) & (iterations == iteration)
            if np.count_nonzero(mask) > 3:
                plt.plot(threads[mask], records[mask], "-o")
                plt.xlabel("Number of Threads")
                plt.ylabel("Time (s)")
                plt.grid(True)
                figs.append(f)
    return figs


def plot_scaling(threads, particles, iterations, records):
    unique_particles = np.unique(particles)
    unique_iterations = np.unique(iterations)
    figs = []

    # Plot strong scaling efficiency (speedup) vs. number of threads for each particle and iteration combination
    for particle in unique_particles:
        for iteration in unique_iterations:
            f = plt.figure()
            plt.title(
                f"Speedup=f(num_thread) for n_particle={particle}, iteration={iteration}"
            )
            mask = (particles == particle) & (iterations == iteration)
            if np.count_nonzero(mask) > 3:
                # Extract the relevant data
                selected_threads = threads[mask]
                selected_times = records[mask]

                # Find the time for single-thread execution (assuming we have this data point)
                single_thread_time = selected_times[selected_threads == 1]

                if single_thread_time.size == 0:
                    continue  # Skip if there's no single-thread data

                single_thread_time = single_thread_time[
                    0
                ]  # There should be exactly one such entry

                # Compute speedup
                speedup = single_thread_time / selected_times

                # Sort by the number of threads for a proper plot
                sorted_indices = np.argsort(selected_threads)
                sorted_threads = selected_threads[sorted_indices]
                sorted_speedup = speedup[sorted_indices]

                # Plot the speedup
                plt.plot(sorted_threads, sorted_speedup, "-o")
                plt.xlabel("Number of Threads")
                plt.ylabel("Speedup (T(1)/T(n))")
                plt.grid(True)
                figs.append(f)

    return figs


def plot_particle_vs_time(threads, particles, iterations, records):
    unique_threads = np.unique(threads)
    unique_iterations = np.unique(iterations)
    figs = []
    for thread in unique_threads:
        for iteration in unique_iterations:
            f = plt.figure()
            mask = thread == threads
            if np.count_nonzero(mask) > 3:
                plt.title(
                    f"Time=f(particles) for n_threads={thread}, n_iteration={iteration}"
                )
                plt.plot(particles[mask], records[mask], "-o")
                plt.xlabel("Number of Particles")
                plt.ylabel("Time (s)")
                plt.grid(True)
                figs.append(f)

    return figs


def add_to_pdf(_pdf, figs):
    for fig in figs:
        _pdf.savefig(fig)
        plt.close(fig)


def plot_csv():
    # Read data from the CSV file
    data = read_csv_to_dict()

    # Extract data for plotting
    threads = np.array([int(row["Thread"]) for row in data])
    particles = np.array([float(row["particles"]) for row in data])
    iterations = np.array([float(row["iteration"]) for row in data])
    records = np.array([float(row["Record"]) for row in data])
    total_figs = []

    pre_mask = particles >= 10000

    with PdfPages(OUTPUT_PDF) as pdf:
        add_to_pdf(
            pdf,
            plot_thread_vs_time(
                threads[pre_mask],
                particles[pre_mask],
                iterations[pre_mask],
                records[pre_mask],
            ),
        )
        #   add_to_pdf(pdf,plot_particle_vs_time(threads,particles,iterations,records))
        add_to_pdf(pdf, plot_scaling(threads, particles, iterations, records))


def do_scale(particles=1000000):
    records = scalling(particles, FINAL_TIME)
    n_iteration = FINAL_TIME / DELTA_TIME
    dump_to_csv(records, particles, n_iteration)


def main(args):
    if args[0] == "plot":
        plot_csv()
    elif args[0] == "scale":
        if len(args) != 4:
            print(
                "Error: Invalid number of arguments for scale. Usage: python script.py scale [particle_n1] [particle_n2] [n_scale]"
            )
            return

        particle_n1 = float(args[1])
        particle_n2 = float(args[2])
        n_scale = int(args[3])
        try:
            raise Exception("")
            # with keep.running():
            #     n_particles = np.linspace(particle_n1, particle_n2, num=n_scale, dtype=np.int32)
            #     for number in n_particles:
            #         do_scale(number)
        except:
            n_particles = np.linspace(
                particle_n1, particle_n2, num=n_scale, dtype=np.int32
            )
            for number in n_particles:
                do_scale(number)
    else:
        print("Error: Invalid argument. Usage: python script.py [plot | scale]")


if __name__ == "__main__":
    if len(sys.argv) not in [2, 5]:
        print(
            "Error: Invalid number of arguments. Usage: python script.py [plot | scale]"
        )
        sys.exit(1)

    args = sys.argv[1:]
    main(args)
