"""
Script to perform bench scaling and profiling of an external running program.

Author: CASALE Benjamin
Date: 10/02/2025
Version: 2.0
"""

from typing import Dict
import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import sys
import glob
import shutil
import json
from matplotlib.backends.backend_pdf import PdfPages
from cycler import cycler


def format_with_k_notation(num):
    print(num)
    if num >= 1000_000:
        return f"{num / 1_000_000:.1f}M"
    elif num >= 1_000:  # For thousands
        return f"{num / 1_000:.1f}K"
    else:
        return str(num)


def add_to_pdf(_pdf, figs):
    for fig in figs:
        _pdf.savefig(fig)
        plt.close(fig)


def default_config():
    bench_config = {
        "OMP_THREADS": [1, 6],
        "mpi": None,
        "path": os.environ.get("BIOMC_ROOT"),
        "executable_name": "biocma_mcst_cli_app",
        "name": "host",
        "folder": "./bench_2",
    }

    case_config = {
        "model_name": "two_meta",
        "tf": 10,
        "dt": 1e-2,
        "cma_path": "/home_pers/casale/Documents/thesis/cfd/sanofi/",
        "cma_init": "./cma_data/sanofi_init.h5",
        "recursive": False,
    }

    return {"bench_config": bench_config, "case_config": case_config}


def read_config(path: str) -> Dict:
    with open(path, "r") as fd:
        data = json.load(fd)
    return data


def annotate_json(filename, bench_config, n_p, n_threads):
    dest = f"bench_{bench_config['name']}_{n_p}_{n_threads}.json"
    dest = f"{bench_config['folder']}/{dest}"
    shutil.move(filename, dest)
    with open(dest, "r") as file:
        data = json.load(file)
    data["n_p"] = str(n_p)
    data["n_threads"] = str(n_threads)
    with open(dest, "w") as file:
        json.dump(data, file, indent=2)


def format_cli(bench_config, case_config, number_particle, nt=None):
    rec = []

    if case_config["recursive"]:
        rec = [
            "-r",
            "1",
        ]

    common_cmd = [
        "-mn",
        case_config["model_name"],
        "-np",
        f"{number_particle}",
        "-d",
        f"{case_config['tf']}",
        "-dt",
        f"{case_config['dt']}",
        *rec,
        "-f",
        case_config["cma_path"],
        "-er",
        f"./tmp_result/{bench_config['name']}",
        "-nex",
        "0",
        "-force",
        "1",
    ]

    if case_config["cma_init"] is not None:
        common_cmd.append("-fi")
        common_cmd.append(case_config["cma_init"])

    if bench_config.get("mpi",None) is None:
        cmd = [f"{bench_config['path']}/{bench_config['executable_name']}", *common_cmd]
    else:
        n_node = bench_config["mpi"]
        cmd = [
            "mpiexec",
            "-np",
            str(n_node),
            "--bind-to",
            "core",
            "--map-by",
            f"slot:PE={nt} ",
            f"{bench_config['path']}/{bench_config['executable_name']}",
            *common_cmd,
        ]
    return cmd


def read_dict(bench_config):
    data = []
    pattern = f"{bench_config['folder']}/bench_{bench_config['name']}_*.json"
    matches = glob.glob(pattern)
    main_dict = {}
    for i, filename in enumerate(matches):
        with open(filename, "r") as fd:
            data = json.load(fd)
            main_dict[i] = data
    return main_dict


def execute(bench_config, n_thread, n_p, command):
    env_var = os.environ.copy()
    env_var["OMP_NUM_THREADS"] = str(n_thread)
    env_var["OMP_PLACES"] = "threads"
    env_var["OMP_PROC_BIND"] = "spread"
    env_var["KOKKOS_TOOLS_LIBS"] = (
        "/home_pers/casale/Documents/code/kokkos-tools/build/lib/libkp_kernel_timer.so"
    )
    env_var["KOKKOS_TOOLS_TIMER_JSON"] = "1"
    commands = command
    os.makedirs(bench_config["folder"], exist_ok=True)
    result = subprocess.Popen(
        commands,
        env=env_var,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    result.wait()
    if result.returncode != 0:
        print(result.stderr)
        raise Exception("error")
    filename = f"{os.uname().nodename}-{result.pid}.json"

    n_err = 0
    while n_err != 10:
        if os.path.exists(filename):
            break
        n_err += 1
        filename = f"{os.uname().nodename}-{result.pid + n_err}.json"

    annotate_json(filename, bench_config, n_p, n_thread)


def scaling(case_config, bench_config, n_p):
    records = []
    command_to_run = format_cli(bench_config, case_config, n_p)
    print("Running: ", *command_to_run, " ....")
    for n_thread in bench_config["OMP_THREADS"]:
        command_to_run = format_cli(bench_config, case_config, n_p, n_thread)
        print("OMP_NUM_THREADS=", n_thread)
        run_time = execute(bench_config, n_thread, n_p, command_to_run)
        records.append(run_time)
    return records


def find_in_data(data, name, key, region=False):
    subd = "kernel-perf-info"
    if region:
        subd = "region-perf-info"
    return np.array(
        [
            next(
                (
                    d[key]
                    for d in data[i]["kokkos-kernel-data"][subd]
                    if d["kernel-name"] == name
                ),
                0,
            )
            for i in data
        ]
    )


def process_results(config):
    # Read data from the file
    bench_config = config["bench_config"]
    data = read_dict(bench_config)

    res = {}
    # Extract data for plotting
    res["total_app_time"] = np.array(
        [int(data[r]["kokkos-kernel-data"]["total-app-time"]) for r in data]
    )
    res["threads"] = np.array([int(data[r]["n_threads"]) for r in data])
    res["particles"] = np.array([int(data[r]["n_p"]) for r in data])
    res["target_kernel"] = "cycleProcess"
    res["call_count"] = find_in_data(data, "cycleProcess", "call-count", True)
    # target_kernel = "cycle_model"
    #    time_model = find_in_data(target_kernel, "total-time", False)
    #  target_kernel = "cycle_move"
    #   time_move = find_in_data(target_kernel, "total-time", False)
    res["total_time"] = find_in_data(
        data, "cycleProcess", "total-time", True
    )  # time_model + time_move
    return res


## FIGURE


def plot_scaling(total_time, threads, particles, iterations, records):
    unique_particles = np.unique(particles)
    unique_iterations = np.unique(iterations)
    font_properties = {"family": "serif", "size": 14}
    fig, (ax_scale, ax_eff) = plt.subplots(1, 2, figsize=(14, 6))
    fig2, ax_percent = plt.subplots(1, 1, figsize=(14, 6))
    figs = [fig, fig2]
    colors = plt.cm.tab20b(np.linspace(0, 1, 10))
    line_widths = [1] * 10
    custom_cycler = cycler(color=colors) + cycler(lw=line_widths)

    ax_scale.set_prop_cycle(custom_cycler)
    ax_eff.set_prop_cycle(custom_cycler)

    # Plot strong scaling efficiency (speedup) vs. number of threads for each particle and iteration combination
    for particle in unique_particles:
        for iteration in unique_iterations:
            # f = plt.figure()
            mask = (particles == particle) & (iterations == iteration)
            condition = True  # np.count_nonzero(mask) >= 3
            if condition:
                # Extract the relevant data
                selected_threads = threads[mask]
                selected_times = records[mask]
                selected_total_times = total_time[mask]

                # Find the time for single-thread execution (assuming we have this data point)
                single_thread_time = selected_times[selected_threads == 1]

                if single_thread_time.size == 0:
                    raise Exception("There's no single-thread-data")
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

                percent = selected_times / selected_total_times
                ax_scale.plot(
                    sorted_threads,
                    sorted_speedup,
                    "-x",
                    label=f"{particle:.0e}",
                    markersize=8,
                    linewidth=2,
                )
                ax_eff.plot(
                    sorted_threads,
                    sorted_speedup / sorted_threads,
                    "-x",
                    label=f"{particle:.0e}",
                    markersize=8,
                    linewidth=2,
                )

                ax_percent.plot(
                    sorted_threads,
                    percent[sorted_indices],
                    "-x",
                    label=f"{particle:.0e}",
                    markersize=8,
                    linewidth=2,
                )

    ax_scale.set_title(
        f"Speedup for {iteration} iterations", fontsize=16, fontdict=font_properties
    )
    ax_eff.set_title(
        f"Efficiency for {iteration} iterations", fontsize=16, fontdict=font_properties
    )
    ax_scale.set_xlabel("Number of Threads", fontsize=14, fontdict=font_properties)
    ax_eff.set_xlabel("Number of Threads", fontsize=14, fontdict=font_properties)
    ax_eff.set_ylabel("Efficiency", fontsize=14, fontdict=font_properties)
    ax_scale.set_ylabel("Speedup", fontsize=14, fontdict=font_properties)

    ax_scale.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax_eff.grid(True, which="both", linestyle="--", linewidth=0.5)

    x = np.linspace(0, np.max(threads))
    ax_scale.plot(x, x, label="ideal", color="k", linewidth=2)
    ax_scale.legend(
        title="n particle",
        bbox_to_anchor=(1.05, 0),
        loc="lower left",
        borderaxespad=0.0,
    )

    ax_percent.set_title(
        f"Percent in kernels for {iteration} iterations",
        fontsize=16,
        fontdict=font_properties,
    )

    ax_percent.set_xlabel("Number of Threads", fontsize=14, fontdict=font_properties)

    ax_percent.legend(
        title="n particle",
        borderaxespad=0.0,
    )

    fig.tight_layout()
    return figs


def plot_results(config, processed_results):
    bench_config = config["bench_config"]
    total_app_time = processed_results["total_app_time"]
    threads = processed_results["threads"]
    particles = processed_results["particles"]
    call_count = processed_results["call_count"]
    total_time = processed_results["total_time"]
    with PdfPages(f"{bench_config['folder']}/{bench_config['name']}.pdf") as pdf:
        add_to_pdf(
            pdf,
            plot_scaling(total_app_time, threads, particles, call_count, total_time),
        )


## MAIN


def _fom2(args):
    case = read_config(args[1])
    bc = case["bench_config"]
    data = read_dict(bc)
    keys = list(data.keys())

    def plot_fom():
        figs = []
        for i in range(0, len(keys), len(keys) // 4):
            datai = data[keys[i]]
            region_perf_info = datai.get("kokkos-kernel-data", {}).get(
                "region-perf-info", []
            )
            kernel_names = [entry["kernel-name"] for entry in region_perf_info]
            total_times = [entry["total-time"] for entry in region_perf_info]
            print(kernel_names)
            sorted_kernels = sorted(
                zip(kernel_names, total_times), key=lambda x: x[1], reverse=True
            )

            sorted_kernel_names, sorted_total_times = zip(*sorted_kernels)

            figs.append(plt.figure(figsize=(10, 6)))
            plt.barh(sorted_kernel_names, sorted_total_times, color="skyblue")
            plt.xlabel("Total Time (s)")

            plt.title(
                f"Most Time Consuming Kernels for {format_with_k_notation(int(datai['n_p']))} particle"
            )
            plt.gca().invert_yaxis()
        return figs

    with PdfPages(f"{bc['folder']}/{bc['name']}_fom.pdf") as pdf:
        add_to_pdf(
            pdf,
            plot_fom(),
        )


def _fom(args):
    cases = [read_config(args[i]) for i in range(1, len(args))]
    names = []
    all_np = []
    all_kt = []
    all_time = []
    all_nt = []
    all_mpi = []

    for i, case in enumerate(cases):
        bc = case["bench_config"]
        data = read_dict(bc)
        names.append([bc["name"]] * len(data))
        res = np.array([int(data[r]["n_p"]) for r in data])
        all_np.append(res)
        kt = find_in_data(data, "cycleProcess", "total-time", True)
        all_time.append(
            np.array([data[i]["kokkos-kernel-data"]["total-app-time"] for i in data])
        )
        all_kt.append(kt)

        all_nt.append(
            np.array( [int(data[i]["n_threads"]) for i in data])
        )
    
        all_mpi.append(
            np.array( [bc.get('mpi',0) for i in data])
        )
        
    all_np = np.concatenate(all_np)
    all_time = np.concatenate(all_time)
    names = np.concatenate(names)
    all_kt = np.concatenate(all_kt)
    all_nt = np.concatenate(all_nt)
    maxnt = np.max(all_nt)
    all_mpi = np.concatenate(all_mpi)

    maxnp = 100e6 #np.max(all_np)
    mask_name = True  # (names == "gpe_3d") | (names == "gpu_3d")
    mask_np = (all_np == maxnp) | (all_np == maxnp)
    mask_mpi = (all_nt == 1) & (all_mpi == 20)
    mask_thread = (all_nt == maxnt) | mask_mpi
    mask = mask_np & mask_name & mask_thread 
    

    fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    ax[0].bar(names[mask], all_kt[mask], color='skyblue')
    ax[0].set_title("Kernel Time (100M Particles)", fontsize=14)
    ax[0].set_ylabel("Time [s]", fontsize=12)
    ax[0].grid(True, which='both', linestyle='--', alpha=0.7)
    ax[1].bar(names[mask], all_time[mask], color='#FF6347')
    ax[1].set_title("Wall Time (100M Particles)", fontsize=14)
    ax[1].set_ylabel("Time [s]", fontsize=12)
    ax[1].set_xlabel("Names", fontsize=12)
    ax[1].grid(True, which='both', linestyle='--', alpha=0.7)
    plt.tight_layout()

    plt.show()

    return 0


def _plot(args):
    if len(args) != 2:
        raise Exception("bad arguments: case ")
    bench_config = read_config(args[1])
    res = process_results(bench_config)
    plot_results(bench_config, res)
    return 0


def _gen(args):
    with open("bench.json", "w") as fd:
        json.dump(default_config(), fd, indent=2)
    return 0


def _add(args):
    if len(args) != 5:
        raise Exception("bad arguments: expected case path np nt ")
    case = args[1]
    path = args[2]
    bc = process_results(read_config(case))
    annotate_json(bc, path, int(args[3]), int(args[4]))
    return 0


def _do_scale(args):
    os.makedirs("tmp_result", exist_ok=True)
    if len(args) != 5:
        print(
            "Error: Invalid number of arguments for scale. Usage: python script.py scale [particle_n1] [particle_n2] [n_scale] [case]"
        )
        return

    particle_n1 = float(args[1])
    particle_n2 = float(args[2])
    n_scale = int(args[3])
    config_path = args[4]
    config = read_config(config_path)
    n_particles = np.linspace(particle_n1, particle_n2, num=n_scale, dtype=np.int32)
    for number in n_particles:
        scaling(
            config["case_config"],
            config["bench_config"],
            number,
        )


if __name__ == "__main__":
    _cb = {
        "plot": _plot,
        "gen": _gen,
        "add": _add,
        "fom": _fom,
        "scale": _do_scale,
        "fom_k": _fom2,
    }
    args = sys.argv[1:]

    if args[0] in _cb:
        _cb[args[0]](args)
        sys.exit(0)
    else:
        print("Error: Invalid argument. Usage: python script.py [plot | scale]")
        sys.exit(1)
