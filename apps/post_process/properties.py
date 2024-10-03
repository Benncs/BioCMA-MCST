from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import numpy as np

from .read_results import Results, PartialResult

from .io import append_resukts_scalar_vtk
from . import FIGURE_TYPE, mkdir
from typing import Dict, List, Optional
from . import RATIO_MASS_LENGTH, get_time


def process_string_title(input_string):
    return " ".join(
        [word.capitalize() for word in input_string.replace("_", " ").split()]
    )


def get_distribution_moment(data):
    mean = np.mean(data)
    variance_population = np.var(data)
    variance_sample = np.var(data, ddof=1)

    return mean, variance_population, variance_sample


def _mk_pdf(data, name, dest: str):
    try:
        if len(data) == 0:
            return
        kde = gaussian_kde(data)
        x_values = np.linspace(min(data), max(data), 1000)
        pdf_values = kde(x_values)
        title_name = process_string_title(name)
        # Plot PDF
        plt.figure()
        plt.plot(x_values, pdf_values, color="blue")
        plt.xlabel("Value")
        plt.ylabel("Density")
        plt.title(f"Probability Density Function for {title_name} (PDF)")
        mkdir(f"{dest}/pdf")
        plt.savefig(f"{dest}/pdf/pdf_{name}{FIGURE_TYPE}")
    except:
        pass


def _mk_histogram(
    partials: List[PartialResult], key: str, index: int, name: str, dest: str
):
    num_bins = 100

    # data = partials[0].extra_bioparam[index][key]
    
    # counts, bin_edges = np.histogram(data, bins=num_bins, density=False)

    counts = np.zeros((num_bins,))
    bin_edges = np.zeros((num_bins + 1,))
    for i in partials:
        init = i.extra_bioparam[index]
        c, e = np.histogram(init[key], num_bins)
        counts += c
        bin_edges += e

    counts_normalized = counts / counts.max()
    mkdir(f"{dest}/histogram")
    plt.figure()
    plt.bar(
        bin_edges[:-1],
        counts_normalized,
        width=np.diff(bin_edges),
        edgecolor="black",
        alpha=0.7,
        color="blue",
    )
    title_name = process_string_title(name)
    plt.xlabel("Value")
    plt.ylabel("Density")
    plt.title(f"Histogram {title_name}")
    plt.savefig(f"{dest}/histogram/histogram_{name}{FIGURE_TYPE}")

    # _mk_pdf(data, name, dest)

    # plt.figure()
    # x = np.arange(1, len(data) + 1)
    # plt.scatter(x, data, color="blue", s=1, label="Temps de division")
    # plt.savefig(f"{dest}/histogram/test_{name}{FIGURE_TYPE}")
    # plt.close()


def property_distribution(
    partials: List[PartialResult],
    index: int,
    prefix: str = "",
    dest: str = "./results/",
    vtk_cma_mesh_path: Optional[str] = None,  # noqa: F821
):
    init = partials[0].extra_bioparam
    keys = [k for k in init[0].keys() if k != "spatial"]
    for key in keys:
        # if isinstance(value, np.ndarray) and np.issubdtype(value.dtype, float):
        # mean, variance_population, variance_sample = get_distribution_moment(value)

        # print(
        #     key,
        #     ": ",
        #     "mean: ",
        #     mean,
        #     "var: ",
        #     variance_population,
        #     "varred: ",
        #     variance_sample,
        # )
        _mk_histogram(partials, key, index, f"{prefix}_{key}", dest)

        # TODO Need to merge into value
        # if vtk_cma_mesh_path is not None:
        #     append_resukts_scalar_vtk(vtk_cma_mesh_path, value, key)


def property_space(i: int, biodict: Dict[str, np.ndarray], key1: str, key2: str):
    if key1 in biodict and key2 in biodict:
        value1 = biodict[key1]
        value2 = biodict[key2]
        MAX_SAMPLE = 1_000
        min_s = min(len(value1), len(value2))
        sample_size = min(min_s, MAX_SAMPLE)
        idx = np.random.choice(min_s, size=sample_size, replace=False)
        if isinstance(value1, np.ndarray) and np.issubdtype(value1.dtype, float):
            if isinstance(value2, np.ndarray) and np.issubdtype(value2.dtype, float):
                plt.scatter(value1[idx], value2[idx], label=f"data {i}", s=1)


def plot_property_space(
    biodicts: List[Dict[str, np.ndarray]],
    key1: str,
    key2: str,
    dest: str = "./results/",
):
    plt.figure()
    for i, d in enumerate(biodicts):
        property_space(i, d, key1, key2)
        plt.xlabel(key1)
        plt.ylabel(key2)

    plt.savefig(f"{dest}/plot_{key1}_{key2}_{0}{FIGURE_TYPE}")


def mean_partial(keys: List[str], partial: PartialResult):
    n_export = len(partial.extra_bioparam)
    sums = np.zeros((n_export, len(keys)))
    for i in range(n_export):
        current_export = partial.extra_bioparam[i]
        for j, key in enumerate(keys):
            sums[i][j] = np.sum(current_export[key])
    return sums


def plot_average(results: Results, dest: str):
    init = results.partial[0].extra_bioparam
    keys = [k for k in init[0].keys() if k != "spatial"]
    n_keys = len(keys)

    average = np.zeros((len(init), len(keys)))
    for i in results.partial:
        average += mean_partial(keys, i)

    average = average / results.total_repartion

    mkdir(dest)
    for i, key in enumerate(keys):
        try:
            plt.figure()
            plt.plot(
                results.time,
                average[:, i],
                "-o",
                label=key,
                markersize=2,
                color="black",
            )
            plt.xlabel(f"Time [{get_time()}]")
            plt.ylabel("Mean Value")
            plt.title(f"Mean Value of {key} Over Time")
            plt.legend()
            plt.savefig(f"{dest}/{key}{FIGURE_TYPE}")
            plt.close()  # Close the plot to free memory
        except Exception as e:
            plt.close()  # Close the plot to free memory
            print(f"ERROR Mean Value {key} Over Time")


def process_particle_data(results: Results, dest_root: str = "./results/"):
    dest = f"{dest_root}/properties"

    if results.partial[0].extra_bioparam is None:
        return

    property_distribution(
        results.partial,
        0,
        "init",
        dest,
    )

    property_distribution(
        results.partial,
        -1,
        "final",
        dest,
    )

    plot_average(results, dest)
