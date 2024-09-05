from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import numpy as np

from apps.post_process.io import append_resukts_scalar_vtk
from . import mkdir
from typing import Dict, List,Optional
from . import RATIO_MASS_LENGTH

def get_distribution_moment(data):
    mean = np.mean(data)
    variance_population = np.var(data)
    variance_sample = np.var(data, ddof=1)

    return mean, variance_population, variance_sample


def _mk_pdf(data, name, dest: str):
    if len(data) == 0:
        return
    kde = gaussian_kde(data)
    x_values = np.linspace(min(data), max(data), 1000)
    pdf_values = kde(x_values)

    # Plot PDF
    plt.figure()
    plt.plot(x_values, pdf_values, color="blue")
    plt.xlabel("Value")
    plt.ylabel("Density")
    plt.title(f"Probability Density Function for {name} (PDF)")
    mkdir(f"{dest}/pdf")
    plt.savefig(f"{dest}/pdf/pdf_{name}")


def _mk_histogram(data, name, dest: str):
    num_bins = 100

    # plt.hist(data, bins=num_bins, density=True,alpha=0.7, color='blue', edgecolor='black')

    counts, bin_edges = np.histogram(data, bins=num_bins, density=False)

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

    plt.xlabel("Value")
    plt.ylabel("Density")
    plt.title(f"Histogram {name}")
    plt.savefig(f"{dest}/histogram//histogram_{name}")

    _mk_pdf(data, name, dest)


def property_distribution(
    biodict: Dict[str, np.ndarray],
    prefix: str = "",
    dest: str = "./results/",
    vtk_cma_mesh_path: Optional[str] = None,  # noqa: F821
):
    if(biodict is None):
        return
    for key in biodict:
        value = biodict[key]

        if key == "lenght":
            mass = np.sum(value) * RATIO_MASS_LENGTH
            print("mass: ", mass)

        if isinstance(value, np.ndarray) and np.issubdtype(value.dtype, float):
            mean, variance_population, variance_sample = get_distribution_moment(value)

            print(
                key,
                ": ",
                "mean: ",
                mean,
                "var: ",
                variance_population,
                "varred: ",
                variance_sample,
            )

            _mk_histogram(value, f"{prefix}_{key}", dest)
            if vtk_cma_mesh_path is not None:
                append_resukts_scalar_vtk(vtk_cma_mesh_path, value, key)


def property_space(i: int, biodict: Dict[str, np.ndarray], key1: str, key2: str):
    if key1 in biodict and key2 in biodict:
        value1 = biodict[key1]
        value2 = biodict[key2]
        MAX_SAMPLE = 1_000
        sample_size = min(len(value2), MAX_SAMPLE)  # or any smaller number
        idx = np.random.choice(sample_size, size=sample_size, replace=False)

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

    plt.savefig(f"{dest}/plot_{key1}_{key2}_{0}")


def process_particle_data(
    biodicts: List[Dict[str, np.ndarray]], dest_root: str = "./results/"
):
    dest = f"{dest_root}/properties"
    for i in biodicts:
        if i is None :
            return

    keys = [k for k in biodicts[0].keys() if k != "spatial"]

    plot_property_space(biodicts, "mu", "lenght", dest_root)

    mean_samples = {k: np.zeros(len(biodicts)) for k in keys}

    for i, bio_dict in enumerate(biodicts):
        for k in keys:
            if len(bio_dict[k]) > 0:
                if not np.isnan(bio_dict[k][0]):
                    mean_samples[k][i] = np.mean(bio_dict[k])
    mkdir(dest)
    for k, values in mean_samples.items():
        plt.figure()
        plt.plot(values, label=k)
        plt.xlabel("Sample Index")
        plt.ylabel("Mean Value")
        plt.title(f"Mean Value of {k} Over Samples")
        plt.legend()
        plt.savefig(f"{dest}/{k}.png")
        plt.close()  # Close the plot to free memory

