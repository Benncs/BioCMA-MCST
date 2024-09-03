
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import numpy as np 
from . import mkdir

def get_distribution_moment(data):
    mean = np.mean(data)
    variance_population = np.var(data)
    variance_sample = np.var(data, ddof=1)

    return mean,variance_population,variance_sample

def mk_pdf(data,name,dest:str):
    if len(data)==0:
        return
    kde = gaussian_kde(data)  
    x_values = np.linspace(min(data), max(data), 1000)
    pdf_values = kde(x_values)

    # Plot PDF
    plt.figure()
    plt.plot(x_values, pdf_values, color='blue')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.title(f"Probability Density Function for {name} (PDF)")
    mkdir(f"{dest}/pdf")
    plt.savefig(f"{dest}/pdf/pdf_{name}")

def mk_histogram(data,name,dest:str):
    num_bins = 100

    # plt.hist(data, bins=num_bins, density=True,alpha=0.7, color='blue', edgecolor='black')
    
    counts, bin_edges = np.histogram(data, bins=num_bins, density=False)

    counts_normalized = counts / counts.max()
    mkdir(f"{dest}/histogram")
    plt.figure()
    plt.bar(bin_edges[:-1], counts_normalized, width=np.diff(bin_edges), edgecolor='black', alpha=0.7, color='blue')

    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.title(f"Histogram {name}")
    plt.savefig(f"{dest}/histogram//histogram_{name}")

    mk_pdf(data,name,dest)




