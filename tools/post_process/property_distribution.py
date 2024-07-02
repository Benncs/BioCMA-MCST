from read_results import import_results
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import numpy as np 


def get_distribution_moment(data):
    mean = np.mean(data)
    variance_population = np.var(data)
    variance_sample = np.var(data, ddof=1)

    return mean,variance_population,variance_sample

def mk_histogram(data,name):
    num_bins = 150
    plt.figure()
    plt.hist(data, bins=num_bins, density=True, alpha=0.7, color='blue', edgecolor='black')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.title('Histogram')
    plt.savefig(f"./results/histogram_{name}")


def mk_pdf(data,name):
    kde = gaussian_kde(data)  
    x_values = np.linspace(min(data), max(data), 1000)
    pdf_values = kde(x_values)

    # Plot PDF
    plt.figure()
    plt.plot(x_values, pdf_values, color='blue')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.title('Probability Density Function (PDF)')
    plt.savefig(f"./results/pdf_{name}")


res = import_results('./results/test_distrib.h5')

for key in res.bioparam:
    
    value = res.bioparam[key]
    
    if isinstance(value, np.ndarray) and np.issubdtype(value.dtype, float):
       
        mean,variance_population,variance_sample = get_distribution_moment(value)

        print(mean,variance_population,variance_sample)

        mk_histogram(value,key)

        # mk_pdf(value,key)