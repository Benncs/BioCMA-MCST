import matplotlib.pyplot as plt 
from read_results import import_results
import numpy as np 
from numpy.polynomial import Polynomial
import scienceplots

cn = "deter"

results = import_results(f"./results/0d/{cn}.h5")


initial_distribution = results.initial_distribution


t2 = np.linspace(0,results.tf+1200)


tau = 20*60 
mu = np.log(2)/tau

N_ref = initial_distribution[0]*np.exp(t2*mu)

log_n_ref = np.log(N_ref)



print(f"reference: log(y) = {initial_distribution[0]:.6f}+{mu:.6f}x")

y_exp = np.log(results.records_distribution[:, 0, 0])

# Fit a linear model to the log-transformed data
fitting = Polynomial.fit(results.t, y_exp, 1)
coeff = fitting.convert().coef

# Calculate the predicted values and residuals
y_pred = coeff[0] + coeff[1] * results.t

# Calculate the exponential predicted values for the original scale
y_pref_exo = np.exp(coeff[0]) * np.exp(coeff[1] * results.t)

# Calculate residuals and R^2
residuals = y_exp - y_pred
ss_res = np.sum(residuals**2)
ss_tot = np.sum((y_exp - np.mean(y_exp))**2)
r2 = 1 - (ss_res / ss_tot)

# Print reference and result equations
print(f"reference: log(y) = {initial_distribution[0]:.6f} + {mu:.6f}x")
print(f"results: log(y) = {coeff[0]:.6f} + {coeff[1]:.6f}x")
print(f"R^2: {r2:.6f}")

# Create the plot for log-transformed data
plt.style.use('ggplot')
fig, ax = plt.subplots()
ax.set_title("Population Growth as a Function of Time (Log Scale)")
ax.plot(t2, log_n_ref, label=f"reference: log(n) = {initial_distribution[0]:.6f} + {mu:.6f}x")
ax.plot(results.t, y_exp, label="results")
ax.plot(results.t, y_pred, label=f"regression: log(n) = {coeff[0]:.6f} + {coeff[1]:.6f}x\n$R^2$ = {r2:.6f}")
ax.set_xlabel("time [s]")
ax.set_ylabel("ln(N)")
ax.legend()
fig.savefig(f"results/0d/show_c_{cn}.svg", bbox_inches="tight")
plt.close()
# Create the plot for the actual population growth
plt.figure()
plt.title("Population Growth as a Function of Time")
plt.plot(results.t, results.records_distribution[:, 0, 0],"*", label="results")
plt.plot(results.t, y_pref_exo, label="predicted")
plt.xlabel("time [s]")
plt.ylabel("N")
plt.legend()
plt.savefig(f"results/0d/show_c_exp_{cn}.svg", bbox_inches="tight")
plt.close()