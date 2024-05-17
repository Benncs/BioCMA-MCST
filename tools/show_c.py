import matplotlib.pyplot as plt 
from read_results import import_results
import numpy as np 

initial_distribution,distribution,cliq,data,tf,dt,npart,records_distribution = import_results("./results/test.h5")


n_t = data.shape[0]
t = np.linspace(0,tf,n_t)



tau = 20*60 
mu = np.log(2)/tau

N = initial_distribution[0]*np.exp(t*mu)

plt.plot(t,N)
plt.plot(t,records_distribution[:,0,0])
plt.show()


# mean = tau
# stddev = tau/2
# # x = np.linspace(mean - 3*stddev, mean + 3*stddev, 1000)
# x = np.linspace(0,tf,1000)
# y = (1/(stddev * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / stddev) ** 2)
# plt.ylim(0,np.max(y))
# plt.plot(x,y)
# plt.show()