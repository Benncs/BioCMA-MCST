import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

import sys 

se = 10
if(len(sys.argv)==2):
    se = float(sys.argv[1])


    
yx = 2    
mumax = 0.77/3600     
v = 20e-3
taum = 1/mumax
Q =   mumax*v
# l = np.linspace(0.1e-6, 5e-6, 1000)
# l0 = 0.8e-6
# l1 = 2e-6

# # Calcul de la probabilité
# p = (l - l0) / (l1 - l0)
# p[l < l0] = 0
# p[l > l1] = 1

# plt.plot(l, p, label='Probabilité')

# plt.show()
# exit(0)
Ks = 0.01
def model(t, y):
    s, x,mu = y  
    mup = mumax*s/(Ks+s)
    mueff = min(mup,mu)
    D = Q/v
    dsdt = D*(se-s)-mueff*yx*x
    dxdt = (mueff-D)*x
    dmudt = 1/taum*(mup-mu)
    return np.array([dsdt, dxdt,dmudt])  

s0 = 0
x0 = 1
mu = mumax/2
initial_conditions = [s0, x0,mu]
# tau =  v/Q
t_span = (0, 75000)  
t_eval = np.linspace(t_span[0], t_span[1], 500)  

solution = solve_ivp(model, t_span, initial_conditions, t_eval=t_eval,method="BDF")
initial_mass_cell = 3.14 * (0.8e-6) * (0.8e-6)/4. * 0.9e-6/2. * 1000



t = solution.t/3600  
s = solution.y[0]  
x = solution.y[1] 
m = solution.y[2] 


mx0 = v*x0
mxf = v*x[-1]

n0=mx0/initial_mass_cell
nf=mxf/initial_mass_cell
# print(nf/n0)


print(s[-1],x[-1],m[-1])

plt.figure(figsize=(12, 6))
plt.plot(t, s, label='Substrate Concentration (s)')
plt.ylabel('Concentration of Substrate')
plt.legend()

plt.figure(figsize=(12, 6))
plt.plot(t, x, label='Biomass Concentration (x)', color='orange')
plt.ylabel('Concentration of Biomass')
plt.xlabel('Time')
plt.legend()

plt.figure(figsize=(12, 6))
plt.plot(t, m, label='Growth rate', color='orange')
plt.ylabel('Growth rate')
plt.xlabel('Time')
plt.legend()


plt.tight_layout()
plt.show()
