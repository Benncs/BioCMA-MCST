import sys
sys.path.insert(0, './builddir/pyrelease/apps/libs/pymodule')

import numpy as np 
import pyBioCMAMCST
import random

avg_mass = 1e-13
max_mass = 1e-11
phi_pts_max = 4.454e-12*1e-3/3600*0.1
std_dev = avg_mass / 10.
k_pts = 1e-3
YXS = 0.3
interdivision_time = 20*60


def cast_from_ptr(particle):
    opaque = particle.getOpaque()
    return opaque.cast()

def rand_division_time():
        return interdivision_time*np.random.normal(0.5,0.1)

class Foo:
    def __init__(self) -> None:
        self.mass = 1
        self.age =0
        self.length=6e-6
        self.phi =0
        self.interdivision_time = 0
    

def init_kernel(particle):
    opaque = particle.getOpaque()
    f = Foo()
    f.interdivision_time = rand_division_time()
    nrandom = random.gauss(avg_mass, std_dev)      
    f.mass = max(0., min(nrandom, max_mass))
    opaque.init(f)
   

def update_kernel(dt, particle, concentrations):
    opaque = particle.getOpaque()
    f = opaque.cast()
    f.age += dt
    s= concentrations[0]

    phi = max(0,(phi_pts_max * s / (k_pts + s)))
    f.phi = phi
    f.length+=dt*phi/YXS
    if f.age >= f.interdivision_time:
         particle.status = pyBioCMAMCST.CellStatus.CYTOKINESIS

   
def contribution_kernel(p, contribution):
    opaque = p.getOpaque()
    f = opaque.cast()
    ic = int(p.current_container)
    contribution[0, ic] -= f.phi * p.weight

def division_kernel(p, child):
    opaque_p = p.getOpaque()
    pmodel = opaque_p.cast()
    pmodel.age =0 
    pmodel.length = pmodel.length/2
    p.status = pyBioCMAMCST.CellStatus.IDLE
    
    f = Foo()
    f.age = 0
    f.length = pmodel.length
    f.interdivision_time =rand_division_time()
    child.status = pyBioCMAMCST.CellStatus.IDLE
    child.weight = p.weight
    opaque_child = child.getOpaque()
    opaque_child.init(f)
   

def get_properties(p):
    model = cast_from_ptr(p)

    return model.__dict__

def __debug(p):
    opaque = p.getOpaque()
    model = opaque.cast()
    print(model.l)