import sys
sys.path.insert(0, './builddir/debug/apps/libs/pymodule')

import pyBioCMAMCST
import numpy as np 

import random;


def rand():
        return random.random()

class Foo:
    a = 1

def init_kernel(particle):
    opaque = particle.getOpaque()
    print(opaque)
    f = Foo()
    f.a = np.random.random()
    opaque.init(f)
   

def update_kernel(dt, particle, concentrations):
    opaque = particle.getOpaque()
    f = opaque.cast()
    print(f.a*100)
   
def contribution_kernel(p, contribution):
    # print("Initializing division_kernel")
    pass

def init_contribution_kernel(particle, contrib):
    print("Initializing contribution_kernel")
