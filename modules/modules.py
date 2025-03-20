import cpp_module
import random


def __start():

    import sys

    # Get the current module
    current_module = sys.modules[__name__]

    if not hasattr(current_module, 'init'):
        raise RuntimeError("update does not exist.")

    if not hasattr(current_module, 'update'):
        raise RuntimeError("update does not exist.")

    if not hasattr(current_module, 'division'):
        raise RuntimeError("division does not exist.")

    print(f"Module {__name__} Loaded")



mu_max = 0.77 / 3600.
ks=0.01
l_0=0.9e-6
ln2=0.69314718056
_init_only_cell_lenghtening = l_0/2. / ln2

tau_metabolism=.1/mu_max

class Model:
    def __init__(self):
        self.l= 0
        self.mu =0
        self.cc=0

def init(p:cpp_module.ParticleDataHolder):
    model =Model()
    model.l=random.normalvariate(0.9e-6/2,0.9e-6/5)
    model.mu = random.normalvariate(mu_max,mu_max/5)
    return model

def update(dt:float,_data:cpp_module.OpaquePointer,p:cpp_module.ParticleDataHolder,concentrations):   
    model =_data.get()
    s = concentrations[0]
    mu_p = mu_max * s / (ks + s)
    
    mu_eff = min(mu_p,model.mu)
    model =_data.get()
    model.l += dt * (mu_eff * _init_only_cell_lenghtening)
    model.mu += dt * (1.0 / tau_metabolism) * (mu_p - model.mu)
    model.cc = mu_eff *2.*1000. * 3.14 * (0.8e-6) * (0.8e-6) / 4.*model.l

    if(model.l>l_0):
        p.signal_division()


def division(_data,p:cpp_module.ParticleDataHolder):
    model =_data.get()
    child = model
    child.l/=2
    model.l/=2

    return model


def contribution(p,_data:cpp_module.OpaquePointer,contribs):
    model =_data.get()
   
    contribs[0,0]-=model.cc*p.weight
    return

def get_properties(_data):
    model =_data.get()
    return {"l":model.l,"mu":model.mu}

def mass(_data):
    model =_data.get()
    return 0.5*1000 * 3.14 * (0.8e-6) * (0.8e-6) / 4.*model.l


def show():
    pass

__start()

if __name__=="__main__":
    Warning("This file canno't be run as main script")