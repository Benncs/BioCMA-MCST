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




class Model:
    def __init__(self,a):
        self.a =a

def init(p:cpp_module.ParticleDataHolder):
    p.id = random.randint(0,10000)
    return Model(p.id)

def update(_data:cpp_module.OpaquePointer,p:cpp_module.ParticleDataHolder,concentrations):   
    model =_data.get()
    model.a+=1
    # print(concentrations)
    # print(model.a)
    # time.sleep(0.2)
    # _data['mass']+=1

def division(p:cpp_module.ParticleDataHolder):
    p.id = random.randint(0,10000)
    return Model(p.id)


def contribution(p,_data:cpp_module.OpaquePointer,contribs):
    contribs[0,0]+=10./p.weight
    contribs[:,1]+=2.
    pass


def show():
    pass

__start()

if __name__=="__main__":
    Warning("This file canno't be run as main script")