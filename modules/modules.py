import time 
import cpp_module



a = cpp_module.a


def init(initmass:float):
    _data={}
    _data['mass']=initmass
    print("initialised")

def update(value):   
    print(value)
    # time.sleep(0.2)
    # _data['mass']+=1


def show():
    print(a)

print("Module Loaded")