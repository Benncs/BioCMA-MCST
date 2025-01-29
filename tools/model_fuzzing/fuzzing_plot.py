import pickle 
import matplotlib.pyplot as plt 
import numpy as np 

# with open("pickled_data.pkl","rb") as file:
#     arr = pickle.load(file)


def add_fig(arr,c):
    print(c)
    # arr = arr.reshape(4,-1)
    plt.figure()
    plt.hist(arr[0,:],bins=100,edgecolor="black",
        alpha=0.7,
        color="blue",)


def show():
    plt.show()
