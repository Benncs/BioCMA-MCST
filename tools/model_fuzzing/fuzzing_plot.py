# import pickle
import matplotlib.pyplot as plt
# import numpy as np

# with open("pickled_data.pkl","rb") as file:
#     arr = pickle.load(file)


def add_fig(arr, c):
    print(c)
    # arr = arr.reshape(4,-1)
    plt.figure()
    plt.hist(
        arr[0, :],
        bins=100,
        edgecolor="black",
        alpha=0.5,
        color="blue",
    )
    plt.hist(arr[1, :], bins=100, edgecolor="black", alpha=0.5, label="child")
    plt.legend()


def show():
    plt.show()
