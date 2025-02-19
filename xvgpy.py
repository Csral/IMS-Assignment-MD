import numpy as np
import matplotlib.pyplot as plt
import os

def read_xvg(file_path):
    """
    Reads an .xvg file and returns the data as a NumPy array.
    Comments (lines starting with # or @) are ignored.
    """
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(('#', '@')):
                continue  # Skip comments and metadata
            data.append(list(map(float, line.split())))
    return np.array(data)

def plot_xvg(file_path, xlabel="X-axis", ylabel="Y-axis", title="Plot"):
    """
    Plots the data from an .xvg file.
    """
    data = read_xvg(file_path)
    x = data[:, 0]  # First column is usually the X-axis
    y = data[:, 1:] # Remaining columns are Y-values

    plt.figure(figsize=(8, 6))
    for i in range(y.shape[1]):
        plt.plot(x, y[:, i], label=f"Column {i+2}")

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.show()

os.chdir("build")
# Example usage
file_path = "rmsd.xvg"  # Replace with your .xvg file path
plot_xvg(file_path, xlabel="Time (ps)", ylabel="Value", title="RMSD")