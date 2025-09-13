"""
Figure: Activation functions

This script creates a visualization of different
activation functions in ML. 
"""

# --- Standard library ---
import numpy as np
import os

# --- Third party ---
from scipy.special import erf
import matplotlib.pyplot as plt

# --- Project modules ---
from thesis_figures.utils import set_style



# ---------------------------------------------------------------------
# Figure-specific plotting function
# ---------------------------------------------------------------------
def relu(x: np.ndarray) -> np.ndarray:
    """
    REctified Linear Unit function.
    """
    relu = np.zeros(len(x))

    for i in range(len(x)):
        relu[i] = max(x[i], 0)
        
    return relu



def gelu(x: np.ndarray) -> np.ndarray:
    """
    Gaussian Error Linear Unit function. 
    """
    return (x/2)*(1 + erf(x/np.sqrt(2)))



def sigmoid(x: np.ndarray) -> np.ndarray:
    """
    Sigmoid activation function.
    """
    return 1/(1 + np.exp(-x))



def plot_activation_function(ax=None, save_path=None):
    """
    Generate the activation functions figure.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Axis to plot on. If None, a new figure and axis are created.
    save_path : str, optional
        File path to save the figure. If None, the figure is not saved.

    Returns
    -------
    ax : matplotlib.axes.Axes
        Axis containing the plot.
    """

    x = np.linspace(-5, 5, 101)

    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(x, sigmoid(x), label="sigmoid")
    ax.plot(x, relu(x), label="ReLU")
    ax.plot(x, gelu(x), label="GELU")
    fig.legend()

    if save_path:
        os.makedirs(save_path, exist_ok=True)
        plt.savefig(os.path.join(save_path, "Activation_function.png"), bbox_inches="tight")

    return ax

# ---------------------------------------------------------------------
# Main wrapper
# ---------------------------------------------------------------------

def main():
    """Run the script to generate and save the figure."""
    set_style()
    plot_activation_function(save_path="thesis_figures/output/04_transformer/")


if __name__ == "__main__":
    main()