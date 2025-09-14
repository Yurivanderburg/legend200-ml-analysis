"""
Figure: Waveform example

This script creates a few waveform examples showing
what we expect for which label.
"""

# --- Standard library ---
import numpy as np
import os
import json

# --- Third party ---
import matplotlib.pyplot as plt

# --- Project modules ---
from utils import set_style



# ---------------------------------------------------------------------
# Figure-specific plotting function
# ---------------------------------------------------------------------
def normalize_signal(waveform: np.ndarray, grad: bool = False):
    """ 
    Returns normalized and scaled charge and current.
    
    Parameters
    ----------
    waveform : np.ndarray
        Waveform to be normalized or derived
    grad : bool
        Whether or not to calculate and return the derivative
    """
    wf_norm = (waveform - np.min(waveform))/(np.max(waveform) - np.min(waveform))
    
    if grad:
        gradient = np.gradient(wf_norm)
        
        return 5*gradient
    else:
        return wf_norm



def plot_waveform_example(ax=None, save_path=None):
    """
    Generate the waveform example figure.

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

    with open("thesis_figures/utils/waveform_example.json", "r") as f:
        data = json.load(f)
    x = np.linspace(0, 1399, 1400)*16e-3


    if ax is None:
        fig, ax = plt.subplots(1,4, figsize=(20,6), sharex=True, sharey=True)

    # SSE Plot
    ax[0].plot(x, normalize_signal(data['SSE'], grad=False), label="Charge pulse", color='C0', linewidth=1.2)
    ax[0].plot(x, normalize_signal(data['SSE'], grad=True), label="Current pulse", color='C1', linewidth=1.2)

    acceptance_min = np.max(normalize_signal(data['SSE'], grad=True)) - 0.025
    acceptance_max = np.max(normalize_signal(data['SSE'], grad=True)) + 0.025


    # MSE Plot
    ax[1].plot(x, normalize_signal(data['MSE'], grad=False), color='C0', linewidth=1.2)
    ax[1].plot(x, normalize_signal(data['MSE'], grad=True), color='C1', linewidth=1.2)

    # ncontact Plot
    ax[2].plot(x, normalize_signal(data['nc'], grad=False), color='C0', linewidth=1.2)
    ax[2].plot(x, normalize_signal(data['nc'], grad=True), color='C1', linewidth=1.2)

    # pcontact Plot
    ax[3].plot(x, normalize_signal(data['pc'], grad=False), color='C0', linewidth=1.2)
    ax[3].plot(x, normalize_signal(data['pc'], grad=True), color='C1', linewidth=1.2)


    ax[0].set_title("SSE")
    ax[1].set_title("MSE")
    ax[2].set_title("n-contact")
    ax[3].set_title("p-contact")


    for i in range(4):
        ax[i].set_xlabel("Time [μs]")
        if i == 0:
            ax[i].axhspan(acceptance_min, acceptance_max, color='skyblue', alpha=0.4, label="Acceptance window")
        else:        
            ax[i].axhspan(acceptance_min, acceptance_max, color='skyblue', alpha=0.4)



    ax[0].set_ylabel("norm. signal")

    plt.xlim(5, 8.5)
    plt.tight_layout()
    fig.legend(loc='center', bbox_to_anchor=(0.5, -0.1), ncols=3)
    plt.show()
    if save_path:
        os.makedirs(save_path, exist_ok=True)
        plt.savefig(os.path.join(save_path, "Waveform_example.png"), bbox_inches="tight")

    return ax

# ---------------------------------------------------------------------
# Main wrapper
# ---------------------------------------------------------------------

def main():
    """Run the script to generate and save the figure."""
    set_style()
    plot_waveform_example(save_path="thesis_figures/output/03_legend/")


if __name__ == "__main__":
    main()