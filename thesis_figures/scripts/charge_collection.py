"""
Figure: Electron-hole collection process in HPGe

This script creates a visualization of electron and hole 
collection process in a HPGe semiconductor.
"""

# --- Standard library ---
import numpy as np
import os

# --- Third party ---
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# --- Project modules ---
from thesis_figures.utils import set_style



# ---------------------------------------------------------------------
# Figure-specific plotting function
# ---------------------------------------------------------------------
def charge_collection(t: np.ndarray, x:np.ndarray):
    """
    Total induced charge on the cathodes of a planar HPGe detector.

    Parameters
    ----------
    t : np.ndarray
    x : np.ndarray
    """

    # Constants
    d = 20          # Detector thickness [cm]
    q0 = 1.0         # Total induced charge [arb. units]
    v_e = 1e7        # Electron drift velocity [cm/s]
    v_h = 1e7      # Hole drift velocity [cm/s]

    # Compute collection times
    t_e = x / v_e           # Electron collection time [s]
    t_h = (d - x) / v_h     # Hole collection time [s]

    Q_vals = np.zeros_like(t)

    # Case 1: both electrons and holes drifting (Eq. 12.20a)
    mask1 = (t < t_e) & (t < t_h)
    Q_vals[mask1] = q0 * ((v_e / d) * t[mask1] + (v_h / d) * t[mask1])

    # Case 2: electrons collected, holes still drifting (Eq. 12.20b)
    mask2 = (t >= t_e) & (t < t_h)
    Q_vals[mask2] = q0 * (x / d + (v_h / d) * t[mask2])

    # Case 3: holes collected, electrons still drifting (Eq. 12.20c)
    mask3 = (t >= t_h) & (t < t_e)
    Q_vals[mask3] = q0 * ((v_e / d) * t[mask3] + (d - x) / d)

    # Case 4: both collected
    Q_vals[t >= max(t_e, t_h)] = q0

    return Q_vals



def plot_charge_collection(ax=None, save_path=None):
    """
    Generate the charge collection figure.

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

    # Parameters
    t = np.linspace(0, 3e-6, 1000)
    d = 20          # Detector thickness [cm]

    if ax is None:
        fig, axes = plt.subplots(3, 1, figsize=(6, 7), sharex=True)

    axes[0].plot(t * 1e6, charge_collection(t, x=0.15*d), color='black')  # Convert time to μs for readability
    axes[0].annotate('Electrons collected', xy=(0.3, 0.29), xytext=(0.2, 0.1),
                    arrowprops=dict(arrowstyle='->'), fontsize=10)
    axes[0].annotate('Holes collected', xy=(1.7, 0.99), xytext=(1.6, 0.8),
                    arrowprops=dict(arrowstyle='->'), fontsize=10)
    axes[1].plot(t * 1e6, charge_collection(t, x=0.5*d), color='black')  # Convert time to μs for readability
    axes[2].plot(t * 1e6, charge_collection(t, x=0.85*d), color='black')  # Convert time to μs for readability
    axes[2].annotate('Holes collected', xy=(0.3, 0.29), xytext=(0.2, 0.1),
                    arrowprops=dict(arrowstyle='->'), fontsize=10)
    axes[2].annotate('Electrons collected', xy=(1.7, 0.99), xytext=(1.6, 0.8),
                    arrowprops=dict(arrowstyle='->'), fontsize=10)

    for i in range(3):
        axes[i].grid(True)
        axes[i].set_ylabel("Q(t) [a.u.]")

    plt.xlabel(r"Time [$\mu$s]")
    plt.tight_layout()

    if save_path:
        os.makedirs(save_path, exist_ok=True)
        plt.savefig(os.path.join(save_path, "Charge_collection.png"), bbox_inches="tight")

    return ax



# ---------------------------------------------------------------------
# Main wrapper
# ---------------------------------------------------------------------

def main():
    """Run the script to generate and save the figure."""
    set_style()
    plot_charge_collection(save_path="thesis_figures/output/03_legend/")


if __name__ == "__main__":
    main()