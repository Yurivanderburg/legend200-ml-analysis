"""
Figure: Charge trapping

This script creates a visualization of charge trapping
in a High Purity Germanium (HPGe) detector. 
"""

# --- Standard library ---
import numpy as np
import os

# --- Third party ---
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# --- Project modules ---
from utils import set_style



# ---------------------------------------------------------------------
# Figure-specific plotting function
# ---------------------------------------------------------------------
def charge_notrap(t_tot, t_e, q0):
    """
    Charge collection for a planar HPGe detector w/o charge trapping.
    """
    q_ind = []
    for t in t_tot:
        if t <= t_e:
            q = q0*(t/t_e)
        else:
            q = q0
        q_ind.append(q)
    
    return q_ind



def charge_trap(t_tot, t_e, q0, tau_T):
    """
    Charge collection for a planar HPGe detector with charge trapping, but
    no detrapping.
    """
    q_ind = []
    for t in t_tot:
        if t <= t_e:
            q = (q0*tau_T/t_e)*(1 - np.exp(-t/tau_T))
        else:
            q = (q0*tau_T/t_e)*(1 - np.exp(-t_e/tau_T))
        q_ind.append(q)
    return q_ind



def charge_detrap(t):
    """
    Charge collection for a planar HPGe detector with trapping and slow
    slow (exponential) detrapping (experimental).
    """

    def f(x, a, b, c, d, e, g):
        """Simple polynomial fit functions for smooth curve. """

        return a + b*x + c*x**2 + d*x**3 + e*x**4 + g*x**5

    x_ = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2])
    y_ = np.array([0, 400, 650, 800, 850, 900, 930, 960, 980, 990, 995, 999, 1000])

    popt, pcov = curve_fit(f=f, xdata=x_, ydata=y_)


    return f(t, *popt)





def plot_charge_trappings(ax=None, save_path=None):
    """
    Generate the charge trapping figure.

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
    t = np.arange(0, 1.2, 1e-3)
    q0 = 1000
    tau_T = 0.15
    t_e = 0.4
    tau_D = 1


    if ax is None:
        fig, ax = plt.subplots()
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=True,      # ticks along the bottom edge are off
        left=True,         # ticks along the top edge are off
        labelbottom=False, 
        labelleft=False) # labels along the bottom edge are off
    ax.set_yticks([0, 1000], labels=["0", "$q_0$"])
    ax.plot(t, charge_notrap(t, t_e, q0), label='no charge trapping')
    ax.plot(t, charge_trap(t, t_e, q0, tau_T), label='permanent trapping')
    ax.plot(t, charge_detrap(t), label='slow detrapping')
    ax.set_xlabel("Time [a.u.]")
    ax.set_ylabel(r"$Q_{\mathrm{ind}}$ [a.u.]")
    ax.legend()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1100)
    fig.tight_layout()

    if save_path:
        os.makedirs(save_path, exist_ok=True)
        plt.savefig(os.path.join(save_path, "Charge_trapping.png"), bbox_inches="tight")

    return ax

# ---------------------------------------------------------------------
# Main wrapper
# ---------------------------------------------------------------------

def main():
    """Run the script to generate and save the figure."""
    set_style()
    plot_charge_trappings(save_path="thesis_figures/output/03_legend/")


if __name__ == "__main__":
    main()