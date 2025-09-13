"""
Figure: Attenuation of Germanium

This script creates the typical photon attentuation figure
for Germanium, based on NIST data.  
"""

# --- Standard library ---
import os

# --- Third party ---
import pandas as pd
import matplotlib.pyplot as plt

# --- Project modules ---
from thesis_figures.utils import set_style



# ---------------------------------------------------------------------
# Figure-specific plotting function
# ---------------------------------------------------------------------
def plot_attenuation(ax=None, save_path=None):
    """
    Generate the photon attenuation figure.

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

    # Load NIST data:
    df = pd.read_csv("thesis_figures/utils/attenuation_ge.csv")


    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(df.Photon_Energy, df.Total_With_Coh, label="Total", linewidth=2)
    ax.plot(df.Photon_Energy, df.Incoherent, label="Compton")
    ax.plot(df.Photon_Energy, df.Photoelectric, label="Photoelectric effect")
    ax.plot(df.Photon_Energy, (df.PP_Nuclear + df.PP_Electron) , label="Pair production")

    ax.axvline(x=2.039, ls='--', color='C4', label=r"$Q_{\beta \beta} = 2039$ keV")

    # plt.axvline(x = 1.11031E-02, ymin=1e-3, ymax=1e4, color="C1", ls='--', label="K")
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlim(df.Photon_Energy.min(), df.Photon_Energy.max()) # np.min(energy), np.max(energy))

    ax.set_xlabel("Energy [MeV]")
    ax.set_ylabel(r"Attentuation $\mu / \rho$ [cm$^2$/g]")
    ax.set_xlim(1e-3, 1e2)
    ax.set_ylim(1e-3, 1e4)
    fig.tight_layout()

    fig.legend()

    if save_path:
        os.makedirs(save_path, exist_ok=True)
        plt.savefig(os.path.join(save_path, "Attenuation_ge.png"), bbox_inches="tight")

    return ax

# ---------------------------------------------------------------------
# Main wrapper
# ---------------------------------------------------------------------

def main():
    """Run the script to generate and save the figure."""
    set_style()
    plot_attenuation(save_path="thesis_figures/output/03_legend/")


if __name__ == "__main__":
    main()