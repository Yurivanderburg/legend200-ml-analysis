"""
Figure: 2vbb and 0vbb energy spectra
spectra

This script creates the energy spectra of the two processes in an 
experiment, where only the electron energies are measured.  
"""

# --- Standard library ---
import os
import numpy as np

# --- Third party ---
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.integrate import simps


# --- Project modules ---
from utils import set_style



# ---------------------------------------------------------------------
# Figure-specific plotting function
# ---------------------------------------------------------------------
def plot_energyspectra(ax=None, save_path=None):
    """
    Generate the energy spectra figure.

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

    # Load the Kotila–Iachello summed spectrum
    df = pd.read_csv('thesis_figures/utils/ge76_sums.csv')

    E_keV = df.E_MeV * 1000  # convert to keV

    # Normalize 2νββ shape
    dGamma = df.dGamma / simps(df.dGamma, E_keV)

    # Constants
    N_A = 6.022e23        # Avogadro
    M_Ge76 = 75.921       # g/mol
    frac = 0.88            # assume 100% enrichment
    T12_2nu = 1.92e21     # yr
    T12_0nu = 1e26        # yr
    Q_bb = 2039  # keV
    FWHM = 2.5   # detector resolution (keV)

    atoms_per_kg = (N_A * frac) / (M_Ge76 / 1000)

    # Scale to physical rate
    rate_2nu = (atoms_per_kg / T12_2nu) * dGamma  # events / keV / kg / yr
    rate_2nu = gaussian_filter1d(rate_2nu, sigma=1)  # simulate detector smearing


    sigma = FWHM / 2.355

    E_2 = np.arange(2000, 2100, 0.1)

    rate_0nu = np.exp(-0.5 * ((E_2 - Q_bb) / sigma) ** 2)
    rate_0nu /= simps(rate_0nu, E_2)
    rate_0nu *= (atoms_per_kg / T12_0nu)


    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(E_keV, rate_2nu, label=r"$2 \nu \beta \beta$", color="C2", linewidth=1.3)
    ax.plot(E_2, rate_0nu, label=r"$0 \nu \beta \beta$", linewidth=1.3)
    ax.set_xlabel('Energy (keV)')
    ax.set_ylabel(r"counts/(keV$\cdot$kg$\cdot$yr)")
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 10)
    ax.set_xlim(800, 2300)
    plt.legend()


    if save_path:
        os.makedirs(save_path, exist_ok=True)
        plt.savefig(os.path.join(save_path, "DoubleBeta_spectrum.png"), bbox_inches="tight")

    return ax

# ---------------------------------------------------------------------
# Main wrapper
# ---------------------------------------------------------------------

def main():
    """Run the script to generate and save the figure."""
    set_style()
    plot_energyspectra(save_path="thesis_figures/output/02_neutrino/")


if __name__ == "__main__":
    main()