"""
Calibration spectrum figure.

This script creates a Th-228 calibration spectrum, as measured
by the LEGEND experiment's HPGe detectors. """

# --- Standard library ---
import numpy as np
import os
import sys

# --- Third party ---
import matplotlib.pyplot as plt

# --- Project modules ---
from utils import set_style
from utils.data_io import load_lh5_file


# Parameters
DATA_PATH = "path/to/calibration_data.lh5"
CHANNEL_ID = 'ch1116803'


# ---------------------------------------------------------------------
# Figure-specific plotting function
# ---------------------------------------------------------------------
def plot_calibration_spectrum(ax=None, save_path=None):
    """
    Generate the calibration spectrum figure.

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

    data = load_lh5_file(DATA_PATH, channel_id=CHANNEL_ID)

    fep_2614 = 2614.5  # Full energy peak from 208Tl decay
    dep_1592 = 1592.5  # Double escape peak (SSEs)
    fep_1620 = 1620.5  # FEP from 212Bi decay
    sep_2103 = 2103.5  # Approximate position of SEP
    qbb = 2039 
    cc_region_low = qbb - 35
    cc_region_high = qbb + 35
    width = 10 # +/- in keV
    bin_keV = 3
    binning_keV = np.arange(-200, 3500 + bin_keV, bin_keV)


    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(14, 5))

    ax.hist(data.cuspEmax_ctc_cal, 
            bins=binning_keV, 
            color='C0', 
            fill=None, 
            histtype='step', 
            label='Calibration data', 
            alpha=1, 
            linewidth=1.5
            )
    ax.set_yscale('log')
    ax.set_ylabel(f"Counts / {bin_keV} keV")
    ax.set_xlabel('Energy [keV]')

    ax.axvspan(xmin=fep_2614 - width, xmax=fep_2614 + width, 
               color='purple', label=f"FEP (2614.5 ± {width} keV)", alpha=0.3)
    ax.axvspan(xmin=dep_1592 - width, xmax=dep_1592 + width, 
               color='skyblue', label=f"DEP (1592.5 ± {width} keV, SSE)", alpha=0.3)
    ax.axvspan(xmin=fep_1620 - width, xmax=fep_1620 + width, 
               color='orange', label=f"FEP (1620.5 ± {width} keV, MSE)",alpha=0.3)
    ax.axvspan(xmin=sep_2103 - width, xmax=sep_2103 + width, 
               color='red', label=f"SEP (2103.5 ± {width} keV, MSE)", alpha=0.3)
    ax.axvspan(cc_region_low, cc_region_high, 
               color='green', alpha=0.3, label="Compton Continuum ($Q_{ββ}$ ± 35 keV)")

    ax.set_xlim(-200, 3500)
    fig.legend(loc='center', bbox_to_anchor=(0.5, -0.1), ncols=3)
    plt.tight_layout()

    if save_path:
        os.makedirs(save_path, exist_ok=True)
        plt.savefig(os.path.join(save_path, "Calibration_spectrum.png"), bbox_inches="tight")

    return ax

# ---------------------------------------------------------------------
# Main wrapper
# ---------------------------------------------------------------------

def main():
    """Run the script to generate and save the figure."""

    if not os.path.exists(DATA_PATH):
        sys.exit("Dataset not available. Please provide 'calibration_data.lh5' to run this script.")

    set_style()
    plot_calibration_spectrum(save_path="thesis_figures/output/03_legend/")


if __name__ == "__main__":
    main()