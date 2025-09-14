"""
Test Transformer and A/E performance on 2vbb events.

This script creates a figure of the 2vbb spectrum, as measured 
by the LEGEND experiment's HPGe detectors, tested on both the A/E
and the Transformer models.  """

# --- Standard library ---
import numpy as np
import os
import sys

# --- Third party ---
import matplotlib.pyplot as plt

# --- Project modules ---
from utils import set_style
from utils.data_io import load_pickle_file


# Parameters
DATA_PATH = "path/to/labelled_dataset.pkl"


# ---------------------------------------------------------------------
# Figure-specific plotting function
# ---------------------------------------------------------------------
def plot_calibration_spectrum(ax=None, save_path=None):
    """
    Generate the 2vbb energy spectrum figure.

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

    data = load_pickle_file(DATA_PATH)
    data_qc = data[ # Apply data quality cuts
        data.geds.psd.is_good &
        ~data.coincident.muon & 
        ~data.coincident.muon_offline & 
        ~data.coincident.spms 
    ]       


    if ax is None:
        fig, ax = plt.subplots(2, 1, figsize=(14, 8), sharex=False)

    # Binning and constants
    x_min, x_max = 0, 5000
    bin_keV = 10
    bin_edges = np.arange(x_min, x_max + bin_keV, bin_keV)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    Qbb = 2039
    Qbb_range = 20

    # --- Subplot 1: Full spectrum and mult. cuts ---
    # All events
    n_before, hist_bins, _ = ax[0].hist(
        data.geds.energy, bins=bin_centers, histtype='step', linewidth=1.2,
        color='C0', label="All events"
    )
    # After mult. cuts
    data_mult = data[~data.coincident.muon & 
                ~data.coincident.muon_offline & 
                ~data.coincident.spms 
                ]
    n_mult, _, _ = ax[0].hist(
        data_mult.geds.energy, bins=bin_centers, histtype='step', linewidth=1.2,
        color='C1', label="After coincidence cuts"
    )
    # After A/E 
    data_aoe = data_mult[data_mult.geds.psd.is_bb_like == True]
    n_AoE, _, _ = ax[0].hist(
        data_aoe.geds.energy, bins=bin_centers, histtype='step', linewidth=1.2,
        color='C2', label="A/E method"
    )

    # --- Subplot 2: A/E vs Transformer ---
    ax[1].hist(
        data_aoe.geds.energy, bins=bin_centers, histtype='step', linewidth=1.2,
        color='C2') #, label="A/E method"
    for i, model in enumerate(['Transformer_v1', 'Transformer_v2']):
        label = f"Transformer {model[-1]}"
        data_sub = data_mult[data_mult['Transformer_v' + str(model[-1])]['Prediction'] == 'prob_SSE']
        n_T, _, _ = ax[1].hist(
            data_sub.geds.energy, bins=bin_centers, histtype='step', linewidth=1.2,
            label=label, color=f"C{i+3}"
        )

    # Highlighted regions
    for a in ax:
        a.axvspan(Qbb - Qbb_range, Qbb + Qbb_range, color='green', alpha=0.2)
        a.axvspan(1000, 1300, color='blue', alpha=0.1)

    # Annotations
    ax[0].annotate(r"$Q_{\beta \beta}$ window", xy=(Qbb + 10, max(n_before)*0.5), color='green')
    ax[0].annotate(r"PSD extr. window", xy=(1000, max(n_before)*0.5), color='blue', alpha=0.75)
    ax[0].set_yscale('log')
    ax[0].set_ylabel(f"Counts / {bin_keV}" + r"$\,\mathrm{{keV}}$")
    ax[1].set_ylabel(f"Counts / {bin_keV}" + r"$\,\mathrm{{keV}}$")
    ax[1].set_xlim(0, 2600)
    ax[1].set_yscale('log')
    ax[1].set_xlabel("Energy [keV]")

    fig.legend(loc='upper center', bbox_to_anchor=(0.5, 0), ncol=3)
    fig.tight_layout()

    if save_path:
        os.makedirs(save_path, exist_ok=True)
        plt.savefig(os.path.join(save_path, "SSE_cuts_2vbb_presentation.png"), bbox_inches="tight")

    return ax

# ---------------------------------------------------------------------
# Main wrapper
# ---------------------------------------------------------------------

def main():
    """Run the script to generate and save the figure."""

    if not os.path.exists(DATA_PATH):
        sys.exit("Dataset not available. Please provide 'calibration_data.lh5' to run this script.")

    set_style()
    plot_calibration_spectrum(save_path="thesis_figures/output/05_psd/")


if __name__ == "__main__":
    main()