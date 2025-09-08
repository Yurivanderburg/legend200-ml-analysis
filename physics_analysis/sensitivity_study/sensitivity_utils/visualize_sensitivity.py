"""
Plot utils for Bayesian sensitivity study.
"""


# --- Standard library ---
import os

# --- Third party library ---
import numpy as np
import matplotlib.pyplot as plt


# --- Project modules ---



def plot_background_index(data: dict) -> None:
    """
    Plot scatter plot and histogram of best-fit background index.

    Parameters
    ----------
    data : dict
        Dict of dict of summary
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), sharey="row", gridspec_kw={'wspace': 0.25}, constrained_layout=True)
    for idx, key in enumerate(data.keys()):
        # Scatter plot on top row
        x = np.arange(1, len(data[key]['global_mode_B']) + 1)
        axes[0, idx].scatter(x, np.array(data[key]['global_mode_B']), s=2)
        axes[0, idx].set_xlabel('Index')
        axes[0, 0].set_ylabel(r'Best-fit $\mathcal{B}$ [$10^{-4}$ counts/(keV$\cdot$kg$\cdot$yr)]')
        if 'AoE' in key:
            axes[0, idx].set_title(f'A/E method')
        else:
            axes[0, idx].set_title(f'Transformer method')


        # Histogram on bottom row
        axes[1, idx].hist(data[key]['global_mode_B'], bins=100)
        axes[1, idx].set_xlabel(r'Best-fit $\mathcal{B}$ [$10^{-4}$ counts/(keV$\cdot$kg$\cdot$yr)]')
        axes[1, 0].set_ylabel('Counts')
        axes[1, idx].set_yscale('log')

        axes[0, idx].grid(True, linestyle='--', alpha=0.5)
        axes[1, idx].grid(True, linestyle='--', alpha=0.5)

    plt.savefig("Plots/SensStudy/Plot_bkg_fit_per_toy.png", dpi=400)

    return None



def plot_signal_half_rate(data: dict) -> None:
    """
    Plot scatter plot and histogram of best-fit background index.

    Parameters
    ----------
    data : dict
        Dict of dict of summary
    """
    # Set up figure and axes (2x2 grid)
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), sharey="row", gridspec_kw={'wspace': 0.1}, constrained_layout=True)

    for idx, key in enumerate(data.keys()):
        x = np.arange(1, len(data[key]['global_mode_S']) + 1)

        # Scatter plot on top row
        axes[0, idx].scatter(x, np.array(data[key]['global_mode_S'])*1e27, s=2)
        axes[0, idx].set_xlabel('Index')
        axes[0, 0].set_ylabel(r'Best-fit $\mathcal{S}$ [$10^{-27}$ yr$^{-1}$]')
        if 'AoE' in key:
            axes[0, idx].set_title(f'A/E method')
        else:
            axes[0, idx].set_title(f'Transformer method')

        # Histogram on bottom row
        axes[1, idx].hist(data[key]['global_mode_S'], bins=100)
        axes[1, idx].set_xlabel(r'Best-fit $\mathcal{S}$ [yr$^{-1}$]')
        axes[1, 0].set_ylabel('Counts')
        axes[1, idx].set_xlim(-0.05e-27, 1.25e-27)
        axes[1, idx].set_yscale('log')

        axes[0, idx].grid(True, linestyle='--', alpha=0.5)
        axes[1, idx].grid(True, linestyle='--', alpha=0.5)

    plt.savefig("Plots/SensStudy/Plot_signal_fit_per_toy.png", dpi=400)
    plt.close()
    
    return None


def plot_roi_counts(data: dict, method=str) -> None:
    """
    Plot a overview of how # events in Qbb ROI effect the 
    best-fit signal half rate signal. 

    Parameters
    ----------
    data : dict
        Dict of summary (inner!)
    method : str
        Method (A/E or Transformer)
    """
    # Categorize roi counts into groups
    s_groups = {}
    for i in data['roi_count']:
        s_groups[i] = []
    for s, cat in zip(data['global_mode_S'], data['roi_count']):
        s_groups[cat].append(s)
    s_groups['All'] = data['global_mode_S']    

    # Bin edges based on full distribution
    bins = np.histogram_bin_edges(data['global_mode_S'], bins=50)

    # Define colors (colorblind-friendly)
    colors = {
        0: '#377eb8',    # blue
        1: '#984ea3',    # purple
        2: '#e41a1c',    # red
        3: '#4daf4a',   # green
        'All': "#ff9900"
    }

    # --- Plot ---
    plt.figure(figsize=(10, 6))
    for label in s_groups.keys():
        plt.hist(
            s_groups[label],
            bins=bins,
            alpha=0.6 if label != 'All' else 0.2,
            label=rf"Toys w. {label} cts. in $Q_{{\beta \beta}} \pm $ FWHM" if label != 'All' else 'All',
            color=colors[label],
            histtype='step' if label != 'All' else 'stepfilled',
            linewidth=1.7
        )

    # Axis labels and formatting
    plt.xlabel(r"Best-fit $\mathcal{S}$ [yr$^{-1}$]", fontsize=14)
    plt.ylabel("Counts", fontsize=14)
    plt.yscale('log')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    plt.legend(
        fontsize=12,
        title_fontsize=13,
        loc='upper right'
    )
    plt.title(r"Distribution of Best-fit $\mathcal{S}$" + f" ({method} method)", fontsize=16)
    plt.tight_layout()

    plt.savefig(f"Plots/SensStudy/Plot_sig_per_roi_count_{method}.png", dpi=400)
    plt.close()

    return None


def plot_sensitivity(data: dict) -> None:
    """
    Plot comparision of A/E and Transformer for sensitivity

    Parameters
    ----------
    data : dict
        Dict of dict of summary
    """
    fig, ax = plt.subplots(len(data.keys()), 1, sharex=True, sharey=True, figsize=(12,6))
    for i, key in enumerate(data.keys()):
 
        ax[i].hist(data[key]['S90'], bins=100, label=r"$\mathcal{S}_{90}$", alpha=0.9)
        ax[i].axvspan(data[key]['summary']['S90']['low'], 
                    data[key]['summary']['S90']['high'], 
                    color='red', alpha=0.2, label=r"68% CI")
        ax[i].axvline(data[key]['summary']['S90']['med'], 
                    color='black', label=r"Median $\mathcal{S}_{90}$")
        ax[i].set_ylabel("Counts")
        ax[i].set_yscale("log")

        if 'AoE' in key:
            ax[i].text(0.9, 0.825, "A/E", transform=ax[0].transAxes, fontsize=18)
        else:
            ax[1].text(0.8, 0.825, "Transformer 1", transform=ax[1].transAxes, fontsize=18)

    ax[1].set_xlabel(r"Signal half-rate [yr$^{-1}$]")
    ax[1].set_xlim(0.5e-27, 2.5e-27)

    handles, labels = ax[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3)
    plt.tight_layout()

    plt.savefig("Plots/SensStudy/Plot_results_S.png", dpi=400, bbox_inches='tight')
    plt.close()