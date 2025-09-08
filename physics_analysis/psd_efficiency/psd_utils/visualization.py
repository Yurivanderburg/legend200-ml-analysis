"""
Plot utils for PSD efficiency at $Q_{\beta \beta}$ analysis.
"""


# --- Standard library ---
import os

# --- Third party library ---
import numpy as np
import matplotlib.pyplot as plt


# --- Project modules ---
from utils.math import linear, gauss_bkg
from .efficiency import mean_and_error_with_delta

fig, ax = plt.subplots(1,3, sharex=True, sharey=True, figsize=(16,6))
x = np.arange(1000, 2500, 1)


def plot_total_eff(
        config: dict, 
        results: dict
    ) -> None:
    """
    Plot the individual PSD efficiencies of 0vbb events, 
    the Th-228 DEP and the Co-56 DEP, including the linear 
    fit and uncertainties. 

    Parameters
    ----------
    config : dict
        Configuration file
    results : dict
        Summary file including the fit and individual efficiencies

    Returns
    -------
        None
    """

    out_dir = os.path.join(config['plot_dir'], 'PSD_eff')
    os.makedirs(out_dir, exist_ok=True)

    for i, model in enumerate(['AoE', 'T1', 'T2']):

        # Plot individual efficiencies as error bars
        ax[i].errorbar(1593.5, results[model]['Th228']['efficiency'], 
                    yerr=results[model]['Th228']['uncertainty'], 
                    fmt=".", ms=10, capsize=5, label="$^{208}$Tl DEP")
        ax[i].errorbar(2231.5, results[model]['Co56']['efficiency'],
                    results[model]['Co56']['uncertainty'],
                    fmt=".", ms=10, capsize=5, label=r"$^{56}$Co DEP")
        ax[i].errorbar(1150, results[model]['2vbb']['eff'], yerr=results[model]['2vbb']['eff_err'], 
                fmt=".", ms=10, capsize=5, label=r"$2 \nu \beta \beta$")


        fit_param = results[model]['fit']
        eff_qbb = results[model]['values']['qbb']
        eff_qbb_err = results[model]['errors']['qbb']

        if 'T' in model:
            ax[i].set_title(f"Transformer {model[-1]}")
        else:
            ax[i].set_title("A/E method")

        # Plot the linear fit
        if model != 'T2':
            ax[i].plot(x, linear(x, *fit_param), ls="--", color='red')
            # ax[i].axvline(x=2039, ymin=0, ymax=1, color='green', ls='--', 
            ax[i].errorbar(x = 2039, y = eff_qbb, yerr = eff_qbb_err, color='purple', fmt=".", ms=10, capsize=5, 
                label=r"$\epsilon_{Q_{\beta \beta}}$ = " + f"{np.round(eff_qbb, 3)} +/- {np.round(eff_qbb_err, 3)}")
            ax[i].set_xlim(1000, 2500)

                    
            h,l = ax[i].get_legend_handles_labels()

            ax[i].legend(h[-1:], l[-1:], loc='lower left', fontsize=16)

    handles, labels = ax[0].get_legend_handles_labels()
    # fig.suptitle(r"PSD efficiency at $Q_{\beta \beta}$")
    fig.supxlabel("Energy [keV]")
    ax[0].set_ylabel("Efficiency")


    fig.tight_layout()
    fig.legend(handles[:-1], labels[:-1], loc='upper center', bbox_to_anchor=(0.5, 0), ncol=3)
    # fig.subplots_adjust(bottom=0.2)  # Make room for legend
    if config['save_figures']:
        plt.savefig(os.path.join(out_dir, 'PSD_eff_qbb.png'), 
                    dpi=400, bbox_inches='tight')

    return None



def plot_per_detector_run(
        config: dict, 
        eff_summary: dict, 
        x_labels: list, 
        delta_name: str, 
        y_lims: tuple=(0.75, 1.02), 
        save_path: str=None
    )-> None:
    """
    Plots the PSD efficiency for different models across detectors / runs.

    Parameters
    ----------

    config : dict
        Configuration file
    eff_summary : dict
        Efficiency summary subset to plot (only ['summary'])
    x_labels : list
        X-labels of detector names
    delta-name : str 
        Either "det" or "run" - determines what to plot
    y_lims : tuple
        Y-axis limts for matplotlib
    save_path : str
        If not None, saves figure to file at said path

    Returns
    -------
    None
    """

    fig, ax = plt.subplots(2, len(config['models']), sharex=True, sharey=False, 
                    figsize=(17, 8), gridspec_kw={'height_ratios': [2.75, 1]})

    for i, model in enumerate(config['models']):
        # Unpack values
        mean_eff, std_eff = mean_and_error_with_delta(
                eps=eff_summary[model]['effs'],
                sig=eff_summary[model]['eff_errs'],
                delta=eff_summary[model]['fit_result'].x[1]
            )
    
        xmin_ = -0.5
        xmax_ = len(x_labels) - 0.5
        if delta_name == 'run':
            delta_label = r"$\langle \delta_{\mathrm{run}} \rangle$ = "
        elif delta_name == 'det':
            delta_label = r"$\langle \delta_{\mathrm{det}} \rangle$ = "
        else:
            delta_label = ""

        #Add standardized residuals
        residuals = (eff_summary[model]['effs'] - mean_eff)/(np.sqrt(eff_summary[model]['eff_errs']**2 + eff_summary[model]['fit_result'].x[1]**2))
        chi2 = np.sum(((eff_summary[model]['effs'] - mean_eff)**2)/(eff_summary[model]['eff_errs']**2 + eff_summary[model]['fit_result'].x[1]**2))


        # Upper plots - mean
        ax[0,i].errorbar(x_labels, eff_summary[model]['effs'], yerr=eff_summary[model]['eff_errs'], fmt=".", capsize=3, markersize=8)
        ax[0,i].hlines(y=mean_eff, xmin=xmin_, xmax=xmax_, linestyle='--', color='black', alpha=0.6,
                     label=rf"$\epsilon$ = {mean_eff:.3f} ± {std_eff:.3f}")
        ax[0,i].fill_between(x=[xmin_, xmax_],
                           y1=mean_eff - std_eff,
                           y2=mean_eff + std_eff,
                           color='C0', alpha=0.2)
        ax[0,i].fill_between(x=[xmin_, xmax_],
                           y1=mean_eff - eff_summary[model]['fit_result'].x[1],
                           y2=mean_eff + eff_summary[model]['fit_result'].x[1],
                           color='C1', alpha=0.2, label=delta_label + str(np.round(eff_summary[model]['fit_result'].x[1], 3)))
        ax[0,i].hlines([mean_eff + std_eff, mean_eff - std_eff], xmin=xmin_, xmax=xmax_, alpha=0.5)

        ax[0,i].set_ylim(*y_lims)
        title = f"Transformer {model[-1]}" if 'T' in model else "A/E method"
        ax[0,i].set_title(title)
        ax[0,i].legend(loc='upper left', fontsize=16)
        if i == 0:
            ax[0,i].set_ylabel("Efficiency")
            ax[1,i].set_ylabel("Norm. res.")


        # Lower plot - residuals 
        ax[1,i].errorbar(x_labels, residuals, yerr=1, fmt='.', ms=8, capsize=3,  
            color='gray', label="Normalized residuals")
        gof = chi2/(len(x_labels)-1)
        ax[1,i].text(0.02, 0.825, fr"$\chi^2/\mathrm{{ndf}} = {np.round(gof, 2)}$",
                    transform=ax[1, i].transAxes, fontsize=12)
        ax[1,i].axhline(0, color='black', linewidth=1)

   
        ax[1,i].set_ylim(-3.5, 3.5)
        ax[1,i].tick_params('x', labelrotation=90, labelsize=12)


    fig.supxlabel("Detectors")
    fig.tight_layout()

    if save_path:
        plt.savefig(os.path.join(config['plot_dir'], save_path), dpi=400, bbox_inches='tight')



def plot_peak_with_residuals(
        config: dict,
        efficiencies: dict,
        bin_range: tuple,
        plot_range: int, 
        res_range: int,
        bin_keV: int,
        peak_energy: float,
        ) -> None:
    """
    Plot the peak fits required for DEP PSD efficiency estimation.

    Parameters
    ----------
    config : dict
        Configuration file
    efficiencies : dict
        Efficiency dictionary containing parameters
    plot_range : int
        Range of the uppermost plot (the fit)
    res_range : int
        Range of the normalized residual plot (y-axis)
    bin_keV : int
        Binning of the histogram, in keV
    peak_energy : float
        Energy of the peak, in keV

    
    Returns
    -------
    None

    """
    # Initial settings
    out_dir = os.path.join(config['plot_dir'], 'PSD_eff')
    os.makedirs(out_dir, exist_ok=True)

    bin_edges = np.arange(bin_range[0], bin_range[1] + bin_keV, bin_keV) 
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2
    x_range = np.arange(peak_energy - 20, peak_energy + 20 + 0.01, 0.01)

    fig, axes = plt.subplots(nrows=3, ncols=len(config['models']), sharex='col', sharey='row', 
        figsize=(5 * len(config['models']), 9), gridspec_kw={'height_ratios': [3, 1, 1]})

    for i, model in enumerate(config['models']):
        ax_main = axes[0, i]
        ax_resid_p = axes[1, i]
        ax_resid_f = axes[2, i]

        # --- Main plot ---
        ax_main.bar(bin_centers, efficiencies['Total']['counts'], width=bin_keV, 
                    align='center', color='C0', alpha=0.5, label='Before cuts')
        ax_main.plot(x_range, gauss_bkg(x_range, *efficiencies['Total']['popt']), 
                    color='C0', linewidth=2)

        ax_main.bar(bin_centers, efficiencies[f"{model}_pass"]['counts'], width=bin_keV, 
                    align='center', color='green', alpha=0.5, label=f"After cuts (passed)")
        ax_main.plot(x_range, gauss_bkg(x_range, *efficiencies[f"{model}_pass"]['popt']), 
                    color='green', linewidth=2)

        ax_main.bar(bin_centers, efficiencies[f"{model}_fail"]['counts'], width=bin_keV, 
                    align='center', color='red', alpha=0.5, label=f"After cuts (failed)")
        ax_main.plot(x_range, gauss_bkg(x_range, *efficiencies[f"{model}_fail"]['popt']), 
                    color='red', linewidth=2)

        ax_main.set_xlim(peak_energy - plot_range, peak_energy + plot_range)
        ax_main.set_title(model)
        if i == 0:
            ax_main.set_ylabel(f"Counts / {bin_keV} keV")

        if 'T' in model:
            ax_main.set_title(f"Transformer {model[-1]}")
        else:
            ax_main.set_title("A/E method")

        # --- Residuals plot 1 ---
        observed_p = efficiencies[f"{model}_pass"]['counts']
        expected_p = gauss_bkg(bin_centers, *efficiencies[f"{model}_pass"]['popt'])
        residuals_p = (observed_p - expected_p)/np.sqrt(observed_p)

        ax_resid_p.errorbar(bin_centers, residuals_p, yerr=1, fmt='.', ms=8, capsize=3,  
                color='gray', label="Normalized residuals")
        ax_resid_p.axhline(0, color='black', linewidth=1)
        if i == 0:
            ax_resid_p.set_ylabel("Res. (passed)")
        ax_resid_p.set_ylim(-res_range, +res_range)

        # --- Residuals plot 2 ---
        observed_f = efficiencies[f"{model}_fail"]['counts']
        expected_f = gauss_bkg(bin_centers, *efficiencies[f"{model}_fail"]['popt'])
        residuals_f = (observed_f - expected_f)/(np.sqrt(observed_f))

        ax_resid_f.errorbar(bin_centers, residuals_f, yerr=1, fmt='.', ms=8, capsize=3,  
                color='gray', label="Normalized residuals")
        ax_resid_f.axhline(0, color='black', linewidth=1)
        if i == 0:
            ax_resid_f.set_ylabel("Res. (failed)")
        ax_resid_f.set_ylim(-res_range, +res_range)


    # fig.suptitle(r"$^{208}$Tl DEP peak fit (detector " + str(detector) + ") for events before and after PSD cut")
    fig.supxlabel("Energy [keV]")


    # Collect handles and labels from one subplot (e.g., the first main axis)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    handles_, labels_ = axes[1, 0].get_legend_handles_labels()

    handles += handles_
    labels += labels_

    # Add a single shared legend below the entire figure
    fig.legend(
        handles, labels,
        loc='lower center',
        bbox_to_anchor=(0.5, -0.075),
        ncol=len(handles)
    )
    fig.subplots_adjust(bottom=0.5)  # Make room for legend

    # Improve layout
    fig.tight_layout()
    if config['save_figures']:
        plt.savefig(os.path.join(out_dir, f"Peakfit_Th228.png"), 
                                 dpi=400, bbox_inches='tight')
