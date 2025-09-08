"""
Determine the PSD efficiency at Qbb, based on 2vbb events, 
Co-56 DEP events and Th-228 DEP events.  
"""

# --- Standard library ---
import datetime as dt
import os


# --- Project modules ---
from physics_analysis.psd_efficiency.psd_utils import psd_eff_2vbb
from physics_analysis.psd_efficiency.psd_utils import psd_eff_dep
from physics_analysis.psd_efficiency.psd_utils import psd_eff_timevar
from physics_analysis.psd_efficiency.psd_utils import psd_eff_qbb
from physics_analysis.psd_efficiency.config import config
from physics_analysis.psd_efficiency.psd_utils import visualization
from utils import log_utils
from utils.data_io import load_pickle_file


def main():
    """
    Calculate the total efficiency at $Q_{\beta \beta}$ and 
    plot some metrics. 
    """

    # ---- Logger ----
    now = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    logger = log_utils.setup_logger(
        filename=f"PSD_eff_{now}.log",
        archive_logs=True,
        log_folder="logs/psd/",
    )
    logger.info("Calculate PSD efficiencies...")


    # Load files if present:
    if not config['calculate_efficiencies']:
        eff_2vbb = load_pickle_file(os.path.join(config['global_path_out'], "Results_2vbb.pkl"))
        eff_Th228 = load_pickle_file(os.path.join(config['global_path_out'], "Results_DEP_Th228.pkl"))
        eff_Co56 = load_pickle_file(os.path.join(config['global_path_out'], "Results_DEP_Co56.pkl"))
        eff_timevar = load_pickle_file(os.path.join(config['global_path_out'], "Results_timevar.pkl"))
    
    else:

        # 1. PSD eff @ 2vbb events
        eff_2vbb = psd_eff_2vbb(config)
    
        if eff_2vbb:
            logger.info("Successfully calculated PSD eff. for 2vbb events.")
        else:
            logger.warning("PSD eff. for 2vbb failed.")

        # 2. PSD eff @ TH228 DEP
        eff_Th228 = psd_eff_dep(config, peak='Th228')
    
        if eff_Th228:
            logger.info("Successfully calculated PSD eff. for Th228 events.")
        else:
            logger.warning("PSD eff. for Th228 failed.")

        # 3. PSD eff @ TH228 DEP
        eff_Co56 = psd_eff_dep(config, peak='Co56')
        if eff_Co56:
            logger.info("Successfully calculated PSD eff. for Co56 events.")
        else:
            logger.warning("PSD eff. for Co56 failed.")

        # 4. PSD eff for timevariation
        eff_timevar = psd_eff_timevar(config)
        if eff_timevar:
            logger.info("Successfully calculated PSD eff. uncertainty due to variation over time.")
        else:
            logger.warning("PSD eff. for time stability failed.")


    # Calculate total efficiency
    efficiency, _, __ = psd_eff_qbb(
        config=config, 
        efficiencies_2vbb=eff_2vbb,
        efficiencies_Th228=eff_Th228, 
        efficiencies_Co56=eff_Co56, 
        efficiencies_Th228_timevar=eff_timevar
    )
    if efficiency:
        logger.info("Successfully calculated PSD eff. at Qbb.")
    else:
        logger.warning("PSD eff. at Qbb failed.")



    # Create plots
    if not config['save_figures']:
        logger.info("Done.")
        exit()
    
    logger.info("Create plots...")

    # Plot overview of total efficiency
    visualization.plot_total_eff(
        config=config, 
        results=efficiency
        )

    # Plot per-detector variance - Th-228
    visualization.plot_per_detector_run(
        config=config,
        eff_summary=eff_Th228['summary'],
        x_labels=sorted(eff_Th228.keys())[:-1],
        delta_name="det",
        y_lims=(0.75, 1.02),
        save_path="PSD_eff/PSD_eff_DEP_Th.png" if config['save_figures'] else None
        )

    # Plot per-detector variance - Co-56
    visualization.plot_per_detector_run(
        config=config,
        eff_summary=eff_Co56['summary'],
        x_labels=sorted(eff_Co56.keys())[:-1],
        delta_name="det",
        y_lims=(0.3, 1.15),
        save_path="PSD_eff/PSD_eff_DEP_Co.png" if config['save_figures'] else None
    )

    # Plot per-run variance - Th-228

    visualization.plot_per_detector_run(
        config=config,
        eff_summary=eff_timevar['summary'],
        x_labels=sorted(eff_timevar.keys())[:-1],
        delta_name="run",
        y_lims=(0.8, 1),
        save_path="PSD_eff/PSD_eff_timevar.png" if config['save_figures'] else None
    )

    # Plot a peak fit - Th228
    visualization.plot_peak_with_residuals(
        config=config,
        efficiencies=eff_Th228['V00048A'],
        bin_range=(1570, 1615),
        plot_range=15, 
        res_range=3,
        bin_keV=1,
        peak_energy=1592.5,
        )
    
    # Plot a peak fit - Co56
    visualization.plot_peak_with_residuals(
        config=config,
        efficiencies=eff_Co56['V00048A'],
        bin_range=(2215, 2247),
        plot_range=15, 
        res_range=3,
        bin_keV=1,
        peak_energy=2231.0,
        )
    
    logger.info("Done.")
    


#######################
# Call main function: #
#######################

if __name__ == '__main__':
    main()