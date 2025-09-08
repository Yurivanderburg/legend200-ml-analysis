"""
Calculate summary based on ZeroNuFit.jl for multiple different fits
and create plots

Filpath structure:
    sensitivity/mcmc_files/*.json
    sensitivity/fake_data/*.json

Returns:
    summary.json
"""

# --- Standard library ---
import datetime as dt
import argparse
import json
import tqdm
import os

# --- Project modules ---
from utils import log_utils
from physics_analysis.sensitivity_study import summarize_results, visualize_sensitivity


def main():
    """
    Calculate summary and create plots
    """

    # ---- Parse path ----
    parser = argparse.ArgumentParser(description="Summarize LEGEND sensitivity results from JSONs.")
    parser.add_argument("--path", 
                        required=True, 
                        type=str, 
                        help="Path to directory containing JSON results"
                        )
    parser.add_argument("--run_calc", 
                        type=str,
                        default="no", 
                        help="Calculate summary file (yes/no)."
                        )

    args = parser.parse_args()

    # ---- Logger ----
    now = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    logger = log_utils.setup_logger(
        filename=f"Sensitivity_{now}.log",
        archive_logs=True,
        log_folder="logs/sens/",
    )

    
    if args.run_calc == 'yes':
        # Calculate summary for every file in given path
        logger.info("Calculate summary..")
        fits = os.listdir(args.path)
        all_results = {}
        for fit in tqdm.tqdm(fits):
            if os.path.isdir(os.path.join(args.path, fit)):
                all_results[str(fit)] = summarize_results(path=os.path.join(args.path, fit))
            else:
                logger.info(f"Skipping {fit} - not a directory")
        # Save results
        with open(os.path.join(args.path, "summary.json"), "w") as f:
            json.dump(all_results, f, indent=2)

    else:
        logger.info("Loading summary from file..")
        with open(os.path.join(args.path, "summary.json"), "r") as f:
            all_results = json.load(f)



    logger.info("Create plots...")
    os.makedirs("Plots/SensStudy/", exist_ok=True)
    visualize_sensitivity.plot_background_index(all_results)

    visualize_sensitivity.plot_signal_half_rate(all_results)

    for name, val_dict in all_results.items():
        if 'AoE' in name:
            visualize_sensitivity.plot_roi_counts(val_dict, method='AoE')
        else:
            visualize_sensitivity.plot_roi_counts(val_dict, method='Transformer')

    visualize_sensitivity.plot_sensitivity(all_results)

    logger.info("Done.")

    

#######################
# Call main function: #
#######################

if __name__ == '__main__':
    main()