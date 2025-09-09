"""
Filter and label LEGEND-200 pre-processed data. Note that this script 
directly loads the output of prepare_dataset.py.

This script:
- Loads merged datasets (RAW, PHT, PSP, PAR) produced by prepare_dataset
- Attaches physics labels (MSE, SSE, ncontact, pcontact)
- Cleans tha data (peak finding and outlier removal)
- Applies cuts and saves filtered data to pickle files
"""

# --- Standard library ---
import os
import glob
import pickle
import traceback
import datetime as dt
import argparse
import numpy as np
import awkward as ak
import h5py
import sys
import pandas as pd

# --- LEGEND libraries ---
from legendmeta import LegendMetadata

# --- Project modules ---
from utils import log_utils
from utils.data_io import load_pickle_file
from data_preparation.data_prep_utils import config_utils, filter_utils

def main():
    """
    Attach labels and apply filtering to LEGEND-200 pre-processed data.
    """
    # ---- CLI arguments ----
    parser = argparse.ArgumentParser(
        prog="prepare_dataset",
        description="Prepare LEGEND-200 data for ML analysis",
        epilog="Written by Yuri van der Burg",
    )
    parser.add_argument("--config", required=True, help="Path to config file (.py)")
    parser.add_argument("--detector-type", required=False, help="Specify a detector type")
    args = parser.parse_args()

    # ---- Load config ----
    config = config_utils.load_config(args.config)

    # ---- Logger ----
    now = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    logger = log_utils.setup_logger(
        filename=f"DataFilterPHT_{now}.log",
        archive_logs=True,
        log_folder="logs/filter/pht",
    )
    logger.info("Starting data filtering...")

    # ---- Detector type ----
    detector_type = args.detector_type or config["detector_type"]
    if detector_type not in config["valid_detector_types"]:
        logger.error(f"Invalid detector type: {detector_type}")
        sys.exit(1)
    logger.info(f"Using detector type: {detector_type}")

    # ---- Validate config ----
    if not config_utils.check_config(config, logger):
        logger.error("Invalid config file, exiting.")
        sys.exit(1)

 
    # Detector types to process
    if config["FILTER"]["detector_types_to_filter"] == "all":
        detector_types = config["valid_detector_types"]
    else:
        detector_types = [
            d for d in config["FILTER"]["detector_types_to_filter"]
            if d in config["valid_detector_types"]
        ]

    # Loop over detector types
    for detector_type in detector_types:
        logger.info(f"Processing detector type: {detector_type}")

        in_dir = os.path.join(config["DIRS"]["output_dir"])
        out_dir = os.path.join(in_dir, config["DIRS"]["filtered_subdir"])

        # Periods to filter
        if config["FILTER"]["periods_to_filter"] == "all":
            period_list = sorted(os.listdir(in_dir))
        else:
            period_list = config["FILTER"]["periods_to_filter"]

        # Only keep valid period names (start with "p")
        period_list = [p for p in period_list if p.startswith("p")]

        # Ensure output directory exists
        os.makedirs(out_dir, exist_ok=True)
        if os.path.exists(out_dir):
            if config["FILTER"]["replace_files"]:
                logger.warning(f"Output dir {out_dir} exists. Files may be overwritten.")

        # Track statistics per detector type
        statistics = {}

        # Loop over periods
        for period in period_list:
            n_events_tot = 0
            n_events_counter = 0
            ct = 0
            detectors = sorted(os.listdir(os.path.join(in_dir, period)))
            statistics[period] = {}
            data_list_final = []

            # Loop over detectors
            for detector in detectors:
                iteration = f"{period}_{detector}"
                logger.info(f"Processing {iteration}")

                # Load pickled files for this detector/period
                file_list = sorted(glob.glob(os.path.join(in_dir, period, detector, "*.pkl")))
                if not file_list:
                    logger.warning("No files found. Skipping detector.")
                    continue

                # Load and merge pickle files
                data = load_pickle_file(file_list)
                logger.info(f"Merged {len(data)} events.")

                filename_out = f"{out_dir}/Data_{period}.pkl"
                filename_meta = f"{out_dir}/peak_param/Findpeaks_param_{period}_{detector}.pkl"

                if not config["FILTER"]["replace_files"] and os.path.isfile(filename_out):
                    logger.warning(f"{filename_out} exists. Skipping.")
                    continue

                # Peak finding
                try:
                    data = filter_utils.peak_finding(data=data, logger=logger, config=config)
                    logger.info("Peak finding applied.")
                except Exception:
                    logger.error("Peak finding failed.")
                    logger.error(traceback.format_exc())
                    continue

                n_prefilter = len(data)
                n_peakval = len(data[data.peak_val == True])

                if config["debug_mode"]:
                    logger.info("Peak validation label attached.")
                else:
                    data = data[data.peak_val == True]
                    logger.info("Applied peak validation cut.")

                if n_peakval == 0:
                    logger.error("No events after peak validation.")
                    continue

                # Outlier removal
                try:
                    data = filter_utils.remove_outliers(
                        data=data, 
                        logger=logger, 
                        filename_meta=filename_meta, 
                        config=config
                    )
                    logger.info("Outliers removed.")
                except Exception:
                    logger.error("Outlier removal failed.")
                    logger.error(traceback.format_exc())
                    continue

                n_filtered = len(data)
                statistics[period][detector] = {
                    "events_pre_filter": n_prefilter,
                    "events_post_peakval": n_peakval,
                    "events_filtered": n_filtered,
                }

                if n_filtered == 0:
                    logger.error("No events after outlier removal.")
                    continue

                # Sort by timestamp
                data = data[ak.argsort(data.timestamp)]

                # Merge surface labels if enabled
                if config["FILTER"]["merge_surface_labels"]:
                    data["label"] = ak.str.replace_substring(data["label"], "pcontact", "surface")
                    data["label"] = ak.str.replace_substring(data["label"], "ncontact", "surface")

                data_list_final.append(data)
                n_events_counter += n_filtered
                n_events_tot += n_filtered

                # Prevent OOM by writing partial results
                if n_events_counter > 5e5:
                    try:
                        data_notfinal = ak.concatenate(data_list_final, mergebool=True)
                        with open(f"{out_dir}/Data_{period}_{ct}.pkl", "wb") as fout:
                            pickle.dump(data_notfinal, fout)
                        logger.info("Stored partial results to prevent OOM.")
                        ct += 1
                        data_list_final = []
                        n_events_counter = 0
                    except Exception:
                        logger.error("Failed to store partial results.")
                        logger.error(traceback.format_exc())

            # Merge and save final dataset
            if data_list_final:
                data_final = ak.concatenate(data_list_final, mergebool=True)
                if config["SAVE"]["as_pickle"]:
                    with open(filename_out, "wb") as fout:
                        pickle.dump(data_final, fout)
                    logger.info(f"Stored final dataset for period {period}.")
            else:
                logger.warning(f"No data to save for period {period}. Skipping.")
                continue

            statistics[period]["events_tot"] = n_events_tot
            with open(f"{in_dir}/N_events_{detector_type}.pkl", "wb") as fout:
                pickle.dump(statistics, fout)

            logger.info(
                f"Finished {period}: {np.round(n_events_tot/1e6, 3)}M events total."
            )

    logger.info("Filtering script finished.")


if __name__ == "__main__":
    main()
