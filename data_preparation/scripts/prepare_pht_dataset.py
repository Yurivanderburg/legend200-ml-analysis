"""
Prepare LEGEND-200 calibration data (PHT/HIT-tier) for ML analysis.

Workflow:
1. Load filenames for RAW, PHT, TCM, PSP, PAR tiers.
2. Select valid detectors based on metadata and config.
3. Load datasets with cuts, masks, and calibration info.
4. Optionally attach labels (MSE/SSE, etc.).
5. Save results as Awkward Arrays in pickle files.
"""

# --- Standard library ---
import os
import sys
import json
import pickle
import argparse
import logging
import datetime as dt
import traceback

# --- Third-party ---
import awkward as ak
from legendmeta import LegendMetadata
import lgdo.lh5 as lh5

# --- Project modules ---
from utils import log_utils
from data_preparation.data_prep_utils import config_utils, prep_utils


def main():
    """
    Prepare LEGEND-200 calibration (HIT-tier) data for ML analysis.

    - Loads metadata and per-detector datasets.
    - Applies energy mode selection, cuts, and masks.
    - Optionally attaches MSE/SSE labels.
    - Stores results as pickle files.
    """

    # ---- CLI args ----
    parser = argparse.ArgumentParser(
        prog="prepare_dataset",
        description="Prepare LEGEND-200 calibration data (HIT-tier)",
        epilog="Written by Yuri van der Burg",
    )
    parser.add_argument("--config", required=True, help="Path to config file (.py)")
    parser.add_argument("--detector_type", type=str, help="Override detector type")
    args = parser.parse_args()

    # ---- Load config ----
    config = config_utils.load_config(args.config)

    # ---- Logger ----
    now = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    logger = log_utils.setup_logger(
        filename=f"DataPreparationPHT_{now}.log",
        archive_logs=True,
        log_folder="logs/prepare/pht",
    )
    logger.info("Starting calibration data preparation...")

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

    # ---- LEGEND Metadata ----
    os.environ["SSH_AUTH_SOCK"] = config["SSH_AUTH"]
    lmeta = LegendMetadata(path=config["DIRS"]["lmeta"])

    # ---- Loop over periods and runs ----
    for period in config["DIRS"]["PERIODS"]:
        period_path = os.path.join(
            config["DIRS"]["LEGEND_RAW_DIR"], config["DIRS"]["TYPE"], period
        )
        runs = sorted(os.listdir(period_path))

        for run in runs:
            period_run = f"{period}_{run}"

            # Skip excluded runs
            if period_run in config["exclude_period_run"]:
                logger.info(f"Skipping excluded run {period_run}.")
                continue

            # ---- Load tier filenames ----
            tier_files, valid = prep_utils.load_filenames(
                period=period, run=run, logger=logger, config=config
            )
            if not valid:
                logger.error(f"Missing tier files for {period_run}, skipping.")
                continue

            # ---- Detector metadata ----
            date_from_raw = tier_files["RAW"][0].split("/")[-1].split("T")[0][-8:]
            timestamp = dt.datetime.fromisoformat(date_from_raw)

            valid_detectors = prep_utils.filter_valid_detectors(
                legend_metadata=lmeta,
                timestamp=timestamp,
                detector_type=detector_type,
            )
            if not valid_detectors:
                logger.error(f"No valid {detector_type} detectors for {period_run}.")
                continue
            logger.info(
                f"Found {len(valid_detectors)} valid {detector_type} detectors in {period_run}."
            )

            # ---- Load PAR file if needed ----
            par_data = {}
            if config["energy_mode"] == "peaks":
                with open(tier_files["PAR"][0], "rb") as f:
                    par_data = json.load(f)

            # ---- Iterate over detectors ----
            for detector_name, item in valid_detectors.items():
                det_id = item["rawid"]
                manufacturer = item["manufacturer"]
                logger.info(
                    f"--- Run {period_run}, detector {detector_name} (ch{det_id}) ---"
                )

                # Output directory
                out_dir = os.path.join(
                    config["DIRS"]["output_dir"], detector_type, period, detector_name
                )
                os.makedirs(out_dir, exist_ok=True)
                if run == "r000":
                    logger.info(f"Output dir: {out_dir}")

                # ---- Manufacturer logic for labels ----
                ortec_type = (
                    manufacturer == "Ortec"
                    if detector_type == "ICPC"
                    else detector_type == "PPC"
                )

                # ---- Validate PAR peaks ----
                if config["energy_mode"] == "peaks":
                    if f"ch{det_id}" not in par_data:
                        logger.error(
                            f"Detector {detector_name} (ch{det_id}) not in PAR file."
                        )
                        continue
                    data_peaks = [
                        k.split(".")[0]
                        for k in par_data[f"ch{det_id}"]["results"]["ecal"][
                            "cuspEmax_ctc_runcal"
                        ]["pk_fits"]
                    ]
                    for config_peak in config["peaks"].values():
                        if str(config_peak) not in data_peaks:
                            logger.warning(
                                f"Peak {config_peak} not in PAR data for ch{det_id}."
                            )

                # ---- Load datasets ----
                # Apply cuts, masks, and merge tiers into ak.Array
                try:
                    dataset = prep_utils.load_data(
                        pht_files=tier_files["PHT"][:200],
                        raw_files=tier_files["RAW"][:200],
                        tcm_files=tier_files["TCM"],
                        psp_files=[]
                        if config["DIRS"]["TYPE"] == "phy"
                        else tier_files["PSP"],
                        par_data=par_data,
                        channel_id=det_id,
                        logger=logger,
                        config=config
                    )
                except Exception:
                    logger.error(f"Data loading failed for ch{det_id}.")
                    logger.error(traceback.format_exc())
                    continue

                dataset["p"] = period
                dataset["r"] = run

                # ---- Attach labels ----
                if config["add_labels"]:
                    try:
                        dataset = prep_utils.attach_labels(dataset, ortec_type=ortec_type)
                        logger.info("Labels attached successfully.")
                    except Exception:
                        logger.error("Failed to attach labels.")
                        logger.error(traceback.format_exc())
                        continue
                else:
                    logger.info("Skipping label attachment.")

                # ---- Save dataset ----
                if not isinstance(dataset, ak.Array):
                    logger.error("Dataset is not an Awkward Array.")
                    continue
                if len(dataset) == 0:
                    logger.warning("Dataset is empty.")
                    continue

                logger.info(f"Saving {len(dataset)} events to pickle...")
                filename = f"Data_{period_run}_ch{det_id}_pre"

                if len(dataset) > 1.5e6:
                    # Split into two files (temporary workaround)
                    with open(f"{out_dir}/{filename}_0.pkl", "wb") as f:
                        pickle.dump(dataset[: int(1.5e6)], f)
                    with open(f"{out_dir}/{filename}_1.pkl", "wb") as f:
                        pickle.dump(dataset[int(1.5e6) :], f)
                else:
                    with open(f"{out_dir}/{filename}.pkl", "wb") as f:
                        pickle.dump(dataset, f)

    logger.info("Finished calibration data preparation.")
    log_utils.close_log_handlers(logger)


if __name__ == "__main__":
    main()
