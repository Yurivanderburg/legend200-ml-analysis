"""
Prepare LEGEND-200 physics data for ML analysis.

Workflow:
1. Load SKM files (summary-level data).
2. Extract per-detector event timestamps.
3. Match events to RAW waveforms using timestamps.
4. Merge SKM attributes with waveform data.
5. Store results as Awkward Arrays serialized in pickle files.
"""

# --- Standard library ---
import os
import sys
import glob
import pickle
import argparse
import logging
import datetime as dt

# --- Third-party ---
import numpy as np
import awkward as ak
from legendmeta import LegendMetadata
import lgdo.lh5 as lh5

# --- Project modules ---
from utils import log_utils
from data_preparation.data_prep_utils import config_utils


def main():
    """
    Prepare LEGEND-200 data for ML analysis.

    - Loads SKM summary data and RAW waveforms.
    - Matches events via timestamps for valid detectors.
    - Merges waveform and metadata into Awkward Arrays.
    - Stores results as pickle files.
    """

    # ---- CLI arguments ----
    parser = argparse.ArgumentParser(
        prog="prepare_dataset",
        description="Prepare LEGEND-200 data for ML analysis",
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
        filename=f"DataPreparationPHY_{now}.log",
        archive_logs=True,
        log_folder="logs/prepare/phy",
    )
    logger.info("Starting data preparation...")

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

    # ---- LH5 Store ----
    store = lh5.LH5Store()

    # ---- Loop over periods ----
    for period in config["DIRS"]["PERIODS"]:
        skm_files = sorted(
            glob.glob(os.path.join(config["DIRS"]["LEGEND_SKM_DIR"], f"*{period}*.lh5"))
        )

        for i, skm_file in enumerate(skm_files):
            run = skm_file.split("/")[-1].split("-")[2]
            period_run = f"{period}_{run}"
            filename_out = f"Data_{period_run}_phy.pkl"

            logger.info(f"Processing {period_run} ({i+1}/{len(skm_files)}).")

            # ---- Load RAW files ----
            raw_files = sorted(
                glob.glob(
                    os.path.join(config["DIRS"]["LEGEND_RAW_DIR"], period, run, "*.lh5")
                )
            )
            time_of_datataking = raw_files[0].split("/")[-1].split("-")[4]
            chmap = lmeta.channelmap(time_of_datataking)

            # ---- Load SKM data ----
            skm_data, _ = store.read(name="skm", lh5_file=skm_file)
            skm_dataset = skm_data.view_as("ak")

            # ---- Group SKM timestamps per detector ----
            skm_timestamps = {}
            for entry in ak.to_list(skm_dataset.geds.rawid):
                skm_timestamps[entry] = []

            for event in skm_dataset:
                skm_timestamps[event.geds.rawid].append(event.trigger.timestamp)

            # ---- Find valid detectors (ICPC, active) ----
            icpc_detectors = [
                v.daq.rawid for k, v in chmap.items() if k[0] == "V"
            ]

            # ---- Match SKM events to RAW waveforms ----
            data_list = []
            for det_id, timestamps in skm_timestamps.items():
                if int(det_id) not in icpc_detectors:
                    continue

                logger.info(f"Matching events for detector {det_id} in {period_run}.")

                # Load RAW timestamps
                raw_timestamps, _ = store.read(
                    name=f"/ch{det_id}/raw",
                    lh5_file=raw_files,
                    field_mask=["timestamp"],
                )
                raw_timestamps = raw_timestamps.view_as("ak")

                # Match RAW ↔ SKM by timestamp
                mask = np.isin(raw_timestamps["timestamp"], timestamps)
                indices = np.where(mask)[0]
                logger.info(f"Found {len(indices)} matches for detector {det_id}.")

                # Load RAW waveforms for matching indices
                raw_data, _ = store.read(
                    name=f"/ch{det_id}/raw",
                    lh5_file=raw_files,
                    idx=indices,
                    field_mask=["waveform_windowed", "baseline", "timestamp"],
                )
                raw_data = raw_data.view_as("ak")

                # Build temporary dataset
                data_temp = ak.Array(
                    {
                        "timestamp": raw_data.timestamp,
                        "waveform": raw_data.waveform_windowed.values,
                        "det_id": det_id,
                        "baseline": raw_data.baseline,
                        "timestamp_r": raw_data.timestamp,
                    }
                )
                data_list.append(data_temp)

            # ---- Merge SKM + RAW ----
            if data_list:
                data_raw = ak.concatenate(data_list)

                # Lookup table: timestamp → waveform data
                lookup = {
                    entry["timestamp"]: {
                        "waveform": entry["waveform"],
                        "det_id": entry["det_id"],
                        "baseline": entry["baseline"],
                        "timestamp_r": entry["timestamp_r"],
                    }
                    for entry in ak.to_list(data_raw)
                }

                merged_dataset = ak.Array(
                    [
                        {
                            **entry,
                            **lookup.get(
                                entry["trigger"]["timestamp"],
                                {
                                    "waveform": None,
                                    "det_id": None,
                                    "baseline": None,
                                    "timestamp_r": None,
                                },
                            ),
                        }
                        for entry in ak.to_list(skm_dataset)
                    ]
                )

                # Keep only ICPC events
                final_dataset = merged_dataset[~ak.is_none(merged_dataset.timestamp_r)]

                if len(final_dataset) == 0:
                    logger.warning(f"No matched events in {period_run}.")
                else:
                    # Save to pickle
                    out_path = os.path.join(config["DIRS"]["output_dir"], filename_out)
                    with open(out_path, "wb") as f:
                        pickle.dump(final_dataset, f)
                    logger.info(
                        f"Saved {len(final_dataset)} events to {filename_out}."
                    )
            else:
                logger.warning(f"No matching detectors in {period_run}.")

    logger.info("Finished executing data preparation script.")
    log_utils.close_log_handlers(logger)


if __name__ == "__main__":
    main()
