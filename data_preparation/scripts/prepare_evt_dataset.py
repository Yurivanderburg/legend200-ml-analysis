"""
Prepare LEGEND-200 calibration data (EVT-tier) for ML analysis.

Workflow:
1. Load EVT files and apply basic event-level filters.
2. Select valid detectors based on metadata and config.
3. Match EVT events to RAW waveforms via timestamps.
4. Merge waveform data into EVT records.
5. Save results as Awkward Arrays serialized in pickle files.
"""

# --- Standard library ---
import os
import sys
import glob
import json
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
    Prepare LEGEND-200 EVT-tier data for ML analysis.

    - Loads EVT-tier data (high-level events).
    - Filters events (no pulsers, single HPGe multiplicity, QC passing).
    - Matches EVT timestamps with RAW waveforms.
    - Produces Awkward Arrays with merged metadata and waveforms.
    """

    # ---- CLI arguments ----
    parser = argparse.ArgumentParser(
        prog="prepare_dataset",
        description="Prepare LEGEND-200 calibration data (EVT-tier)",
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
        filename=f"DataPreparationEVT_{now}.log",
        archive_logs=True,
        log_folder="logs/prepare/evt",
    )
    logger.info("Starting EVT-tier data preparation...")

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

    # ---- Restriction ----
    if config["energy_mode"] != "continuum":
        raise NotImplementedError(
            "EVT-tier processing currently supports only `continuum` mode. "
            "Please use a manually defined continuous energy region."
        )

    # ---- Loop over periods ----
    for period in config["DIRS"]["PERIODS"]:
        run_dirs = sorted(
            os.listdir(os.path.join(config["DIRS"]["LEGEND_SKM_DIR"], period))
        )

        for run in run_dirs:
            period_run = f"{period}_{run}"
            filename_out = f"Data_{period_run}_cal.pkl"
            logger.info(f"Processing {period_run}...")

            # ---- Locate SKM/EVT and RAW files ----
            skm_files = sorted(
                glob.glob(
                    os.path.join(config["DIRS"]["LEGEND_SKM_DIR"], period, run, "*.lh5")
                )
            )
            raw_files = sorted(
                glob.glob(
                    os.path.join(config["DIRS"]["LEGEND_RAW_DIR"], period, run, "*.lh5")
                )
            )
            if not skm_files or not raw_files:
                logger.warning(f"Missing files for {period_run}, skipping.")
                continue

            time_of_datataking = raw_files[0].split("/")[-1].split("-")[4]
            chmap = lmeta.channelmap(time_of_datataking)

            # ---- Load EVT data ----
            evt_data, _ = store.read(name="evt", lh5_file=skm_files[0])
            evt_dataset = evt_data.view_as("ak")

            # ---- Event-level filtering ----
            evt_dataset = evt_dataset[
                ~evt_dataset.coincident.puls  # exclude pulsers
                & (evt_dataset.geds.multiplicity == 1)  # single HPGe signal
                & evt_dataset.geds.quality.is_bb_like  # QC filter
            ]

            # ---- Energy window selection ----
            low, high = config["continuous_region"]
            evt_dataset = evt_dataset[
                (evt_dataset.geds.energy[:, 0] > low)
                & (evt_dataset.geds.energy[:, 0] < high)
            ]

            if len(evt_dataset) == 0:
                logger.warning(f"No events remain after filtering for {period_run}.")
                continue

            # ---- Group timestamps by detector ----
            evt_timestamps = {}
            for entry in evt_dataset:
                det_id = entry.geds.rawid[0]
                evt_timestamps.setdefault(det_id, []).append(entry.trigger.timestamp)

            # ---- Valid detectors (ICPC) ----
            icpc_detectors = [
                v.daq.rawid for k, v in chmap.items() if k[0] == "V"
            ]

            data_list = []
            for det_id, timestamps in evt_timestamps.items():
                if int(det_id) not in icpc_detectors:
                    continue

                logger.info(f"Matching RAW data for detector {det_id}.")

                # Load RAW timestamps
                raw_timestamps, _ = store.read(
                    name=f"/ch{det_id}/raw",
                    lh5_file=raw_files,
                    field_mask=["timestamp"],
                )
                raw_timestamps = raw_timestamps.view_as("ak")

                # Match EVT timestamps with RAW
                mask = np.isin(raw_timestamps["timestamp"], timestamps)
                indices = np.where(mask)[0]
                logger.info(
                    f"Found {len(indices)} matches for detector {det_id} in {period_run}."
                )

                if len(indices) == 0:
                    continue

                # Load waveforms for matching events
                raw_data, _ = store.read(
                    name=f"/ch{det_id}/raw",
                    lh5_file=raw_files,
                    idx=indices,
                    field_mask=["waveform_windowed", "baseline", "timestamp"],
                )
                raw_data = raw_data.view_as("ak")

                # Temporary dataset for merging
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

            if not data_list:
                logger.warning(f"No matching detectors for {period_run}.")
                continue

            # ---- Merge EVT + RAW ----
            data_out = ak.concatenate(data_list)
            lookup = {
                entry["timestamp"]: {
                    "waveform": entry["waveform"],
                    "det_id": entry["det_id"],
                    "baseline": entry["baseline"],
                    "timestamp_r": entry["timestamp_r"],
                }
                for entry in ak.to_list(data_out)
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
                    for entry in ak.to_list(evt_dataset)
                ]
            )

            # ---- Final validation and save ----
            final_dataset = merged_dataset[
                ~ak.is_none(merged_dataset.timestamp_r)
            ]
            if len(final_dataset) == 0:
                logger.warning(f"No merged events for {period_run}.")
                continue

            out_path = os.path.join(config["DIRS"]["output_dir"], filename_out)
            with open(out_path, "wb") as f:
                pickle.dump(final_dataset, f)
            logger.info(f"Saved {len(final_dataset)} events to {filename_out}.")

    logger.info("Finished EVT-tier data preparation.")
    log_utils.close_log_handlers(logger)


if __name__ == "__main__":
    main()
