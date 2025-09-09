"""
Convert LEGEND-200 pickle datasets into an HDF5 file.

Usage
-----
python -m data_preparation.scripts.pickle_to_hdf5 input.pkl output.h5

- Loads one or more Awkward Arrays from pickle files
- Converts fields into HDF5 datasets
- Handles numeric and string fields, plus fixed-length waveforms
"""

# --- Standard library ---
import os
import glob
import pickle
import datetime as dt
import h5py
import awkward as ak
import numpy as np
from tqdm import tqdm
import argparse

# --- Project modules ---
from utils.data_io import load_pickle_file
from utils import log_utils


# --- Output directory ---
OUT_DIR = "output/dir/"
os.makedirs(OUT_DIR, exist_ok=True)



def pickle_to_hdf5(input_files: list[str], output_file: str, logger):
    """
    Merge pickle files into a single HDF5 file.
    """

    with h5py.File(output_file, "a") as h5f:
        for file in tqdm(input_files, desc="Processing pickle files"):
            data = load_pickle_file(file)

            # Drop fields not present everywhere
            data = data[[x for x in ak.fields(data) if x != "tail_pars"]]

            # Fix formatting of 'p' field if possible
            if "p" in ak.fields(data):
                try:
                    data["p"] = np.array(data["p"][:, 1:]).astype(int)
                except Exception:
                    pass

            # Process all fields
            for field in ak.fields(data):
                try:
                    field_data = data[field]
                    field_data_numpy = ak.to_numpy(field_data)

                    # Skip empty fields
                    if field_data_numpy.size == 0:
                        logger.info(f"Skipping empty field '{field}' in {file}")
                        continue

                    # Handle string fields
                    if np.issubdtype(field_data_numpy.dtype, np.str_):
                        dtype = h5py.string_dtype(encoding="utf-8")
                        if field not in h5f:
                            h5f.create_dataset(
                                field,
                                shape=(0,),
                                maxshape=(None,),
                                dtype=dtype,
                                chunks=True,
                                compression="gzip",
                            )

                        dset = h5f[field]
                        chunk_size = len(field_data_numpy)
                        dset.resize(dset.shape[0] + chunk_size, axis=0)
                        dset[-chunk_size:] = field_data_numpy

                    # Handle numeric fields
                    else:
                        if field not in h5f:
                            if field == "waveform":
                                # Fixed-length waveforms
                                h5f.create_dataset(
                                    field,
                                    shape=(0, 1400),
                                    maxshape=(None, 1400),
                                    dtype=field_data_numpy.dtype,
                                    chunks=True,
                                    compression="gzip",
                                )
                            else:
                                # Generic numeric field
                                h5f.create_dataset(
                                    field,
                                    shape=(0,),
                                    maxshape=(None,),
                                    dtype=field_data_numpy.dtype,
                                    chunks=True,
                                    compression="gzip",
                                )

                        dset = h5f[field]
                        chunk_size = len(field_data_numpy)

                        # Validate waveform shape
                        if field == "waveform" and field_data_numpy.shape[1] != 1400:
                            logger.info(f"Skipping waveform with invalid shape in {file}")
                            continue

                        # Extend dataset
                        dset.resize(dset.shape[0] + chunk_size, axis=0)
                        dset[-chunk_size:] = field_data_numpy

                except Exception as e:
                    logger.info(f"Error processing field '{field}' in {file}: {e}")


def main():
    # ---- CLI arguments ----
    parser = argparse.ArgumentParser(
        description="Convert LEGEND-200 pickle files into HDF5."
    )
    parser.add_argument("input", nargs="+", help="Input pickle file(s)")
    parser.add_argument("output", help="Output HDF5 filename")
    args = parser.parse_args()

    # ---- Logger ----
    now = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    logger = log_utils.setup_logger(
        filename=f"DataConversion_{now}.log",
        archive_logs=True,
        log_folder="logs/convert/",
    )
    logger.info("Start converting pickle to hdf5 ...")

    pickle_to_hdf5(args.input, args.output, logger)



if __name__ == "__main__":
    main()