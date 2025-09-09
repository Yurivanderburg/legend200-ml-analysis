"""
Data I/O utilities for LEGEND-200 ML Analysis.

Includes:
- Pickle helpers
- HDF5 helpers
- LH5 helpers
"""
import awkward as ak
import pickle
from tqdm import tqdm
import numpy as np
import h5py
from lgdo import lh5 



# ------------------------------------------------------------
# Pickle I/O
# ------------------------------------------------------------
def load_pickle_file(files: str | list) -> ak.Array:
    """
    Load data from one or multiple pickle files.

    Parameters
    ----------
    files : str or list of str
        Path to a `.pkl` file or a list of file paths.

    Returns
    -------
    ak.Array
        Loaded data. If a list is given, files are concatenated.

    Raises
    ------
    TypeError
        If input is not str or list.
    """
    if isinstance(files, str):
        with open(files, 'rb') as file_loader:
            data_out = pickle.load(file_loader)
    elif isinstance(files, list):
        files_to_append = []
        for file in files:
            with open(file, 'rb') as file_loader:
                data_ = pickle.load(file_loader)
                files_to_append.append(data_)
        data_out = ak.concatenate(files_to_append)
    else:
        raise TypeError("File must be either str or list.")

    return data_out



def save_pickle_file(obj, path: str) -> None:
    """
    Save an object to a pickle file.

    Parameters
    ----------
    obj : any
        Object to be pickled.
    path : str
        Output file path.
    """
    with open(path, "wb") as f:
        pickle.dump(obj, f)



# ------------------------------------------------------------
# HDF5 I/O
# ------------------------------------------------------------
def load_hdf5_file(filenames: list, load_waveforms: bool=False) -> ak.Array:
    """
    Load datasets from an HDF5 file into an awkward array.

    Parameters
    ----------
    path : str
        Path to the `.hdf5` file.

    Returns
    -------
    data
        Awkward Array
    """
    data_dict = {}
    data_list = []

    if not isinstance(filenames, list):
        raise TypeError("Filenames must be a list.")

    for file in tqdm(filenames):

        with h5py.File(file, "r") as hdf_file:
            for field in hdf_file.keys():
                if field == 'waveform':
                    if load_waveforms:
                        data_dict[field] = np.array(hdf_file[field])
                else:
                    try:
                        data_dict[field] = np.array(hdf_file[field])
                        if isinstance(data_dict[field][0], bytes):
                            data_dict[field][:] = [x.decode('utf-8') for x in data_dict[field]]
                        data_dict[field] = data_dict[field].tolist()
                    

                    except:
                        raise ValueError("Field {} not present in file.".format(field))


        data_list.append(ak.Array(data_dict))
    
    # Concatanate all arrays
    data = ak.concatenate(data_list, axis=0)

    return data



# ------------------------------------------------------------
# LH5 I/O
# ------------------------------------------------------------
def load_lh5_file(data: str | list, channel_id: str, filenames):
    """
    Load data from one or multiple LH5 files for a given channel.

    Parameters
    ----------
    data : str or list of str
        Path to a `.lh5` file or a list of file paths.
    channel_id : str
        Channel identifier (e.g. "ch001") used to select the data group.

    Returns
    -------
    ak.Array
        Loaded data. If multiple files are given, arrays are concatenated.

    Raises
    ------
    NotImplementedError
        If `data` is not a str or list.
    """
    if isinstance(data, str):
        data_out = lh5.core.read_as(
            name = f"{channel_id}/hit",
            lh5_file = data,
            library = "ak" 
        )
    elif isinstance(data, list):
        data_list = []
        for ele in data:
            data_list.append(lh5.core.read_as(
                name = f"{channel_id}/hit",
                lh5_file = ele,
                library = "ak" 
                ))
        data_out = ak.concatenate(data_list)
    else:
        raise NotImplementedError("Argument data must be list or str.")
    
    return data_out



def load_remage_lh5_files(files: str | list) -> ak.Array:
    """ Load .lh5 files into an awkward array specifically for ReMaGe"""

    if isinstance(files, str):
        for det_id in lh5.tools.ls(files, lh5_group='/stp/')[:-1]:

            data = lh5.core.read_as(
                        name = det_id,
                        lh5_file = files,
                        library = "ak"
            )
            data['detid'] = det_id[3:] # remove the 'det'

    elif isinstance(files, list):
        data_to_append = []
        for file in tqdm(files):
            for detector_id in lh5.tools.ls(file, lh5_group='/stp/')[:-1]: # Last one is called 'vertices', so not a detector
                                
                data_ = lh5.core.read_as(
                            name = detector_id,
                            lh5_file = file,
                            library = "ak" 
                        )
                data_['detid'] = detector_id[3:] # remove the 'det'
                
                data_to_append.append(data_)
        data = ak.concatenate(data_to_append)

    else:
        raise TypeError("Parameters files must be either string or list!")

    return data
