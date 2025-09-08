"""
Convert ReMaGe `.lh5` simulation files into a PET-format CSV.

"""

# --- Standard library ---
import os
import glob
import argparse

# --- Third-party ---
import awkward as ak
import pandas as pd

# --- Project-specific ---
from utils.data_io import load_remage_lh5_files


# Parameters
input_files = sorted(glob.glob("path/to/lh5_files*.lh5"))
output_file = "path/to/output.csv"
det_pos = [109.0000, -188.7935, 439.6895] # x,y,z coordinates for corresponding detector



def main():
    """ Create PET file from lh5 files."""

    # Load data
    data = load_remage_lh5_files(input_files)

    # Create csv file
    header_text = "#class tools::wcsv::ntuple\n#title Energy and time\n#separator 44\n#vector_separator 59\n#column double Event#\n#column double Detector_ID\n#column double Time\n#column double EnergyDeposit\n#column double X\n#column double Y\n#column double Z\n"


    print("Writing csv file")
    with open(output_file,'w') as f:
        f.write(header_text)


    # Transform array to pandas and bring in correct format
    print("Transform awkward array to Dataframe")
    df = ak.to_dataframe(data)
    print("Format dataframe")
    df.rename(columns=
            {
                'evtid': 'Event#',
                'detid': 'Detector_ID',
                'edep': 'EnergyDeposit',
                'time': 'Time',
                'xloc': 'X',
                'yloc': 'Y',
                'zloc': 'Z',
            }, inplace=True)
    df['Detector_ID'] = df['Detector_ID'].str[4:]

    # Transform coordinates
    for col, pos in zip(['X','Y','Z'], det_pos):
        df[col] = (df[col] * 1000) - pos

    # Perform coordinate transformation - because LegendGeSim differes from ReMaGe
    df.sort_values(by='Event#', inplace=True)
    df.drop(columns='particle', inplace=True)
    df.reset_index(inplace=True)


    # Append to previously created csv file
    print("Store PET file.")

    df.to_csv(output_file, 
            mode='a',
            index=False, 
            header=False,
            columns=['Event#', 'Detector_ID', 'Time', 'EnergyDeposit', 'X', 'Y', 'Z']
            )
    

if __name__ == "__main__":
    main()
