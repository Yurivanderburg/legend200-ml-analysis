"""
Test Transformer model performance on test data.

This script create confusion matrices and estimates the data imbalance
for a transformer model trained in this work.
"""

# --- Standard library ---
import numpy as np
import os
import sys

# --- Third party ---
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
import awkward as ak
from sklearn.metrics import confusion_matrix
import seaborn as sns

# --- Project modules ---
from utils import set_style
from utils.data_io import load_pickle_file


# Parameters
DATA_PATH = "path/to/labelled_dataset.pkl"


# ---------------------------------------------------------------------
# Figure-specific plotting function
# ---------------------------------------------------------------------
def plot_calibration_spectrum(ax=None, save_path=None):
    """
    Generate the confusion and data imbalance plots.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Axis to plot on. If None, a new figure and axis are created.
    save_path : str, optional
        File path to save the figure. If None, the figure is not saved.

    Returns
    -------
    ax : matplotlib.axes.Axes
        Axis containing the plot.
    """

    data = load_pickle_file(DATA_PATH)
    

    # Convert to pandas dataframe
    df = ak.to_dataframe(data[['Prediction', 'label_true', 'peak']])
    df = df.replace({
        'DEP_Tl': 'DEP (Tl)',
        'SEP_Tl': 'SEP (Tl)',
        'FEP_Tl': 'FEP (Tl)',
        'FEP_Bi': 'FEP (Bi)'
    })

    # Data imbalance plot
    order = ['DEP (Tl)', 'FEP (Bi)', 'SEP (Tl)', 'FEP (Tl)']
    ct = pd.crosstab(df.peak, df.label_true, normalize=True)
    ct = ct.reindex(index=order)

    if ax is None:
        fig, ax = ct.plot(kind='bar', figsize=(10, 5))
    
    ax.set_xlabel("Calibration Peak")
    ax.set_ylabel("Normalized counts")

    ax.legend(['MSE (39%)', 'SSE (27%)', "n-contact (31%)", 'p-contact (3%)'], title="Label", fontsize=13)

    if save_path:
        os.makedirs(save_path, exist_ok=True)
        plt.savefig(os.path.join(save_path, "Data_balance_labels.png"), bbox_inches="tight")


    # Confusion matrix
    cm = confusion_matrix(df.label_true, df.Prediction, normalize='pred')
    labels=np.array(['MSE', 'SSE', 'ncontact', 'pcontact'])
    plt.figure(figsize=(10,8))
    sns.heatmap(cm, annot=True,cmap='Blues', xticklabels=labels, yticklabels=labels, norm=LogNorm())#  fmt='d', 
    plt.xlabel('Predicted Label')
    plt.ylabel('True Label')
    plt.title(f"Transformer 1")
    plt.tight_layout()

    if save_path:
        plt.savefig(os.path.join(save_path, "Results_confusionmatrix.png"), bbox_inches="tight")


    return ax

# ---------------------------------------------------------------------
# Main wrapper
# ---------------------------------------------------------------------

def main():
    """Run the script to generate and save the figure."""

    if not os.path.exists(DATA_PATH):
        sys.exit("Dataset not available. Please provide 'calibration_data.lh5' to run this script.")

    set_style()
    plot_calibration_spectrum(save_path="thesis_figures/output/05_psd/")


if __name__ == "__main__":
    main()