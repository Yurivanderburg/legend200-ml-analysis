"""
Unsupervised learning - UMAP clustering of waveforms

Here, I tried some unsupervised learning (UMAP) to see whether
the different labels are somewhat detected in unsupervised learning.
"""

# --- Standard library ---
import numpy as np
import os
import sys

# --- Third party ---
import matplotlib.pyplot as plt
import awkward as ak
import pandas as pd
from scipy.optimize import curve_fit
import umap.umap_ as umap
from sklearn.cluster import KMeans
from hdbscan import HDBSCAN


# --- Project modules ---
from utils import set_style
from utils.data_io import load_pickle_file
from utils.math import linear


# Parameters
DATA_PATH = "path/to/labelled_dataset.pkl"
CHANNEL_ID = 'ch1116803'

# ---------------------------------------------------------------------
# Figure-specific plotting function
# ---------------------------------------------------------------------
def weighted_moving_average(arr, half_win_samps:int, gauss_sigma_samps:float):
    """
    This function produce a smoothed array of the input array.
    """
    ser = pd.Series(np.array(arr))

    output = ser.rolling(window=2*half_win_samps+1, 
                         center=True, 
                         closed='both', 
                         min_periods=half_win_samps+1, 
                         win_type='gaussian'
                         ).mean(std=gauss_sigma_samps)
    return output



def normalize(data: ak.Array, 
                  gradient=False, align_tp0=False, gaussian_smearing=False) -> ak.Array:
    """ Normalize waveforms and optionally calculate gradient (w or wo smearing) and align waveforms."""
    
    # First, waveform alignment
    wfs = data.waveform
    tp_0s = (data.geds.t0[:,0] - 42000)/16
    x_bsl = np.linspace(0, 249, 250)
    lis = []
    for wf, tp0 in zip(wfs, tp_0s):
        if align_tp0:
            wf_new = wf[(tp0-250):-(300 - (tp0-250))]
            if len(wf_new) != 1100:
                print("Error")
                wf_new = wf[150:-150]
        else:
            wf_new = wf

        # Secondly, smear the signal slightly (if necessary)
        if gaussian_smearing:
            wf_new = weighted_moving_average(wf_new, half_win_samps=2, gauss_sigma_samps=2)

        # Thirdly, calculate gradient (if necessary)
        if gradient:
            wf_new = np.gradient(wf_new)

        # Normalization must always happen
        y = wf_new[:250]
        popt, _ = curve_fit(linear, x_bsl, y)
        wf_new_subtr = wf_new - popt[1]   
        wf_norm = (wf_new_subtr - np.min(wf_new_subtr))/(np.max(wf_new_subtr) - np.min(wf_new_subtr))
        lis.append(wf_norm)

    data['wf_norm'] = np.array(lis)
    
    return data



def create_masks(data_out, total_data, type=None):
    """
    Create a mask of the embedding.
    """
    embeddings = {}
    
    if type == 'label':
        for label in ['SSE', 'MSE', 'surface']:
            mask = (total_data.label == label)

            embeddings[label]  = data_out[mask]

    elif type == 'peak':
        for peak in ['DEP_Tl', 'FEP_Tl', 'SEP_Tl', 'FEP_Bi']:
            mask = (total_data.peak == peak)

            embeddings[peak]  = data_out[mask]
    else:
        raise NotImplementedError()
        
    return embeddings   



def create_umap(data_in: ak.Array, 
                total_data: ak.Array, 
                type: str=None, 
                n_neighbors: int=15, 
                min_dist: float=0.1, 
                n_components: int=2, 
                metric: str='euclidean',
                plot_title: str=None,
                save_path: str=None
                ):
    """
    Create a UMAP instance for given data
    """
    
    print("Create UMAP instance")
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, 
                        n_components=n_components, metric=metric)
    reducer.fit(data_in)
    embedding = reducer.transform(data_in)
    assert(np.all(embedding == reducer.embedding_))
    
    print("Prepare plots")
    embeddings = create_masks(embedding, total_data, type=type)


    fig = plt.figure(figsize=(8,6))
    for peak_name, peak_data in embeddings.items():
        if n_components == 1:
            plt.scatter(peak_data[:, 0], range(len(peak_data)), s=1, label=peak_name)
        elif n_components == 2: 
            plt.scatter(peak_data[:, 0], peak_data[:, 1], s=1, label=peak_name)
        elif n_components == 3:
            ax = fig.add_subplot(111, projection='3d')
            plt.scatter(peak_data[:,0], peak_data[:,1], peak_data[:,2], label=peak_name)
    plt.gca().set_aspect('equal', 'datalim')
    plt.title(plot_title, fontsize=20)
    plt.legend(fontsize=12, markerscale=10)
    plt.tight_layout()
    if save_path:
        plt.savefig(os.path.join(save_path, "Umap.png"), bbox_inches="tight")


    return embedding



def create_kmeans(data, umap, save_path):
    """
    For a given UMAP clustering, try to find clusters with kmeans
    """
    kmeans_wf_norm = KMeans(n_clusters=5, random_state=42, n_init='auto', max_iter=int(1e5)).fit(umap)

    # Save a sketch
    fig, ax = plt.subplots(1,1)
    ax.scatter(umap[:,0], umap[:,1], s=0.5, c=kmeans_wf_norm.labels_)
    ax.scatter(kmeans_wf_norm.cluster_centers_[:,0], kmeans_wf_norm.cluster_centers_[:,1], color="r")
    ax.set_title("UMAP kmeans clustering (wfs, norm)")
    if save_path:
        plt.savefig(os.path.join(save_path, "Umap_kmeans.png"), bbox_inches="tight")

    kmeans_wf_norm = kmeans_wf_norm.predict(umap)

    # Plot of single clusters:
    for i in np.arange(0, max(kmeans_wf_norm)+1):
        fig, ax = plt.subplots(1,2, figsize=(12,5))

        data_temp = data[(kmeans_wf_norm == i)]
        print(f"label: {i}")

        ax[0].scatter(umap[~(kmeans_wf_norm == i), 0],
                umap[~(kmeans_wf_norm == i), 1],
                color=(0.5, 0.5, 0.5),
                s=0.1,
                alpha=0.5)
        ax[0].scatter(umap[(kmeans_wf_norm == i), 0],
                umap[(kmeans_wf_norm == i), 1],
                c='red',
                s=0.1)
        ax[1].hist(data_temp.label)

        ax[0].set_title(f"kmeans cluster {i}: UMAP")
        ax[1].set_title(f"kmeans cluster {i}: distribution")

        ax[1].set_ylabel("Counts")
        ax[1].set_yscale("log")
        ax[1].set_ylim(1,2*len(data_temp))

        if save_path:
                plt.savefig(os.path.join(save_path, f"Umap_kmeans_clusters_{i}.png"), bbox_inches="tight")

    return kmeans_wf_norm



def create_hdbscan(data, umap, save_path):
    """
    For a given UMAP clustering, try to find clusters with hdbscan
    """
    hdbscan_wf_norm = HDBSCAN(
        min_cluster_size=15,
        min_samples=40,
        max_cluster_size=len(data[data.label == 'MSE']),
    ).fit_predict(umap)

    fig, ax = plt.subplots(1,1)
    ax.scatter(umap[~(hdbscan_wf_norm >= 0), 0],
                umap[~(hdbscan_wf_norm >= 0), 1],
                color=(0.5, 0.5, 0.5),
                s=0.1,
                alpha=0.5)
    ax.scatter(umap[(hdbscan_wf_norm >= 0), 0],
                umap[(hdbscan_wf_norm >= 0), 1],
                c=hdbscan_wf_norm[(hdbscan_wf_norm >= 0)],
                s=0.1)

    for i in np.arange(0, max(hdbscan_wf_norm)+1):
        fig, ax = plt.subplots(1,2, figsize=(12,5))

        data_temp = data[(hdbscan_wf_norm == i)]
        ax[0].scatter(umap[~(hdbscan_wf_norm == i), 0],
        umap[~(hdbscan_wf_norm == i), 1],
            color=(0.5, 0.5, 0.5),
            s=0.1,
            alpha=0.5)
        ax[0].scatter(umap[(hdbscan_wf_norm == i), 0],
        umap[(hdbscan_wf_norm == i), 1],
            c=hdbscan_wf_norm[(hdbscan_wf_norm == i)],
            s=0.1)
        ax[1].hist(data_temp.label)
    
        ax[0].set_title(f"hdbscan cluster {i}: UMAP")
        ax[1].set_title(f"hdbscan cluster {i}: distribution")

        ax[1].set_ylabel("Counts")
        ax[1].set_yscale("log")
        ax[1].set_ylim(1,2*len(data_temp))

        if save_path != None:
                plt.savefig(os.path.join(save_path, f"Umap_kmeans_clusters_{i}.png"), bbox_inches="tight")



def plot_umaps(ax=None, save_path=None):
    """
    Generate the umaps.

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
    data = data[data.geds.rawid[:,0] == CHANNEL_ID]

    # Create labels
    label = []
    for ele in data:
        if ele.geds.psd.is_bb_like[0] == True:
            label.append('SSE')
        else:
            if ele.geds.psd.high_aoe.is_bulk[0] == True:
                label.append('MSE')
            else:
                label.append('surface')


    # ---- Normalized, aligned gradients ----
    data_norm_align = normalize(data, 
                                gradient=True, 
                                align_tp0=True, 
                                gaussian_smearing=True
                                )
    umap_norm_align = create_umap(data_in=data_norm_align['wf_norm'], 
                                  total_data=data, 
                                  type='label',
                                  save_path=save_path
                                  )
    
    # ---- Normalized, not-aligned waveforms ----

    data_norm = normalize(data, 
                          gradient=False, 
                          align_tp0=False, 
                          gaussian_smearing=True
                          )
    umap_wf_norm = create_umap(data_in=data_norm['wf_norm'], 
                                total_data=data, 
                                type='label',
                                save_path=save_path
                                )

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(14, 5))

    return ax



# ---------------------------------------------------------------------
# Main wrapper
# ---------------------------------------------------------------------
def main():
    """Run the script to generate and save the figure."""

    if not os.path.exists(DATA_PATH):
        sys.exit("Dataset not available. Please provide 'calibration_data.lh5' to run this script.")

    set_style()
    plot_umaps(save_path="thesis_figures/output/misc/")


if __name__ == "__main__":
    main()