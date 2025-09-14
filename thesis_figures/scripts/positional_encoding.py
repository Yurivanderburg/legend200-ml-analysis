"""
Figure: Positional encoding used in Transformer

This script creates a visualization of exact type of positional
encoding of the embeddings, as used in the Transformer of this
work.
"""

# --- Standard library ---
import numpy as np
import os

# --- Third party ---
import matplotlib.pyplot as plt
import torch
from torch import nn


# --- Project modules ---
from utils import set_style



# ---------------------------------------------------------------------
# Figure-specific plotting function
# ---------------------------------------------------------------------
class PositionalEncoding(nn.Module):
    def __init__(self, d_model: int, dropout: float = 0.1, max_len: int = 3800):
        """
        Initializes the PositionalEncoding module.

        Parameters:
        d_model (int): The dimension of the model (embedding size).
        dropout (float): The dropout rate. Default is 0.1.
        max_len (int): The maximum length of the input sequences. Default is 3800.
        """
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)

        # Compute the positional encodings once in log space.
        position = torch.arange(max_len).unsqueeze(1)  # Shape: [max_len, 1]
        div_term = torch.exp(torch.arange(0, d_model, 2) * (-np.log(10000.0) / d_model))

        pe = torch.zeros(1, max_len, d_model)  # Shape: [1, max_len, d_model]
        pe[0, :, 0::2] = torch.sin(position * div_term)  # Apply sin to even indices in the array
        pe[0, :, 1::2] = torch.cos(position * div_term)  # Apply cos to odd indices in the array

        self.register_buffer('pe', pe)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Forward pass for adding positional encoding to the input tensor.

        Parameters:
        x (torch.Tensor): The input tensor with shape [seq_len, batch_size, embedding_dim].

        Returns:
        torch.Tensor: The tensor with added positional encodings and dropout applied.
        """
        x = x + self.pe[:, :x.size(1)]
        return self.dropout(x)
        



def plot_positional_encoding(ax=None, save_path=None):
    """
    Generate the positional encoding figure.

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

    # Parameters
    pe = PositionalEncoding(d_model=128, dropout=0.1, max_len=140).pe
        # Detector thickness [cm]

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 5))

    for pos, ls in zip([33, 36, 39, 42], ['solid', 'dashed', 'dashdot', 'dotted']):
        ax.plot(np.arange(140), pe[0, :, pos], label=pos, linestyle=ls, linewidth=2)

    ax.set_xlabel("Row (i)")
    ax.set_ylabel("Positional encoding")
    fig.legend(title="Col (j)", bbox_to_anchor=(1.03, 0.7))

    # ax.set_title("Sinusoidal positional encoding")    
    ax.grid(True, which='both', linestyle='--', linewidth=0.3)


    if save_path:
        os.makedirs(save_path, exist_ok=True)
        plt.savefig(os.path.join(save_path, "Positional_encoding.png"), bbox_inches="tight")

    return ax



# ---------------------------------------------------------------------
# Main wrapper
# ---------------------------------------------------------------------

def main():
    """Run the script to generate and save the figure."""
    set_style()
    plot_positional_encoding(save_path="thesis_figures/output/04_transformer/")


if __name__ == "__main__":
    main()