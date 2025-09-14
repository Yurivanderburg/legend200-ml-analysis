"""
Utilities for creating figures for the thesis. 
"""

# --- Standard library ---


# --- Third party ---
import matplotlib as mpl
import matplotlib.pyplot as plt

# --- Project modules ---



def set_style():
    """Apply global matplotlib style for thesis figures."""
    mpl.rcParams.update({
        "figure.figsize": (14, 8),
        "font.size": 16,
        "axes.labelsize": 14,
        "lines.linewidth": 1,
        "image.cmap": "magma",
        "axes.titlesize": 16,
        "legend.fontsize": 12,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "figure.dpi": 400,
    })