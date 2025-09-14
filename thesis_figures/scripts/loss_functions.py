"""
Figure: Loss functions for machine learning

This script creates a visualization of different loss functions 
converging an ill-conditioned quadratic function. 
"""

# --- Standard library ---
import numpy as np
import os

# --- Third party ---
import torch
import matplotlib.pyplot as plt

# --- Project modules ---
from utils import set_style



# ---------------------------------------------------------------------
# Figure-specific plotting function
# ---------------------------------------------------------------------
def quad_loss(w: torch.tensor,  lambda1: float = 1.0, lambda2: float = 0.05):
    """
    Quadratic loss function, preferrably very ill-conditioned.
    """
    x, y = w  # w is a torch tensor of shape (2,)
    return 0.5*(lambda1*x**2 + lambda2*y**2)



def run_optimizer(optimizer_cls, w_start, steps, **opt_kwargs):
    """
    Optimizer function.
    """
    w = torch.nn.Parameter(torch.tensor(w_start, dtype=torch.float32))
    opt = optimizer_cls([w], **opt_kwargs)
    path = [w.detach().numpy().copy()]
    for _ in range(steps):
        opt.zero_grad(set_to_none=True)
        loss = quad_loss(w)
        loss.backward()
        opt.step()
        path.append(w.detach().numpy().copy())
    return np.array(path)


def plot_charge_collection(ax=None, save_path=None):
    """
    Generate the charge collection figure.

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
    # Settings
    w_start = [1.5, 1.5]
    steps = 50
    gd_path   = run_optimizer(torch.optim.SGD,  w_start, steps, lr=0.05)               # SGD
    mom_path  = run_optimizer(torch.optim.SGD,  w_start, steps, lr=0.05, momentum=0.9) # Momentum
    adam_path = run_optimizer(torch.optim.Adam, w_start, steps, lr=0.05)               # Adam

    # Contour grid
    x = np.linspace(-2, 2, 500)
    y = np.linspace(-2, 2, 500)
    X, Y = np.meshgrid(x, y)
    Z = quad_loss([X,Y])

    # Plot
    if ax is None:
        fig, axes = plt.subplots(figsize=(9,7.5))

    levels = np.linspace(0, np.sqrt(Z).max(), 20)**2  # expand near min
    cs = plt.contour(X, Y, Z, levels=levels, alpha=0.5)
    plt.clabel(cs, inline=True, fontsize=8)

    plt.plot(gd_path[:,0],   gd_path[:,1],   'o-', ms=3, label='SGD')
    plt.plot(mom_path[:,0],  mom_path[:,1],  'o-', ms=3, label='Momentum')
    plt.plot(adam_path[:,0], adam_path[:,1], 'o-', ms=3, label='Adam')

    plt.scatter(*w_start, marker='x', s=100, label='Start', color='C3')
    plt.scatter(0, 0, marker='*', s=150, label='Minimum', color='C3')

    plt.xlabel(r'$\theta^{(0)}$')
    plt.ylabel(r'$\theta^{(1)}$')
    plt.legend(loc='center', bbox_to_anchor=(0.5, -0.15), fontsize=14, ncol = 5)
    plt.grid(True)
    plt.tight_layout()



    if save_path:
        os.makedirs(save_path, exist_ok=True)
        plt.savefig(os.path.join(save_path, "Loss_function.png"), bbox_inches="tight")

    return ax



# ---------------------------------------------------------------------
# Main wrapper
# ---------------------------------------------------------------------

def main():
    """Run the script to generate and save the figure."""
    set_style()
    plot_charge_collection(save_path="thesis_figures/output/04_transformer/")


if __name__ == "__main__":
    main()