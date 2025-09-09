"""
Math utilities for LEGEND-200 analysis.

Provides various math functions used for the analysis, 
e.g. gaussian fit functions, activation functions, 
linear functions, etc.
"""

# --- Standard library ---

# # --- Third party library ---
import numpy as np
from scipy.special import erfc

# --- Project modules ---



def gauss_bkg(x: np.ndarray | float,
              A: float,
              mu: float,
              sig: float,
              a: float,
              b: float,
              d: float) -> np.ndarray | float:
    """
    Gaussian + linear background + step function.

    Parameters
    ----------
    x : array_like
        Independent variable.
    A : float
        Gaussian amplitude.
    mu : float
        Gaussian mean.
    sig : float
        Gaussian standard deviation.
    a : float
        Slope of linear background.
    b : float
        Intercept of linear background.
    d : float
        Step amplitude.

    Returns
    -------
    array_like
        Evaluated model values.
    """
    gauss = A * np.exp(-((x - mu) ** 2) / (2 * sig**2))
    linear = a * x + b
    step_term = 0.5 * d * erfc((x - mu) / (np.sqrt(2) * sig))
    return gauss + linear + step_term


def gauss(x: np.ndarray | float, A: float, mu: float, sig: float) -> np.ndarray | float:
    """
    Pure Gaussian function.

    Parameters
    ----------
    x : array_like
        Independent variable.
    A : float
        Gaussian amplitude.
    mu : float
        Gaussian mean.
    sig : float
        Gaussian standard deviation.

    Returns
    -------
    array_like
        Evaluated Gaussian.
    """
    return A * np.exp(-((x - mu) ** 2) / (2 * sig**2))


def linear(x: np.ndarray | float, a: float, b: float) -> np.ndarray | float:
    """
    Linear function y = a*x + b.

    Parameters
    ----------
    x : array_like
        Independent variable.
    a : float
        Slope.
    b : float
        Intercept.

    Returns
    -------
    array_like
        Evaluated line.
    """
    return a * x + b


def step(x: np.ndarray | float, mu: float, sig: float, d: float) -> np.ndarray | float:
    """
    Step function modeled with error function complement.

    Parameters
    ----------
    x : array_like
        Independent variable.
    mu : float
        Center of the step.
    sig : float
        Controls smoothness of the step edge.
    d : float
        Step amplitude.

    Returns
    -------
    array_like
        Evaluated step.
    """
    return 0.5 * d * erfc((x - mu) / (np.sqrt(2) * sig))


def linear_with_err(x: float,
                    fit_param: tuple[float, float],
                    covariance_matrix: np.ndarray) -> tuple[float, float]:
    """
    Linear function with error propagation.

    Parameters
    ----------
    x : float
        Independent variable (scalar).
    fit_param : (float, float)
        Fitted parameters (a, b).
    covariance_matrix : (2,2) ndarray
        Covariance matrix of the fit parameters.

    Returns
    -------
    y : float
        Evaluated linear function at x.
    y_err : float
        Propagated uncertainty at x.
    """
    y = linear(x, *fit_param)

    # Jacobian of f(x, a, b) = a*x + b wrt [a, b]
    J = np.array([x, 1.0])

    y_err = float(np.sqrt(J @ covariance_matrix @ J.T))
    return y, y_err



