import numpy as np
from typing import Union
from scipy.special import erf as scipy_erf, wofz

ArrayLike = Union[np.ndarray, float]

def gaussian(x: ArrayLike, A: float, x0: float, sigma: float) -> ArrayLike:
    """
    Standard Gaussian function.

    Parameters
    ----------
    x : array-like
        Input x values.
    A : float
        Amplitude.
    x0 : float
        Center position.
    sigma : float
        Standard deviation.

    Returns
    -------
    array-like
        Gaussian profile.
    """
    return A * np.exp(-((x - x0)**2) / (2 * sigma**2))


def lorentzian(x: ArrayLike, A: float, x0: float, gamma: float) -> ArrayLike:
    """
    Lorentzian function.

    Parameters
    ----------
    x : array-like
        Input x values.
    A : float
        Amplitude.
    x0 : float
        Center position.
    gamma : float
        Half-width at half-maximum (HWHM).

    Returns
    -------
    array-like
        Lorentzian profile.
    """
    return A * gamma**2 / ((x - x0)**2 + gamma**2)


def erf(x: ArrayLike, A: float, x0: float, sigma: float) -> ArrayLike:
    """
    Error function (erf) shaped step function.

    Parameters
    ----------
    x : array-like
        Input x values.
    A : float
        Amplitude (scaling of step).
    x0 : float
        Step center position.
    sigma : float
        Slope scale parameter (related to width).

    Returns
    -------
    array-like
        Step-like function based on erf.
    """
    return A * (1 + scipy_erf((x - x0) / (sigma * np.sqrt(2))))


def voigt(x: ArrayLike, A: float, x0: float, sigma: float, gamma: float) -> ArrayLike:
    """
    Voigt profile (Gaussian + Lorentzian convolution).

    Parameters
    ----------
    x : array-like
        Input x values.
    A : float
        Amplitude.
    x0 : float
        Center.
    sigma : float
        Gaussian standard deviation.
    gamma : float
        Lorentzian half-width at half-maximum (HWHM).

    Returns
    -------
    array-like
        Voigt profile.
    """
    z = ((x - x0) + 1j * gamma) / (sigma * np.sqrt(2))
    return A * np.real(wofz(z)) / (sigma * np.sqrt(2 * np.pi))


def gaussian_fwhm(x: ArrayLike, A: float, x0: float, fwhm: float) -> ArrayLike:
    """
    Gaussian function using FWHM parameter.

    Parameters
    ----------
    x : array-like
        Input x values.
    A : float
        Amplitude.
    x0 : float
        Center.
    fwhm : float
        Full width at half maximum.

    Returns
    -------
    array-like
        Gaussian profile.
    """
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return A * np.exp(-((x - x0)**2) / (2 * sigma**2))


def lorentzian_fwhm(x: ArrayLike, A: float, x0: float, fwhm: float) -> ArrayLike:
    """
    Lorentzian function using FWHM parameter.

    Parameters
    ----------
    x : array-like
        Input x values.
    A : float
        Amplitude.
    x0 : float
        Center.
    fwhm : float
        Full width at half maximum.

    Returns
    -------
    array-like
        Lorentzian profile.
    """
    gamma = fwhm / 2
    return A * gamma**2 / ((x - x0)**2 + gamma**2)


def erf_fwhm(x: ArrayLike, A: float, x0: float, fwhm: float) -> ArrayLike:
    """
    Error function shaped step using FWHM-style width parameter.

    Parameters
    ----------
    x : array-like
        Input x values.
    A : float
        Amplitude.
    x0 : float
        Step center.
    fwhm : float
        Effective transition width (FWHM-style).

    Returns
    -------
    array-like
        Step-like function.
    """
    alpha = 2 * np.sqrt(np.log(2)) / fwhm
    return A * (1 + scipy_erf(alpha * (x - x0)))


def voigt_fwhm(x: ArrayLike, A: float, x0: float, fwhm_g: float, fwhm_l: float) -> ArrayLike:
    """
    Voigt profile using FWHM parameters for Gaussian and Lorentzian.

    Parameters
    ----------
    x : array-like
        Input x values.
    A : float
        Amplitude.
    x0 : float
        Center.
    fwhm_g : float
        Gaussian FWHM.
    fwhm_l : float
        Lorentzian FWHM.

    Returns
    -------
    array-like
        Voigt profile.
    """
    sigma = fwhm_g / (2 * np.sqrt(2 * np.log(2)))
    gamma = fwhm_l / 2
    z = ((x - x0) + 1j * gamma) / (sigma * np.sqrt(2))
    return A * np.real(wofz(z)) / (sigma * np.sqrt(2 * np.pi))
