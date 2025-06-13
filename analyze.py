# analyze.py

from .signal1d import SignalData
from .map2d import Map2D
from .funcs import gaussian, lorentzian, voigt
from typing import Optional, Tuple, Callable, Union, Literal
import numpy as np

__all__ = [
    # SignalData
    "load_signal", "save_signal", "set_cursors", "get_segment", "plot_signal",
    "fit_segment", "fit_peak", "fft_segment", "find_extrema_segment",
    "find_peaks_segment", "integrate_segment", "detrend_segment",
    "baseline_subtract_segment", "smooth_segment", "resample_segment",

    # Map2D
    "load_map", "save_map", "set_map_cursors", "get_map_segment", "plot_map",
    "contour_map", "extract_map_line", "smooth_map", "track_map_peak",
    "integrate_map_axis", "project_map_axis"
]

# SignalData 関連関数

def load_signal(path: str) -> SignalData:
    """Load 1D signal from CSV."""
    return SignalData.from_csv(path)

def save_signal(signal: SignalData, path: str) -> None:
    """Save 1D signal to CSV."""
    signal.to_csv(path)

def set_cursors(signal: SignalData, a: float, b: float) -> None:
    """Set cursor range for 1D signal."""
    signal.set_cursors(a, b)

def get_segment(signal: SignalData) -> Tuple[np.ndarray, np.ndarray]:
    """Get x, y data in the cursor range."""
    return signal.get_segment()

def plot_signal(signal: SignalData, ax=None):
    """Plot the 1D signal."""
    return signal.plot(ax=ax)

def fit_segment(
    signal: SignalData,
    func: Callable,
    p0: Optional[Tuple] = None,
    bounds: Optional[Tuple] = (-np.inf, np.inf),
    return_curve: bool = False
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Fit signal segment with given function."""
    return signal.fit_segment(func, p0=p0, bounds=bounds, return_curve=return_curve)

def fit_peak(
    signal: SignalData,
    model: Literal["gaussian", "lorentzian", "voigt"] = "gaussian"
) -> np.ndarray:
    """Fit the peak in signal segment using the specified model."""
    from numpy import argmax
    x, y = signal.get_segment()
    model_func = {
        "gaussian": gaussian,
        "lorentzian": lorentzian,
        "voigt": voigt,
    }.get(model)
    if model_func is None:
        raise ValueError(f"Unsupported model: {model}")

    amp0 = y.max()
    x0 = x[argmax(y)]
    fwhm0 = (x[-1] - x[0]) / 4
    p0 = (amp0, x0, fwhm0)
    return signal.fit_segment(model_func, p0=p0)

def fft_segment(signal: SignalData) -> Tuple[np.ndarray, np.ndarray]:
    """Return FFT and frequency array for the signal segment."""
    return signal.fft_segment()

def find_extrema_segment(signal: SignalData) -> dict:
    """Return extrema (min/max) info in the signal segment."""
    return signal.find_extrema_segment()

def find_peaks_segment(signal: SignalData, **kwargs) -> np.ndarray:
    """Find peaks in the signal segment."""
    return signal.find_peaks_segment(**kwargs)

def integrate_segment(signal: SignalData, method: str = "simpson") -> float:
    """Integrate the signal segment using the specified method."""
    return signal.integrate_segment(method=method)

def detrend_segment(signal: SignalData, inplace: bool = False) -> np.ndarray:
    """Remove linear trend from the signal segment."""
    return signal.detrend_segment(inplace=inplace)

def baseline_subtract_segment(
    signal: SignalData,
    method: str = "mean",
    value: Optional[float] = None,
    inplace: bool = False
) -> np.ndarray:
    """Subtract baseline from the signal segment."""
    return signal.baseline_subtract_segment(method=method, value=value, inplace=inplace)

def smooth_segment(
    signal: SignalData,
    mode: str = "moving_average",
    window_size: int = 5,
    polyorder: int = 2,
    inplace: bool = False
) -> np.ndarray:
    """Smooth the signal segment using moving average or Savitzky-Golay."""
    return signal.smooth_segment(mode=mode, window_size=window_size, polyorder=polyorder, inplace=inplace)

def resample_segment(
    signal: SignalData,
    num: Optional[int] = None,
    step: Optional[float] = None,
    method: str = "linear"
) -> Tuple[np.ndarray, np.ndarray]:
    """Resample the signal segment using interpolation."""
    return signal.resample_segment(num=num, step=step, method=method)


# Map2D 関連関数

def load_map(path: str) -> Map2D:
    """Load 2D map from CSV."""
    return Map2D.from_csv(path)

def save_map(map2d: Map2D, path: str) -> None:
    """Save 2D map to CSV."""
    map2d.to_csv(path)

def set_map_cursors(map2d: Map2D, xa: float, xb: float, ya: float, yb: float) -> None:
    """Set cursor range for 2D map."""
    map2d.set_cursors(xa, xb, ya, yb)

def get_map_segment(map2d: Map2D) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Get x, y, z segment within cursor range."""
    return map2d.get_segment()

def plot_map(map2d: Map2D, ax=None, cmap: str = 'viridis', aspect: str = 'auto', **kwargs):
    """Plot 2D map with imshow."""
    return map2d.plot(ax=ax, cmap=cmap, aspect=aspect, **kwargs)

def contour_map(map2d: Map2D, ax=None, levels: int = 10, cmap: str = 'viridis', **kwargs):
    """Plot 2D map with contour lines."""
    return map2d.contour_plot(ax=ax, levels=levels, cmap=cmap, **kwargs)

def extract_map_line(map2d: Map2D, axis: str = "x", value: float = None) -> Tuple[np.ndarray, np.ndarray]:
    """Extract 1D cross-section along given axis at specified coordinate."""
    return map2d.extract_line(axis=axis, value=value)

def smooth_map(map2d: Map2D, sigma: float = 1.0, inplace: bool = False) -> np.ndarray:
    """Smooth 2D map with Gaussian filter."""
    return map2d.smooth(sigma=sigma, inplace=inplace)

def track_map_peak(map2d: Map2D, axis: str = 'x') -> Tuple[np.ndarray, np.ndarray]:
    """Track peak positions along specified axis."""
    return map2d.track_peak(axis=axis)

def integrate_map_axis(map2d: Map2D, axis: str = 'x', method: str = 'trapz') -> Tuple[np.ndarray, np.ndarray]:
    """Integrate map along specified axis using the given method."""
    return map2d.integrate_axis(axis=axis, method=method)

def project_map_axis(map2d: Map2D, axis: str = 'x', method: str = 'mean') -> Tuple[np.ndarray, np.ndarray]:
    """Project map along specified axis (mean, max, etc)."""
    return map2d.project_axis(axis=axis, method=method)
