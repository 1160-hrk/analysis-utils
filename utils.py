from typing import Any, Callable, Optional, Union, Tuple, List
import csv
import math as mt
import pickle
import numpy as np
from numpy.typing import NDArray
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.signal import find_peaks


def get_ind(array: NDArray[np.float64], values: Union[List[float], NDArray[np.float64]]) -> NDArray[np.int64]:
    """
    Find the indices of the nearest values in the array for each value in values.

    Parameters
    ----------
    array : NDArray[np.float64]
        The input array to search within.
    values : list[float] or NDArray[np.float64]
        Values for which the nearest indices in the array are to be found.

    Returns
    -------
    NDArray[np.int64]
        Indices of the closest values in `array` corresponding to each entry in `values`.
    """
    array = np.asarray(array)
    values = np.asarray(values)
    sorted_indices = np.argsort(array)
    insertion_points = np.searchsorted(array[sorted_indices], values)
    insertion_points = np.clip(insertion_points, 1, len(array) - 1)
    left = sorted_indices[insertion_points - 1]
    right = sorted_indices[insertion_points]
    indices = np.where(np.abs(values - array[left]) <= np.abs(values - array[right]), left, right)
    return indices


def get_subarray_1D(x: NDArray[np.float64], xi: Optional[float] = None, xf: Optional[float] = None,
                    return_ind: bool = True) -> Union[NDArray[np.float64], Tuple[NDArray[np.float64], List[int]]]:
    """
    Extract a subarray from a 1D array within the specified bounds.

    Parameters
    ----------
    x : NDArray[np.float64]
        Input array from which to extract a subarray.
    xi : float, optional
        Start value for the subarray.
    xf : float, optional
        End value for the subarray.
    return_ind : bool, default=False
        If True, also return the index range [i, f].

    Returns
    -------
    NDArray[np.float64] or (NDArray[np.float64], list[int])
        The extracted subarray, optionally with the index range.
    """
    if xi is None:
        i = 0
    else:
        i = int(get_ind(x, [xi])[0])
    if xf is None:
        xx = x[i:]
        f = len(x)
    else:
        f = int(get_ind(x, [xf])[0]) + 1
        xx = x[i:f]
    if return_ind:
        return xx, [i, f]
    else:
        return xx


def get_subarray_2D(x: NDArray[np.float64], y: NDArray[np.float64], xi: Optional[float] = None,
                    xf: Optional[float] = None, return_ind: bool = True) -> Union[
                        Tuple[NDArray[np.float64], NDArray[np.float64]],
                        Tuple[NDArray[np.float64], NDArray[np.float64], List[int]]
                    ]:
    """
    Extract corresponding subarrays from two 1D arrays within a specified range on x.

    Parameters
    ----------
    x : NDArray[np.float64]
        Horizontal axis array.
    y : NDArray[np.float64]
        Vertical axis array.
    xi : float, optional
        Start value on the x-axis.
    xf : float, optional
        End value on the x-axis.
    return_ind : bool, default=False
        Whether to return the index bounds [i, f].

    Returns
    -------
    (x_sub, y_sub) or (x_sub, y_sub, [i, f])
        Subarrays and optionally the index range.
    """
    if xi is None:
        i = 0
    else:
        i = int(get_ind(x, [xi])[0])
    if xf is None:
        xx = x[i:]
        yy = y[i:]
        f = len(x)
    else:
        f = int(get_ind(x, [xf])[0]) + 1
        xx = x[i:f]
        yy = y[i:f]
    if return_ind:
        return xx, yy, [i, f]
    else:
        return xx, yy


def get_ind_max_xi_xf(x: NDArray[np.float64], y: NDArray[np.float64],
                      xi: Optional[float] = None, xf: Optional[float] = None) -> int:
    """
    Find the index of the maximum y-value within a specified x-range.

    Parameters
    ----------
    x : NDArray[np.float64]
        Horizontal data array.
    y : NDArray[np.float64]
        Vertical data array.
    xi : float, optional
        Lower bound on x.
    xf : float, optional
        Upper bound on x.

    Returns
    -------
    int
        Index of the maximum y-value in the range.
    """
    _, yy, ind = get_subarray_2D(x, y, xi=xi, xf=xf, return_ind=True)
    return int(np.argmax(yy) + ind[0])


def get_ind_xi_xf(v: float, x: NDArray[np.float64], y: NDArray[np.float64],
                  xi: Optional[float] = None, xf: Optional[float] = None) -> int:
    """
    Find the index of the y-value closest to a target value v within a range.

    Parameters
    ----------
    v : float
        Target value in y.
    x : NDArray[np.float64]
        Horizontal data.
    y : NDArray[np.float64]
        Vertical data.
    xi : float, optional
        Start x-value.
    xf : float, optional
        End x-value.

    Returns
    -------
    int
        Index of closest y to v within the range.
    """
    _, yy, ind = get_subarray_2D(x, y, xi=xi, xf=xf, return_ind=True)
    ii = int(get_ind(yy, [v])[0])
    return ii + ind[0]


def get_inds_peak_xi_xf(v: NDArray[np.float64], x: NDArray[np.float64], y: NDArray[np.float64],
                        xi: Optional[float] = None, xf: Optional[float] = None,
                        **kwargs: Any) -> NDArray[np.int64]:
    """
    Detect peaks in an array within a specified x-range.

    Parameters
    ----------
    v : NDArray[np.float64]
        Input signal for peak detection.
    x, y : NDArray[np.float64]
        Coordinate arrays (used for slicing range).
    xi, xf : float, optional
        Start and end x-values.
    **kwargs : dict
        Additional arguments to scipy.signal.find_peaks.

    Returns
    -------
    NDArray[np.int64]
        Indices of detected peaks in original array.
    """
    _, _, ind = get_subarray_2D(x, y, xi=xi, xf=xf, return_ind=True)
    ii = np.array(find_peaks(v, **kwargs)[0])
    return ii + ind[0]


def fitting_w_range(
    fit: Callable[[NDArray[np.float64], Any], NDArray[np.float64]],
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    xi: Optional[float] = None,
    xf: Optional[float] = None,
    p0: Optional[Union[List[float], Tuple[float, ...]]] = None,
    bounds: Optional[Tuple[Any, Any]] = None,
    return_x: bool = False,
    return_pcov: bool = False
) -> Any:
    """
    Fit a function to data within a horizontal x-range.

    Parameters
    ----------
    fit : Callable
        Fitting function.
    x, y : NDArray[np.float64]
        Input data.
    xi, xf : float, optional
        Start and end x-values.
    p0 : list or tuple
        Initial fitting parameters.
    bounds : tuple, optional
        Bounds for parameters.
    return_x : bool, default=False
        Return the x-values used in fitting.
    return_pcov : bool, default=False
        Return covariance matrix.

    Returns
    -------
    Depends on return_x and return_pcov.
    """
    xx, yy, _ = get_subarray_2D(x, y, xi=xi, xf=xf, return_ind=True)
    if bounds is None:
        popt, pcov = curve_fit(fit, xx, yy, p0=p0)
    else:
        popt, pcov = curve_fit(fit, xx, yy, p0=p0, bounds=bounds)
    if return_pcov and return_x:
        return popt, pcov, xx
    elif return_pcov:
        return popt, pcov
    elif return_x:
        return popt, xx
    else:
        return popt


def ndarray_from_txtfile(fullname: str, manner: str) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Load 2-column numerical data from a text file.

    Parameters
    ----------
    fullname : str
        Path to the text file.
    manner : str
        Delimiter used in the file.

    Returns
    -------
    a, b : NDArray[np.float64]
        Arrays containing the two columns.
    """
    a = np.array([], dtype=np.float64)
    b = np.array([], dtype=np.float64)
    with open(fullname, 'r') as fid:
        for line in fid:
            currentrow = line.strip().split(manner)
            if currentrow[0] != '':
                cr = np.array(list(map(float, currentrow)))
                a = np.append(a, cr[0])
                b = np.append(b, cr[1])
    return a, b


def ndarray_from_csvfile(path: str) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Load 2-column numerical data from a CSV file.

    Parameters
    ----------
    path : str
        File path.

    Returns
    -------
    a, b : NDArray[np.float64]
        Arrays of the first and second columns.
    """
    a = np.array([], dtype=np.float64)
    b = np.array([], dtype=np.float64)
    with open(path, newline='', encoding='utf-8-sig') as csvfile:
        rows = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in rows:
            if row[0] != '':
                cr = np.array(list(map(float, row)))
                a = np.append(a, cr[0])
                b = np.append(b, cr[1])
    return a, b


def pickle_dump(obj: Any, path: str) -> None:
    """
    Save a Python object using pickle.

    Parameters
    ----------
    obj : Any
        Object to serialize.
    path : str
        Output file path.
    """
    with open(path, mode='wb') as f:
        pickle.dump(obj, f)


def pickle_load(path: str) -> Any:
    """
    Load a pickled Python object.

    Parameters
    ----------
    path : str
        Path to the pickled file.

    Returns
    -------
    Any
        Deserialized object.
    """
    with open(path, mode='rb') as f:
        return pickle.load(f)
