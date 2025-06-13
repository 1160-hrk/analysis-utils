from importlib.resources import files

try:
    version_file = files(__package__).joinpath("../../../VERSION.txt")
    __version__ = version_file.read_text(encoding="utf-8").strip()
except Exception:
    __version__ = "unknown"



from .fft_utils import (
    fft_with_freq,
    fft_positive_freq,
    ifft_with_time,
    zero_padding,
    zero_filling,
    spectrogram,
    spectrogram_fast,
    spectrogram_vectorized,
    spectrogram_scipy
)

from .utils import (
    get_ind,
    get_subarray_1D,
    get_subarray_2D,
    get_ind_max_xi_xf,
    get_ind_xi_xf,
    get_inds_peak_xi_xf,
    fitting_w_range,
    ndarray_from_txtfile,
    ndarray_from_csvfile,
    pickle_dump,
    pickle_load
)


from .funcs import (
    gaussian,
    lorentzian,
    erf,
    voigt,
    gaussian_fwhm,
    lorentzian_fwhm,
    erf_fwhm,
    voigt_fwhm
)

from . import constants as C  # constants モジュールは名前空間ごと "C" で公開
from .signal1d import SignalData
from .map2d import Map2D

__all__ = [
    # fft_utils
    "fft_with_freq", "fft_positive_freq", "ifft_with_time",
    "zero_padding", "zero_filling",
    "spectrogram", "spectrogram_fast", "spectrogram_loop_fast",
    "spectrogram_vectorized", "spectrogram_scipy",

    # utils
    "get_ind", "get_subarray_1D", "get_subarray_2D",
    "get_ind_max_xi_xf", "get_ind_xi_xf", "get_inds_peak_xi_xf",
    "fitting_w_range",
    "ndarray_from_txtfile", "ndarray_from_csvfile",
    "pickle_dump", "pickle_load",

    # data
    "SignalData",

    # funcs
    "gaussian", "lorentzian", "erf", "voigt",
    "gaussian_fwhm", "lorentzian_fwhm", "erf_fwhm", "voigt_fwhm",

    # constants as namespace
    "C",
    
    # signal1d
    "SignalData",
    
    # map2d
    "Map2D",
    
]
