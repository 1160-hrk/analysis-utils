import numpy as np
from typing import Tuple, Union, List
from scipy.signal import stft

def fft_with_freq(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform FFT on a signal and return the corresponding frequency array and FFT result.

    Parameters:
    x : np.ndarray
        Time or spatial axis of the signal.
    y : np.ndarray
        Signal to be transformed.

    Returns:
    freq : np.ndarray
        Frequency components.
    fft_y : np.ndarray
        FFT result of the signal.
    """
    dx = x[1] - x[0]
    N = len(x)
    freq = np.fft.fftfreq(N, d=dx)
    fft_y = np.fft.fft(y)
    return freq, fft_y

def fft_positive_freq(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform FFT and return only the positive frequency components.

    Parameters:
    x : np.ndarray
        Time or spatial axis of the signal.
    y : np.ndarray
        Signal to be transformed.

    Returns:
    freq : np.ndarray
        Positive frequency components.
    fft_y : np.ndarray
        FFT result corresponding to positive frequencies.
    """
    dx = x[1] - x[0]
    N = len(x)
    freq = np.fft.fftfreq(N, d=dx)
    fft_y = np.fft.fft(y)
    mask = freq >= 0
    return freq[mask], fft_y[mask]

def zero_padding(x: np.ndarray, len0pd: int, connection: str = 'lr') -> np.ndarray:
    """
    Pad the input array with zeros.

    Parameters:
    x : np.ndarray
        Input array.
    len0pd : int
        Total number of zeros to pad.
    connection : str
        Where to apply the padding:
        'lr' - split equally on left and right,
        'l'  - all on the left,
        'r'  - all on the right.

    Returns:
    np.ndarray
        Zero-padded array.
    """
    pad_len = len0pd
    if pad_len < 0:
        raise ValueError("len0pd must be non-negative")
    if connection == 'lr':
        pad_left = pad_len // 2
        pad_right = pad_len - pad_left
    elif connection == 'l':
        pad_left = pad_len
        pad_right = 0
    elif connection == 'r':
        pad_left = 0
        pad_right = pad_len
    else:
        raise ValueError("connection must be 'lr', 'l', or 'r'")
    return np.pad(x, (pad_left, pad_right), mode='constant')

def zero_filling(x: np.ndarray, len0fl: int, connection: str = 'lr') -> np.ndarray:
    """
    Zero-fill the input array so that the output has a specified total length.

    Parameters:
    x : np.ndarray
        Input array.
    len0fl : int
        Desired total output length.
    connection : str
        Padding direction: 'lr', 'l', or 'r'.

    Returns:
    np.ndarray
        Zero-filled array of length len0fl.
    """
    if len0fl < len(x):
        raise ValueError("len0fl must be greater than or equal to the length of x")
    len0pd = len0fl - len(x)
    return zero_padding(x, len0pd, connection)

def ifft_with_time(freq: np.ndarray, fft_y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform inverse FFT and return the corresponding time axis.

    Parameters:
    freq : np.ndarray
        Frequency axis.
    fft_y : np.ndarray
        Frequency domain signal.

    Returns:
    time : np.ndarray
        Time axis corresponding to the inverse transform.
    y : np.ndarray
        Signal in time domain after inverse FFT.
    """
    N = len(freq)
    df = freq[1] - freq[0]
    T = 1 / df
    dt = T / N
    time = np.arange(N) * dt
    y = np.fft.ifft(fft_y)
    return time, y

def spectrogram(
    x: np.ndarray,
    y: np.ndarray,
    T: Union[int, float],
    unit_T: str = 'index',
    window_type: str = 'hamming',
    return_max_index: bool = False,
    step: int = 1
) -> Union[Tuple[np.ndarray, np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """
    Compute the spectrogram of signal y using sliding window FFT.

    Parameters:
    x : np.ndarray
        Time axis.
    y : np.ndarray
        Signal to analyze.
    T : int or float
        Window size. If unit_T is 'index', interpreted as number of samples.
        If unit_T is 'x', interpreted as range in units of x.
    unit_T : str, default='index'
        Unit of T. Either 'index' or 'x'.
    window_type : str, default='hamming'
        Type of window function to apply. Options: 'rectangle', 'triangle', 'hamming', 'han', 'blackman'.
    return_max_index : bool, default=False
        Whether to return the frequency index with the maximum amplitude for each window.
    step : int, default=1
        Step size for sliding the window.

    Returns:
    Tuple containing:
    - x_spec : np.ndarray
        Center positions of each window on the x-axis.
    - freq : np.ndarray
        Array of positive frequency components.
    - spec : np.ndarray
        2D spectrogram array of shape (len(freq), number of windows).
    - max_indices : np.ndarray, optional
        Indices of the frequency component with the maximum amplitude per window.
        Only returned if return_max_index is True.
    """
    dx = x[1] - x[0]
    N = len(x)

    if unit_T == 'x':
        T_index = int(np.round(T / dx))
    elif unit_T == 'index':
        T_index = int(T)
    else:
        raise ValueError("unit_T must be 'index' or 'x'")

    if window_type == 'rectangle':
        window = np.ones(T_index)
    elif window_type == 'triangle':
        window = 1 - np.abs((np.arange(T_index) - T_index / 2) / (T_index / 2))
    elif window_type == 'hamming':
        window = np.hamming(T_index)
    elif window_type == 'han':
        window = np.hanning(T_index)
    elif window_type == 'blackman':
        window = np.blackman(T_index)
    else:
        raise ValueError("Unknown window type")

    num_windows = (N - T_index) // step + 1
    spec = []
    max_indices = []
    x_spec = []

    for i in range(num_windows):
        start = i * step
        end = start + T_index
        if end > N:
            break
        segment = y[start:end] * window
        freq, fft_segment = fft_with_freq(x[start:end], segment)
        mask = freq >= 0
        fft_segment = fft_segment[mask]
        freq = freq[mask]
        spec.append(np.abs(fft_segment))
        if return_max_index:
            max_indices.append(np.argmax(np.abs(fft_segment)))
        x_spec.append(x[start + T_index // 2])

    spec = np.array(spec).T
    if return_max_index:
        return np.array(x_spec), freq, spec, np.array(max_indices)
    else:
        return np.array(x_spec), freq, spec

def spectrogram_fast(
    x: np.ndarray,
    y: np.ndarray,
    T: Union[int, float],
    unit_T: str = 'index',
    window_type: str = 'hamming',
    return_max_index: bool = False,
    step: int = 1
) -> Union[Tuple[np.ndarray, np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """
    Fast version of spectrogram using precomputed FFT frequency and preallocated arrays.

    This function computes the spectrogram of a signal using a sliding window approach.
    It applies a window function to each segment of the signal and performs a Fast Fourier Transform (FFT)
    on the windowed segment. Only the non-negative frequency components are kept. The output is a 2D array
    where each column represents the FFT magnitude spectrum at a specific time point.

    Parameters:
    x : np.ndarray
        Time axis (1D array).
    y : np.ndarray
        Signal to analyze (same shape as x).
    T : int or float
        Window size. If unit_T is 'index', interpreted as number of samples.
        If unit_T is 'x', interpreted as range in units of x.
    unit_T : str, default='index'
        Unit of T. Either 'index' or 'x'.
    window_type : str, default='hamming'
        Window function to apply. Options: 'rectangle', 'triangle', 'hamming', 'han', 'blackman'.
    return_max_index : bool, default=False
        Whether to return the frequency index with the maximum amplitude for each window.
    step : int, default=1
        Step size for sliding the window.

    Returns:
    Tuple containing:
    - x_spec : np.ndarray
        Center positions of each window on the x-axis.
    - freq : np.ndarray
        Array of positive frequency components.
    - spec : np.ndarray
        2D spectrogram array of shape (len(freq), number of windows).
    - max_indices : np.ndarray, optional
        Indices of the frequency component with the maximum amplitude per window.
        Only returned if return_max_index is True.
    """
    dx = x[1] - x[0]
    N = len(x)

    if unit_T == 'x':
        T_index = int(np.round(T / dx))
    elif unit_T == 'index':
        T_index = int(T)
    else:
        raise ValueError("unit_T must be 'index' or 'x'")

    if T_index > N:
        raise ValueError("Window size T is larger than input signal.")

    if window_type == 'rectangle':
        window = np.ones(T_index)
    elif window_type == 'triangle':
        window = 1 - np.abs((np.arange(T_index) - T_index / 2) / (T_index / 2))
    elif window_type == 'hamming':
        window = np.hamming(T_index)
    elif window_type == 'han':
        window = np.hanning(T_index)
    elif window_type == 'blackman':
        window = np.blackman(T_index)
    else:
        raise ValueError("Unknown window type")

    freq_full = np.fft.fftfreq(T_index, d=dx)
    mask = freq_full >= 0
    freq = freq_full[mask]
    num_freq = len(freq)

    num_windows = (N - T_index) // step + 1
    spec = np.empty((num_freq, num_windows), dtype=np.float64)
    max_indices = np.empty(num_windows, dtype=int) if return_max_index else None
    x_spec = np.empty(num_windows, dtype=np.float64)

    for i in range(num_windows):
        start = i * step
        end = start + T_index
        segment = y[start:end] * window
        fft_segment = np.fft.fft(segment)[mask]
        abs_fft = np.abs(fft_segment)
        spec[:, i] = abs_fft
        if return_max_index:
            max_indices[i] = np.argmax(abs_fft)
        x_spec[i] = x[start + T_index // 2]

    if return_max_index:
        return x_spec, freq, spec, max_indices
    else:
        return x_spec, freq, spec


def spectrogram_vectorized(
    x: np.ndarray,
    y: np.ndarray,
    T: Union[int, float],
    unit_T: str = 'index',
    window_type: str = 'hamming',
    return_max_index: bool = False,
    step: int = 1
) -> Union[Tuple[np.ndarray, np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """
    Fast version of spectrogram using vectorized batch FFT computation (Numba-free).
    """
    dx = x[1] - x[0]
    N = len(x)

    if unit_T == 'x':
        T_index = int(np.round(T / dx))
    elif unit_T == 'index':
        T_index = int(T)
    else:
        raise ValueError("unit_T must be 'index' or 'x'")

    if T_index > N:
        raise ValueError("Window size T is larger than input signal.")

    if window_type == 'rectangle':
        window = np.ones(T_index)
    elif window_type == 'triangle':
        window = 1 - np.abs((np.arange(T_index) - T_index / 2) / (T_index / 2))
    elif window_type == 'hamming':
        window = np.hamming(T_index)
    elif window_type == 'han':
        window = np.hanning(T_index)
    elif window_type == 'blackman':
        window = np.blackman(T_index)
    else:
        raise ValueError("Unknown window type")

    freq_full = np.fft.fftfreq(T_index, d=dx)
    mask = freq_full >= 0
    freq = freq_full[mask]
    num_freq = len(freq)

    num_windows = (N - T_index) // step + 1
    indices = np.arange(0, num_windows * step, step)[:, None] + np.arange(T_index)
    segments = y[indices] * window  # shape: (num_windows, T_index)

    fft_segments = np.fft.fft(segments, axis=1)[:, mask]  # shape: (num_windows, num_freq)
    abs_fft = np.abs(fft_segments).T  # shape: (num_freq, num_windows)

    x_spec = x[indices[:, T_index // 2]]

    if return_max_index:
        max_indices = np.argmax(abs_fft, axis=0)
        return x_spec, freq, abs_fft, max_indices
    else:
        return x_spec, freq, abs_fft


def spectrogram_scipy(
    x: np.ndarray,
    y: np.ndarray,
    T: Union[int, float],
    unit_T: str = 'index',
    window_type: str = 'hamming',
    return_max_index: bool = False,
    step: int = 1
) -> Union[Tuple[np.ndarray, np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """
    Compute a spectrogram using SciPy's optimized short-time Fourier transform (STFT).

    This function leverages scipy.signal.stft to perform STFT with efficient overlapping
    and windowing, returning the magnitude spectrogram and optionally the index of
    the dominant frequency at each time step.

    Parameters
    ----------
    x : np.ndarray
        Time axis of the input signal.
    y : np.ndarray
        Signal to analyze.
    T : int or float
        Window size. Interpreted as number of samples if unit_T is 'index', or as time range if 'x'.
    unit_T : {'x', 'index'}, default='index'
        Unit of T.
    window_type : str, default='hamming'
        Type of window passed to scipy.signal.stft.
    return_max_index : bool, default=False
        If True, returns the index of the peak amplitude in each time slice.
    step : int, default=1
        Step size for sliding window.

    Returns
    -------
    x_spec : np.ndarray
        Time axis of the spectrogram.
    freq : np.ndarray
        Frequency axis.
    spec : np.ndarray
        Spectrogram magnitude array (frequency x time).
    max_indices : np.ndarray, optional
        Indices of the frequency bin with maximum amplitude at each time point.
    """
    dx = x[1] - x[0]
    fs = 1 / dx

    T_index = int(np.round(T / dx)) if unit_T == 'x' else int(T)
    noverlap = T_index - step

    f, t_spec, Zxx = stft(
        y,
        fs=fs,
        window=window_type,
        nperseg=T_index,
        noverlap=noverlap,
        nfft=T_index,
        return_onesided=True,
        boundary=None,
        padded=False
    )

    x_spec = x[0] + t_spec
    spec = np.abs(Zxx).T  # 転置して (freq x time) → (time x freq) を修正
    f = f[:-1]
    spec = spec[:, :-1]
    
    if return_max_index:
        max_indices = np.argmax(spec.T, axis=0)
        return x_spec, f, spec.T, max_indices
    else:
        return x_spec, f, spec.T


if __name__ == "__main__":
    import time
    import matplotlib.pyplot as plt
     # --- GDD（Group Delay Dispersion）によるチャープ ---
    t = np.linspace(-2000, 2000, 2**15)  # 時間軸 [fs]
    tau0 = 30  # 初期パルス幅 [fs]
    f0 = 0.07  # 中心周波数 [1/fs]
    GDD = -50000  # fs^2（群遅延分散）

    # 初期ガウシアンパルス（時間領域）
    E0 = np.exp(-t**2 / (2 * tau0**2)) * np.cos(2 * np.pi * f0 * t)

    # 周波数領域に変換
    dt = t[1] - t[0]
    freq = np.fft.fftfreq(len(t), d=dt)
    E_freq = np.fft.fft(E0)

    # GDDによる位相（freq>=0: -1jπGDD(f-f0)^2, freq<0: +1jπGDD(f+f0)^2）
    phase = np.where(
        freq >= 0,
        np.exp(-1j * np.pi * GDD * (freq - f0)**2),
        np.exp(+1j * np.pi * GDD * (freq + f0)**2)
    )
    E_freq_chirped = E_freq * phase

    # 時間領域に戻す
    E = np.real(np.fft.ifft(E_freq_chirped))

    # --- FFT ---
    freq_pos, E_fft = fft_positive_freq(t, E)

    # --- スペクトログラム（通常版） ---
    T_window = 1000
    step = 1
    start_time = time.time()
    x_spec1, freq1, spec1, max_idx1 = spectrogram(t, E, T=T_window, unit_T='x', step=step, return_max_index=True)
    elapsed1 = time.time() - start_time
    print(f"spectrogram      : {elapsed1:.4f} sec, shape = {spec1.shape}")

    # --- スペクトログラム（高速版） ---
    start_time = time.time()
    x_spec2, freq2, spec2, max_idx2 = spectrogram_scipy(t, E, T=T_window, unit_T='x', step=step, return_max_index=True)
    elapsed2 = time.time() - start_time
    print(f"spectrogram_fast : {elapsed2:.4f} sec, shape = {spec2.shape}")

    # --- 差の確認 ---
    diff = np.max(np.abs(spec1 - spec2))
    print(f"max(abs(diff))   = {diff:.4e}")

    # %%
    # --- プロット ---
    fig, axs = plt.subplots(3, 1, figsize=(10, 10))

    axs[0].plot(t, E)
    axs[0].set_title('GDD-Applied Gaussian Pulse')
    axs[0].set_xlabel('Time [fs]')
    axs[0].set_ylabel('Amplitude')

    axs[1].plot(freq_pos, np.abs(E_fft))
    axs[1].set_title('FFT Spectrum')
    axs[1].set_xlabel('Frequency [1/fs]')
    axs[1].set_ylabel('Amplitude')

    X, Y = np.meshgrid(x_spec2, freq2)
    cf = axs[2].contourf(X, Y, spec2, levels=100, cmap='viridis')
    axs[2].plot(x_spec2, freq2[max_idx2], color='white', linewidth=1.5, label='Max freq')
    axs[2].set_title('Spectrogram (Fast)')
    axs[2].set_xlabel('Time [fs]')
    axs[2].set_ylabel('Frequency [1/fs]')
    axs[2].legend()
    axs[2].set_ylim(0.05, 0.09)
    fig.colorbar(cf, ax=axs[2])

    plt.tight_layout()
    plt.show()
