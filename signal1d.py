import numpy as np
from typing import Optional, Tuple, Callable, Union
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, detrend, savgol_filter
from scipy.integrate import simpson, trapezoid
from scipy.fft import fft, fftfreq
from scipy.interpolate import interp1d

class SignalData:
    """
    1次元の信号データとカーソル操作・解析を扱うクラス。
    """

    def __init__(self, x: np.ndarray, y: np.ndarray):
        """
        Parameters
        ----------
        x : np.ndarray
            X軸のデータ（時間、周波数など）
        y : np.ndarray
            Y軸の信号データ
        """
        self.x = np.array(x)
        self.y = np.array(y)
        self.cursor_a: Optional[float] = None
        self.cursor_b: Optional[float] = None

    def set_cursors(self, a: float, b: float):
        """カーソル A, B の範囲を設定"""
        self.cursor_a = min(a, b)
        self.cursor_b = max(a, b)

    def get_segment(self) -> Tuple[np.ndarray, np.ndarray]:
        """カーソル範囲内の x, y サブデータを返す"""
        if self.cursor_a is None or self.cursor_b is None:
            raise ValueError("カーソルAとBが設定されていません")
        mask = (self.x >= self.cursor_a) & (self.x <= self.cursor_b)
        return self.x[mask], self.y[mask]

    def plot(self, ax=None):
        """簡易プロット（matplotlibが必要）"""
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(self.x, self.y, label='Signal')
        if self.cursor_a is not None and self.cursor_b is not None:
            ax.axvline(self.cursor_a, color='r', linestyle='--', label='Cursor A')
            ax.axvline(self.cursor_b, color='b', linestyle='--', label='Cursor B')
        ax.legend()
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        return ax

    def fit_segment(
        self,
        func: Callable,
        p0: Optional[Tuple] = None,
        bounds: Optional[Tuple] = (-np.inf, np.inf),
        return_curve: bool = False
    ) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """
        カーソル範囲内のデータに対してフィッティングを行う。

        Parameters
        ----------
        func : Callable
            フィッティング関数（x, *params）
        p0 : tuple, optional
            初期推定値
        bounds : 2-tuple, optional
            パラメータの制約
        return_curve : bool
            フィット曲線も返すかどうか

        Returns
        -------
        popt : np.ndarray
            フィット後のパラメータ
        pcov : np.ndarray
            共分散行列（return_curve=True のときのみ）
        y_fit : np.ndarray
            フィット曲線（return_curve=True のときのみ）
        """
        xseg, yseg = self.get_segment()
        popt, pcov = curve_fit(func, xseg, yseg, p0=p0, bounds=bounds)
        if return_curve:
            y_fit = func(xseg, *popt)
            return popt, pcov, y_fit
        return popt

    def fft_segment(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        カーソル範囲内のFFTと対応する周波数配列を返す。

        Returns
        -------
        freq : np.ndarray
            周波数軸
        fft_y : np.ndarray
            フーリエ変換された信号（複素数）
        """
        xseg, yseg = self.get_segment()
        dx = xseg[1] - xseg[0]
        N = len(xseg)
        freq = fftfreq(N, d=dx)
        fft_y = fft(yseg)
        return freq, fft_y
    
    def find_extrema_segment(self) -> dict:
        """
        カーソル範囲内の最大・最小値のインデックスと値を返す。

        Returns
        -------
        dict
            {
                'max_index': int,
                'max_value': float,
                'max_x': float,
                'min_index': int,
                'min_value': float,
                'min_x': float
            }
        """
        xseg, yseg = self.get_segment()
        i_max = np.argmax(yseg)
        i_min = np.argmin(yseg)
        return {
            "max_index": i_max,
            "max_value": yseg[i_max],
            "max_x": xseg[i_max],
            "min_index": i_min,
            "min_value": yseg[i_min],
            "min_x": xseg[i_min],
        }


    def find_peaks_segment(self, **kwargs) -> np.ndarray:
        """
        カーソル範囲内でピーク検出を行う。

        Returns
        -------
        peak_indices : np.ndarray
            セグメント内のピークのインデックス
        """
        _, yseg = self.get_segment()
        peaks, _ = find_peaks(yseg, **kwargs)
        return peaks

    def integrate_segment(self, method: str = "simpson") -> float:
        """
        カーソル範囲内の積分を行う。

        Parameters
        ----------
        method : str
            "simpson" or "trapezoidal"

        Returns
        -------
        area : float
            積分値
        """
        xseg, yseg = self.get_segment()
        if method == "simpson":
            return simpson(yseg, xseg)
        elif method == "trapezoidal":
            return trapezoid(yseg, xseg)
        else:
            raise ValueError("method must be 'simpson' or 'trapezoidal'")
    

    def detrend_segment(self, inplace: bool = False) -> np.ndarray:
        """
        カーソル範囲内の信号から線形トレンドを除去。

        Parameters
        ----------
        inplace : bool, default=False
            True の場合、元のデータを置き換える。

        Returns
        -------
        y_detrended : np.ndarray
            トレンド除去後の信号（カーソル範囲）
        """
        xseg, yseg = self.get_segment()
        y_detrended = detrend(yseg)
        if inplace:
            mask = (self.x >= self.cursor_a) & (self.x <= self.cursor_b)
            self.y[mask] = y_detrended
        return y_detrended

    def baseline_subtract_segment(
        self, method: str = "mean", value: Optional[float] = None, inplace: bool = False
    ) -> np.ndarray:
        """
        ベースラインを引く（カーソル範囲内）処理。

        Parameters
        ----------
        method : str
            "mean" or "custom"
        value : float, optional
            method="custom" の場合の基準値
        inplace : bool, default=False
            True の場合、元のデータを置き換える

        Returns
        -------
        y_corrected : np.ndarray
            ベースライン補正後の信号（カーソル範囲）
        """
        xseg, yseg = self.get_segment()

        if method == "mean":
            baseline = np.mean(yseg)
        elif method == "custom":
            if value is None:
                raise ValueError("method='custom' のときは value を指定してください")
            baseline = value
        else:
            raise ValueError("method must be 'mean' or 'custom'")

        y_corrected = yseg - baseline
        if inplace:
            mask = (self.x >= self.cursor_a) & (self.x <= self.cursor_b)
            self.y[mask] = y_corrected
        return y_corrected

    def smooth_segment(
        self,
        mode: str = "moving_average",
        window_size: int = 5,
        polyorder: int = 2,
        inplace: bool = False
    ) -> np.ndarray:
        """
        カーソル範囲内の信号を平滑化する。

        Parameters
        ----------
        mode : str
            'moving_average' または 'savgol'
        window_size : int
            窓幅（奇数推奨）
        polyorder : int
            Savitzky-Golay の場合の多項式次数
        inplace : bool
            元データに適用するか

        Returns
        -------
        y_smooth : np.ndarray
            平滑化された y データ（セグメント範囲）
        """
        xseg, yseg = self.get_segment()
        if mode == "moving_average":
            kernel = np.ones(window_size) / window_size
            y_smooth = np.convolve(yseg, kernel, mode="same")
        elif mode == "savgol":
            if window_size % 2 == 0:
                window_size += 1  # 必要に応じて奇数化
            y_smooth = savgol_filter(yseg, window_length=window_size, polyorder=polyorder)
        else:
            raise ValueError("mode must be 'moving_average' or 'savgol'")

        if inplace:
            mask = (self.x >= self.cursor_a) & (self.x <= self.cursor_b)
            self.y[mask] = y_smooth
        return y_smooth

    def resample_segment(
        self,
        num: Optional[int] = None,
        step: Optional[float] = None,
        method: str = "linear"
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        カーソル範囲内のデータを等間隔に補間し直す。

        Parameters
        ----------
        num : int, optional
            出力点数（優先）
        step : float, optional
            x 軸のステップ幅（num が None のとき有効）
        method : str
            補間法：'linear', 'cubic', 'nearest', 'quadratic' など

        Returns
        -------
        x_new : np.ndarray
            等間隔な x 軸
        y_new : np.ndarray
            補間された y 値
        """
        xseg, yseg = self.get_segment()

        # x 軸の新しい分割点を生成
        if num is not None:
            x_new = np.linspace(xseg[0], xseg[-1], num)
        elif step is not None:
            x_new = np.arange(xseg[0], xseg[-1], step)
        else:
            raise ValueError("num または step のいずれかを指定してください")

        # 補間器で補完
        interpolator = interp1d(xseg, yseg, kind=method, bounds_error=False, fill_value="extrapolate")
        y_new = interpolator(x_new)

        return x_new, y_new


    def to_csv(self, path: str):
        """CSVファイルに保存"""
        np.savetxt(path, np.column_stack((self.x, self.y)), delimiter=',')

    @classmethod
    def from_csv(cls, path: str) -> "SignalData":
        """CSVファイルから読み込み"""
        data = np.loadtxt(path, delimiter=',')
        return cls(data[:, 0], data[:, 1])
