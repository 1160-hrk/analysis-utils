import numpy as np
from typing import Optional, Tuple
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import csv


class Map2D:
    """
    2次元マップデータ（x, y, z）を扱うクラス。
    """

    def __init__(self, x: np.ndarray, y: np.ndarray, z: np.ndarray):
        """
        Parameters
        ----------
        x : np.ndarray
            横軸データ（shape = (N,)）
        y : np.ndarray
            縦軸データ（shape = (M,)）
        z : np.ndarray
            2Dデータ本体（shape = (M, N)） ← z[y, x] 形式
        """
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)

        self.cursor_xa: Optional[float] = None
        self.cursor_xb: Optional[float] = None
        self.cursor_ya: Optional[float] = None
        self.cursor_yb: Optional[float] = None

    def set_cursors(self, xa: float, xb: float, ya: float, yb: float) -> None:
        """
        カーソルで選択範囲を設定（x方向とy方向）

        Parameters
        ----------
        xa, xb : float
            x方向の範囲（任意の順序でOK）
        ya, yb : float
            y方向の範囲（任意の順序でOK）
        """
        self.cursor_xa, self.cursor_xb = sorted([xa, xb])
        self.cursor_ya, self.cursor_yb = sorted([ya, yb])

    def get_segment(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        カーソル範囲の部分マップを取得

        Returns
        -------
        xseg : np.ndarray
            選択された x 軸データ
        yseg : np.ndarray
            選択された y 軸データ
        zseg : np.ndarray
            選択された z マップ（2D）
        """
        if None in (self.cursor_xa, self.cursor_xb, self.cursor_ya, self.cursor_yb):
            raise ValueError("すべてのカーソル（x, y）を設定してください")
        xmask = (self.x >= self.cursor_xa) & (self.x <= self.cursor_xb)
        ymask = (self.y >= self.cursor_ya) & (self.y <= self.cursor_yb)
        return self.x[xmask], self.y[ymask], self.z[np.ix_(ymask, xmask)]

    def plot(self, ax=None, cmap: str = 'viridis', aspect: str = 'auto', **kwargs) -> plt.Axes:
        """
        2Dマップのimshow表示

        Parameters
        ----------
        ax : plt.Axes, optional
            描画に使うAxes（Noneなら新規作成）
        cmap : str
            カラーマップ
        aspect : str
            表示の縦横比（'auto', 'equal'など）
        **kwargs : dict
            plt.imshow に渡す追加引数

        Returns
        -------
        ax : plt.Axes
        """
        if ax is None:
            fig, ax = plt.subplots()
        extent = [self.x[0], self.x[-1], self.y[0], self.y[-1]]
        im = ax.imshow(self.z, extent=extent, origin='lower', cmap=cmap,
                       aspect=aspect, interpolation='nearest', **kwargs)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        plt.colorbar(im, ax=ax)
        return ax

    def contour_plot(self, ax=None, levels: int = 10, cmap: str = 'viridis', **kwargs) -> plt.Axes:
        """
        等高線（contour）プロットを行う

        Parameters
        ----------
        ax : plt.Axes, optional
            描画に使うAxes（Noneなら新規作成）
        levels : int
            等高線の数
        cmap : str
            カラーマップ
        **kwargs : dict
            plt.contour に渡す追加引数

        Returns
        -------
        ax : plt.Axes
        """
        if ax is None:
            fig, ax = plt.subplots()
        X, Y = np.meshgrid(self.x, self.y)
        cs = ax.contour(X, Y, self.z, levels=levels, cmap=cmap, **kwargs)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        plt.colorbar(cs, ax=ax)
        return ax

    def extract_line(self, axis: str = "x", value: float = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        指定軸に沿った1D断面を取得

        Parameters
        ----------
        axis : str
            "x" → y固定、"y" → x固定
        value : float
            断面を取る位置（最近傍を選ぶ）

        Returns
        -------
        coord : np.ndarray
            軸の座標（xまたはy）
        profile : np.ndarray
            z値の断面（xまたはy軸に沿った値）
        """
        if axis == "x":
            idx = np.argmin(np.abs(self.x - value))
            return self.y, self.z[:, idx]
        elif axis == "y":
            idx = np.argmin(np.abs(self.y - value))
            return self.x, self.z[idx, :]
        else:
            raise ValueError("axis must be 'x' or 'y'")

    def smooth(self, sigma: float = 1.0, inplace: bool = False) -> np.ndarray:
        """
        2Dガウシアン平滑化を行う

        Parameters
        ----------
        sigma : float
            平滑化の標準偏差
        inplace : bool
            Trueならself.zを上書き

        Returns
        -------
        z_smooth : np.ndarray
            平滑化後のマップ
        """
        z_smooth = gaussian_filter(self.z, sigma=sigma)
        if inplace:
            self.z = z_smooth
        return z_smooth

    def track_peak(self, axis: str = 'x') -> Tuple[np.ndarray, np.ndarray]:
        """
        各行または各列の最大値位置（ピークトラッキング）

        Parameters
        ----------
        axis : str
            'x' → 行ごと（y方向に対して）最大xを探す
            'y' → 列ごと（x方向に対して）最大yを探す

        Returns
        -------
        coord : np.ndarray
            スキャン軸（行方向なら y、列方向なら x）
        peak_pos : np.ndarray
            ピーク位置（最大値の x または y）
        """
        if axis == "x":
            coord = self.y
            peak_pos = np.array([self.x[np.argmax(row)] for row in self.z])
        elif axis == "y":
            coord = self.x
            peak_pos = np.array([self.y[np.argmax(col)] for col in self.z.T])
        else:
            raise ValueError("axis must be 'x' or 'y'")
        return coord, peak_pos
    
    def integrate_axis(self, axis: str = 'x', method: str = 'trapz') -> Tuple[np.ndarray, np.ndarray]:
        """
        x軸またはy軸方向に沿って積分（各行または列に対して）

        Parameters
        ----------
        axis : str
            'x'：x方向に沿って（各行）積分 → y vs 積分値
            'y'：y方向に沿って（各列）積分 → x vs 積分値
        method : str
            'trapz'（台形積分）または 'simpson'

        Returns
        -------
        coord : np.ndarray
            積分結果に対応する x または y 軸
        values : np.ndarray
            各行または列に対する積分値
        """
        from scipy.integrate import simps, trapz

        if axis == 'x':
            coord = self.y
            if method == 'trapz':
                values = np.array([trapz(row, self.x) for row in self.z])
            elif method == 'simpson':
                values = np.array([simps(row, self.x) for row in self.z])
            else:
                raise ValueError("method must be 'trapz' or 'simpson'")
        elif axis == 'y':
            coord = self.x
            if method == 'trapz':
                values = np.array([trapz(col, self.y) for col in self.z.T])
            elif method == 'simpson':
                values = np.array([simps(col, self.y) for col in self.z.T])
            else:
                raise ValueError("method must be 'trapz' or 'simpson'")
        else:
            raise ValueError("axis must be 'x' or 'y'")

        return coord, values
    
    def project_axis(self, axis: str = 'x', method: str = 'mean') -> Tuple[np.ndarray, np.ndarray]:
        """
        zマップを指定軸方向に沿って射影（平均、最大、最小など）

        Parameters
        ----------
        axis : str
            'x': x方向に沿って（各行に対して）→ y vs 平均値など
            'y': y方向に沿って（各列に対して）→ x vs 平均値など
        method : str
            'mean'（平均）, 'max'（最大）, 'min'（最小）

        Returns
        -------
        coord : np.ndarray
            x または y 軸の座標
        values : np.ndarray
            各行または列に対する代表値（平均など）
        """
        if method not in ['mean', 'max', 'min']:
            raise ValueError("method must be 'mean', 'max', or 'min'")

        if axis == 'x':
            coord = self.y
            if method == 'mean':
                values = np.mean(self.z, axis=1)
            elif method == 'max':
                values = np.max(self.z, axis=1)
            elif method == 'min':
                values = np.min(self.z, axis=1)
        elif axis == 'y':
            coord = self.x
            if method == 'mean':
                values = np.mean(self.z, axis=0)
            elif method == 'max':
                values = np.max(self.z, axis=0)
            elif method == 'min':
                values = np.min(self.z, axis=0)
        else:
            raise ValueError("axis must be 'x' or 'y'")

        return coord, values



    def to_csv(self, path: str) -> None:
        """
        マップデータをCSV保存（1行目: x軸、1列目: y軸）

        Parameters
        ----------
        path : str
            保存先のファイルパス
        """
        with open(path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([''] + list(self.x))
            for i, row in enumerate(self.z):
                writer.writerow([self.y[i]] + list(row))

    @classmethod
    def from_csv(cls, path: str) -> "Map2D":
        """
        CSVからMap2Dオブジェクトを生成

        Parameters
        ----------
        path : str
            読み込むファイルパス

        Returns
        -------
        Map2D
            インスタンス
        """
        with open(path, newline='') as f:
            reader = list(csv.reader(f))
        x = np.array([float(v) for v in reader[0][1:]])
        y = np.array([float(row[0]) for row in reader[1:]])
        z = np.array([[float(val) for val in row[1:]] for row in reader[1:]])
        return cls(x, y, z)
