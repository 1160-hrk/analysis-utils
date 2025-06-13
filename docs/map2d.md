# map2d.py

## 🌐 Map2D クラス

`Map2D` は、2次元のスキャン・スペクトルマップ・時間周波数マップなどの信号を扱うクラスです。座標指定による部分抽出やプロット、積分・ピークトラッキングなどの機能を備えています。

---

## 🧱 初期化

```python
Map2D(x: np.ndarray, y: np.ndarray, z: np.ndarray)
```

- `x`: 横軸配列（長さ N）
- `y`: 縦軸配列（長さ M）
- `z`: 信号マップ（shape = (M, N)）

---

## 🎯 カーソル選択・抽出

### `set_cursors(xa, xb, ya, yb)`
x, y 両軸に対してカーソル範囲を指定。

### `get_segment()`
現在のカーソルに基づく部分マップ `(xseg, yseg, zseg)` を返します。

---

## 📊 表示・可視化

### `plot(ax=None, cmap='viridis')`
2D カラーマップ（imshow）として表示します。

### `contour_plot(levels=10)`
等高線表示（matplotlib.contour）を行います。

### `extract_line(axis='x', value=...)`
ある固定値（x または y）における1D断面データを抽出します。

---

## ✨ 平滑化・変換・演算

### `smooth(sigma=1.0, inplace=False)`
2Dガウスフィルターで平滑化します。

### `track_peak(axis='x')`
各行または列方向の最大位置をピークトラッキング。

---

## 🧮 統計・演算

### `integrate_axis(axis='x', method='trapz')`
各行または列方向で積分。結果は `(coord, values)` 形式で返されます。

### `project_axis(axis='x', method='mean')`
各軸に沿って平均・最大・最小を投影します（`mean`, `max`, `min`）。

---

## 💾 データの保存と読み込み

### `to_csv(path)`
マップデータを CSV 形式で保存（x, y, z 全てを含む）。

### `from_csv(path)`
CSV 形式から `Map2D` オブジェクトを生成します。

---

## ✅ 使用例

```python
from pyana.map2d import Map2D

x = np.linspace(0, 10, 100)
y = np.linspace(0, 5, 50)
z = np.random.rand(50, 100)

mapdata = Map2D(x, y, z)
mapdata.set_cursors(2, 6, 1, 4)
mapdata.plot()
```

---

## 🔗 関連機能

- `pyana.signal1d.SignalData`: 1次元信号の解析用
- `pyana.fft_utils`: フーリエ解析・スペクトログラム処理