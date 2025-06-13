# fft_utils.py

FFT（高速フーリエ変換）関連のユーティリティをまとめたモジュール。  
時系列信号やスペクトルの処理・可視化に用います。

---

## 🔁 基本FFT関数

### `fft_with_freq(x, y)`
x軸と信号yに対してFFTを実行し、対応する周波数軸とフーリエ変換結果を返す。

```python
freq, fft_y = fft_with_freq(t, E)
```

---

### `fft_positive_freq(x, y)`
正の周波数成分のみを抽出したFFT結果を返す。

---

### `ifft_with_time(freq, fft_y)`
逆FFTと同時に時間軸を生成。

---

## ➕ ゼロ埋め

### `zero_padding(x, len0pd, connection='lr')`
信号xの前後（または片方）に指定数のゼロを追加。

- `'lr'`: 左右均等に追加
- `'l'`: 左のみ
- `'r'`: 右のみ

---

### `zero_filling(x, len0fl, connection='lr')`
最終的な長さを指定して、ゼロを追加（`zero_padding` を内部利用）。

---

## 📊 スペクトログラム解析

### `spectrogram(...)`
短時間フーリエ変換により、信号の時間変化を解析。

| 引数 | 説明 |
|------|------|
| `T` | 窓幅（サンプル数または物理長） |
| `unit_T` | `'index'` または `'x'` |
| `window_type` | 窓関数（ハミング、ハニングなど） |
| `step` | スライド間隔 |
| `return_max_index` | 各窓で最大の周波数インデックスも返す |

```python
x_spec, freq, spec = spectrogram(t, E, T=100, unit_T='x')
```

---

### `spectrogram_fast(...)`
`numpy`のプリアロケーションで高速化されたバージョン。

---

### `spectrogram_vectorized(...)`
より高速なベクトル化実装（forループを排除）。

---

### `spectrogram_scipy(...)`
scipy.signal.stftを利用した高速化。

---

## 🧪 使用例（GDDのあるチャープパルス）

```python
from pyana.fft_utils import fft_positive_freq, spectrogram_fast
freq, fft_y = fft_positive_freq(t, E)
x_spec, freq, spec = spectrogram_fast(t, E, T=1000, unit_T='x')
```

![Spectrogram](path/to/spectrogram_plot.png)

---

## ⚠️ 注意点

- 窓サイズ `T` は信号長より短く設定すること。
- `spectrogram_fast` / `vectorized` ではステップ数が多すぎると時間がかかる。
- `fft_with_freq` と `ifft_with_time` は信号解析で逆変換や周波数マッチングに便利。

---

## 📚 関連項目

- [`signal1d.py`](./signal1d.md): FFTを利用した解析をラップしたクラス
- [`func.py`](./func.md): ガウシアンなどのフィッティング関数群