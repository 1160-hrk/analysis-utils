# func.py

## 📐 概要

このモジュールは、フィッティングでよく使われる関数（ガウシアン・ローレンツ・Voigt・エラーファンクションなど）を提供します。  
それぞれに標準偏差またはFWHMベースのバージョンが含まれています。

---

## 🌟 標準プロファイル関数

### `gaussian(x, A, x0, sigma)`
ガウシアン関数（正規分布）

```python
from pyana.func import gaussian
y = gaussian(x, A=1.0, x0=0.0, sigma=1.0)
```

---

### `lorentzian(x, A, x0, gamma)`
ローレンツ関数

---

### `erf(x, A, x0, sigma)`
エラーファンクション型ステップ

---

### `voigt(x, A, x0, sigma, gamma)`
Voigtプロファイル（ガウシアン＋ローレンツ畳み込み）

---

## 🧾 FWHMベース関数

FWHM（Full Width at Half Maximum：半値全幅）でパラメータを指定できる関数。

### `gaussian_fwhm(x, A, x0, fwhm)`
FWHM指定のガウシアン関数

---

### `lorentzian_fwhm(x, A, x0, fwhm)`
FWHM指定のローレンツ関数

---

### `erf_fwhm(x, A, x0, fwhm)`
FWHM風の傾き制御によるエラーファンクション

---

### `voigt_fwhm(x, A, x0, fwhm_g, fwhm_l)`
FWHM指定のVoigtプロファイル（ガウスとローレンツを独立指定）

---

## 📌 使用例

```python
import numpy as np
import matplotlib.pyplot as plt
from pyana.func import gaussian_fwhm, voigt_fwhm

x = np.linspace(-5, 5, 500)
y1 = gaussian_fwhm(x, A=1.0, x0=0.0, fwhm=1.0)
y2 = voigt_fwhm(x, A=1.0, x0=0.0, fwhm_g=1.0, fwhm_l=0.5)

plt.plot(x, y1, label='Gaussian')
plt.plot(x, y2, label='Voigt')
plt.legend()
plt.show()
```

---

## 📚 関連モジュール

- [`signal1d.py`](./signal1d.md): フィッティングを含む信号処理用クラス
- [`utils.py`](./utils.md): fitting_w_range など汎用フィッティング関数群
