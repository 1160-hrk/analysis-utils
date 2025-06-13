# 🧪 pyana

Python-based analysis toolkit for scientific data processing.

---

## 🔰 Overview

`pyana` は、実験データやシミュレーションデータの解析に役立つ汎用的なツール群を集めた Python ライブラリです。スペクトル解析、2次元マップ解析、物理定数、ユーティリティ関数などを包括的に提供します。

- 🧮 FFTとスペクトログラム
- 📈 ピーク検出・フィッティング・積分
- 📊 信号データ (`SignalData`) / 2Dマップ (`Map2D`) の扱い
- 📚 高度な可視化と前処理支援

---

## 📦 Installation

### ✅ Standard installation

You can install `pyana` directly from the GitHub repository using pip:

```bash
pip install git+https://github.com/yourname/pyana.git
```

> Note: This project uses a `pyproject.toml`-based build. No `setup.py` is required.

---

### 🧪 Editable (development) installation

To install the project in "editable" mode (for local development):

```bash
git clone https://github.com/yourname/pyana.git
cd pyana
pip install -e .
```

---

## 📁 Project structure (src layout)

```
pyana/
├── pyproject.toml
├── src/
│   └── pyana/
│       ├── __init__.py
│       ├── utils.py
│       ├── map2d.py
│       └── ...
├── README.md
├── VERSION.txt
└── docs/
```

---

## 📦 Dependencies

This project depends on:

- `numpy`
- `matplotlib`

These will be automatically installed via `pyproject.toml`.

---

## 🧠 Python version

Requires **Python 3.8 or higher**.

---

## 🗂 Module Structure

| モジュール名         | 説明                                                                 |
|----------------------|----------------------------------------------------------------------|
| `utils.py`           | 汎用ユーティリティ（部分抽出・ファイル入出力・フィッティングなど） |
| `fft_utils.py`       | FFT・Zero padding・スペクトログラムなど                              |
| `constants.py`       | 物理定数と単位の辞書                                                 |
| `func.py`            | Gaussian, Lorentzian, Voigt などの代表的関数群                       |
| `signal1d.py`        | 1次元信号クラス `SignalData`（カーソル・解析機能）                   |
| `map2d.py`           | 2次元マップクラス `Map2D`（平滑化・等高線・トラッキング）            |

---

## 🔧 Usage Example

### 1D Signal Analysis

```python
from pyana.signal1d import SignalData
import numpy as np

x = np.linspace(0, 10, 1000)
y = np.sin(x) + np.random.normal(0, 0.1, x.shape)

sig = SignalData(x, y)
sig.set_cursors(2, 8)
peak_idx = sig.find_peaks_segment()
```

### 2D Map Processing

```python
from pyana.map2d import Map2D
import numpy as np

x = np.linspace(0, 1, 100)
y = np.linspace(0, 2, 200)
z = np.random.rand(len(y), len(x))

mp = Map2D(x, y, z)
mp.plot()
```

---

## 📄 License

MIT License

---

## ✨ Author

Hiroki Tsusaka @1160-hrk  
The University of Tokyo
 → commit

# analysis-utils
