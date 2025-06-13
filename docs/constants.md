# constants.py

## 🔬 概要

このモジュールは、物理・化学で頻繁に使用される **物理定数** や **単位辞書** をまとめた定義ファイルです。  
`scipy.constants` にある標準定数にエイリアスを与え、加えて FWHM や単位変換等の用途で使えるカテゴリ別辞書を提供します。

---

## 📏 主な定数

| 名前 | 値の概要 |
|------|-----------|
| `PLANCK_CONSTANT` | プランク定数 `h` |
| `PLANCK_CONSTANT_REDUCED` | ディラック定数 `ℏ` |
| `BOLTZMANN_CONSTANT` | ボルツマン定数 `k` |
| `AVOGADRO_CONSTANT` | アボガドロ定数 `N_A` |
| `GAS_CONSTANT` | 気体定数 `R` |
| `SPEED_OF_LIGHT` | 光速 `c` |
| `ELECTRON_MASS` | 電子質量 `m_e` |
| `PROTON_MASS` | 陽子質量 `m_p` |
| `NEUTRON_MASS` | 中性子質量 `m_n` |
| `ELEMENTARY_CHARGE` | 電子の電荷 `e` |
| `PERMITTIVITY_VACUUM` | 真空の誘電率 `ε₀` |
| `PERMEABILITY_VACUUM` | 真空の透磁率 `μ₀` |
| `GRAVITATIONAL_CONSTANT` | 万有引力定数 `G` |
| `BOLTZMANN_CONSTANT_JOULE` | 単位: J/K のボルツマン定数 |
| `BOLTZMANN_CONSTANT_EV` | 単位: eV/K のボルツマン定数 |

---

## 📦 単位辞書一覧

### 🔹 `SI_UNIT`
SI基本単位（m, kg, s, K など）

### 🔸 `CGS_UNIT`
CGS単位系（cm, g, s）

### ⚛ `EV_UNIT`
電子ボルト関連：`"energy": "eV"` など

---

## 💡 専門領域別の単位カテゴリ

- `ENERGY_UNIT`: エネルギー単位（J, eV, cal, kWh など）
- `PRESSURE_UNIT`: 圧力（Pa, atm, torr, bar）
- `TEMPERATURE_UNIT`: 温度（K, °C, °F）
- `ANGLE_UNIT`: 角度（rad, °）
- `SPEED_UNIT`: 速度（m/s, km/h, mph）
- `ACCELERATION_UNIT`: 加速度（m/s² など）
- `FORCE_UNIT`: 力（N, kgf, lbf）
- `POWER_UNIT`: 電力（W, kW）
- `VOLTAGE_UNIT`: 電圧（V, mV）
- `CURRENT_UNIT`: 電流（A, mA）
- `RESISTANCE_UNIT`: 抵抗（Ω, kΩ）
- `CAPACITANCE_UNIT`: 静電容量（F, μF）
- `INDUCTANCE_UNIT`: 誘導（H, mH）
- `MAGNETIC_FIELD_UNIT`: 磁場（T, G）
- `ELECTRIC_FIELD_UNIT`: 電場（V/m）
- `ELECTRIC_CHARGE_UNIT`: 電気量（C, mAh）
- `ELECTRIC_CAPACITY_UNIT`: 電気容量（F など）
- `ELECTRIC_INDUCTANCE_UNIT`: 電気インダクタンス
- `ELECTRIC_MAGNETIC_UNIT`: 磁束（Wb など）

---

## ✅ 使用例

```python
from pyana import constants as C

energy_ev = 0.1  # eV
energy_joule = energy_ev / C.ELEMENTARY_CHARGE
print(f"Energy in J = {energy_joule:.3e} J")
```

---

## 📚 関連

- `scipy.constants`: 基礎定数の元定義（単位：SI）
- `pyana.constants`: 上記にわかりやすい名前と単位辞書を付加