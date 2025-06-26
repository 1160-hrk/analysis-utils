#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分子スペクトロスコピー解析用ユーティリティ関数

Created on Fri May  5 10:02:50 2023
Modified: 2024

@author: hirokitsusaka
"""

import math
from typing import List, Tuple
from dataclasses import dataclass
import numpy as np
import numpy.typing as npt


# 物理定数
@dataclass(frozen=True)
class MolecularConstants:
    """分子定数クラス"""
    vibrational_wavenumbers: npt.NDArray[np.float64]
    rotational_constant_b: float
    rotational_constant_db: float


# CO分子の定数
CO_CONSTANTS = MolecularConstants(
    vibrational_wavenumbers=np.array([
        2349.1, 2324.2, 2299.3, 2274.4, 2249.5, 2224.7,
        2200.0, 2175.3, 2150.7, 2126.1, 2101.7, 2077.4
    ]),
    rotational_constant_b=0.39022,
    rotational_constant_db=0.003076
)

# プロット用色配列
PLOT_COLORS = [
    'brown', 'orangered', 'orange', 'mediumseagreen',
    'mediumblue', 'darkorchid', 'fuchsia', 'crimson',
    'brown', 'orangered', 'orange', 'mediumseagreen'
]


@dataclass
class TransitionResult:
    """遷移結果を格納するデータクラス"""
    wavenumber: float
    vibrational_level: int
    rotational_level: int
    delta_j: int
    difference: float
    color: str = ""


def _calculate_rotational_term(j_quantum_number: int,
                               vibrational_level: int,
                               constants: MolecularConstants) -> float:
    """回転項の値を計算
    
    Args:
        j_quantum_number: 回転量子数
        vibrational_level: 振動準位
        constants: 分子定数
        
    Returns:
        回転項の値 [cm^-1]
    """
    b_eff = (constants.rotational_constant_b -
             constants.rotational_constant_db * vibrational_level)
    return b_eff * j_quantum_number * (j_quantum_number + 1)


def _get_rotational_quantum_number(branch_index: int,
                                   vibrational_level: int) -> int:
    """枝のインデックスから回転量子数を計算
    
    Args:
        branch_index: 枝のインデックス
        vibrational_level: 振動準位
        
    Returns:
        回転量子数
    """
    if vibrational_level % 2 == 0:
        return branch_index * 2
    else:
        return branch_index * 2 + 1


def find_best_transition(frequency: float,
                         max_j: int,
                         max_v: int,
                         constants: MolecularConstants = CO_CONSTANTS
                         ) -> TransitionResult:
    """
    PPスペクトル上のピーク信号に最も近い遷移を特定
    
    Args:
        frequency: 観測された周波数 [cm^-1]
        max_j: 探索する最大回転量子数
        max_v: 探索する最大振動準位
        constants: 分子定数
        
    Returns:
        最適な遷移情報
        
    Raises:
        ValueError: 入力パラメータが不正な場合
        IndexError: 振動準位が範囲外の場合
    """
    if max_j <= 0 or max_v <= 0:
        raise ValueError("max_j and max_v must be positive integers")
    
    if max_v > len(constants.vibrational_wavenumbers):
        raise IndexError(
            f"max_v ({max_v}) exceeds available vibrational levels "
            f"({len(constants.vibrational_wavenumbers)})")
    
    best_result = TransitionResult(
        wavenumber=0.0,
        vibrational_level=0,
        rotational_level=0,
        delta_j=0,
        difference=float('inf')
    )
    
    for v_level in range(max_v):
        for branch_idx in range(max_j):
            j_lower = _get_rotational_quantum_number(branch_idx, v_level)
            
            # P枝 (ΔJ = -1) とR枝 (ΔJ = +1) の両方を計算
            for delta_j in [-1, 1]:
                j_upper = j_lower + delta_j
                if j_upper < 0:  # 負の量子数は物理的に意味がない
                    continue
                
                # 下準位と上準位の回転項を計算
                lower_rotational_term = _calculate_rotational_term(
                    j_lower, v_level, constants)
                upper_rotational_term = _calculate_rotational_term(
                    j_upper, v_level + 1, constants)
                
                # 遷移波数を計算
                transition_wavenumber = (
                    constants.vibrational_wavenumbers[v_level] +
                    upper_rotational_term - lower_rotational_term)
                
                # 観測値との差を計算
                difference = abs(transition_wavenumber - frequency)
                
                if difference < best_result.difference:
                    best_result = TransitionResult(
                        wavenumber=transition_wavenumber,
                        vibrational_level=v_level,
                        rotational_level=j_lower,
                        delta_j=delta_j,
                        difference=difference,
                        color=PLOT_COLORS[v_level % len(PLOT_COLORS)]
                    )
    
    return best_result


def find_transitions_in_range(frequency: float,
                              tolerance: float,
                              max_j: int,
                              max_v: int,
                              constants: MolecularConstants = CO_CONSTANTS
                              ) -> List[TransitionResult]:
    """
    指定した許容範囲内の全ての遷移を検索
    
    Args:
        frequency: 観測された周波数 [cm^-1]
        tolerance: 許容範囲 [cm^-1]
        max_j: 探索する最大回転量子数
        max_v: 探索する最大振動準位
        constants: 分子定数
        
    Returns:
        範囲内の全遷移のリスト
        
    Raises:
        ValueError: 入力パラメータが不正な場合
    """
    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    
    if max_j <= 0 or max_v <= 0:
        raise ValueError("max_j and max_v must be positive integers")
    
    if max_v > len(constants.vibrational_wavenumbers):
        raise IndexError(
            f"max_v ({max_v}) exceeds available vibrational levels")
    
    matching_transitions = []
    
    for v_level in range(max_v):
        for branch_idx in range(max_j):
            j_lower = _get_rotational_quantum_number(branch_idx, v_level)
            
            for delta_j in [-1, 1]:
                j_upper = j_lower + delta_j
                if j_upper < 0:
                    continue
                
                lower_rotational_term = _calculate_rotational_term(
                    j_lower, v_level, constants)
                upper_rotational_term = _calculate_rotational_term(
                    j_upper, v_level + 1, constants)
                
                transition_wavenumber = (
                    constants.vibrational_wavenumbers[v_level] +
                    upper_rotational_term - lower_rotational_term)
                
                difference = abs(transition_wavenumber - frequency)
                
                if difference <= tolerance:
                    matching_transitions.append(TransitionResult(
                        wavenumber=transition_wavenumber,
                        vibrational_level=v_level,
                        rotational_level=j_lower,
                        delta_j=delta_j,
                        difference=difference,
                        color=PLOT_COLORS[v_level % len(PLOT_COLORS)]
                    ))
    
    # 差の小さい順にソート
    matching_transitions.sort(key=lambda x: x.difference)
    return matching_transitions


def calculate_transition_wavenumber(v_lower: int,
                                    j_lower: int,
                                    delta_j: int,
                                    constants: MolecularConstants = CO_CONSTANTS
                                    ) -> float:
    """
    指定した準位間の遷移波数を計算
    
    Args:
        v_lower: 下準位の振動量子数
        j_lower: 下準位の回転量子数
        delta_j: 回転量子数の変化 (+1 for R-branch, -1 for P-branch)
        constants: 分子定数
        
    Returns:
        遷移波数 [cm^-1]
        
    Raises:
        ValueError: 入力パラメータが不正な場合
        IndexError: 振動準位が範囲外の場合
    """
    if v_lower < 0 or j_lower < 0:
        raise ValueError("Quantum numbers must be non-negative")
    
    if delta_j not in [-1, 1]:
        raise ValueError("delta_j must be -1 (P-branch) or +1 (R-branch)")
    
    if v_lower >= len(constants.vibrational_wavenumbers):
        raise IndexError(
            f"v_lower ({v_lower}) exceeds available vibrational levels")
    
    j_upper = j_lower + delta_j
    if j_upper < 0:
        raise ValueError("Upper rotational quantum number cannot be negative")
    
    # 回転項を計算
    lower_rotational_term = _calculate_rotational_term(
        j_lower, v_lower, constants)
    upper_rotational_term = _calculate_rotational_term(
        j_upper, v_lower + 1, constants)
    
    # 遷移波数を計算
    return (constants.vibrational_wavenumbers[v_lower] +
            upper_rotational_term - lower_rotational_term)


def correct_wavelength_ihr320(wavelength: float,
                              angle_correction: float = 0.027 * math.pi / 180,
                              grating_lines_per_meter: float = 300e-6,
                              blaze_angle: float = 21.26 * math.pi / 180
                              ) -> float:
    """
    iHR320分光器での波長補正
    
    Args:
        wavelength: 補正前の波長 [m]
        angle_correction: 角度補正値 [rad]
        grating_lines_per_meter: 回折格子の溝密度 [lines/m]
        blaze_angle: ブレーズ角 [rad]
        
    Returns:
        補正後の波長 [m]
        
    Raises:
        ValueError: 物理的に不正な値の場合
    """
    if wavelength <= 0:
        raise ValueError("Wavelength must be positive")
    
    if grating_lines_per_meter <= 0:
        raise ValueError("Grating lines per meter must be positive")
    
    try:
        # 回折格子の角度を計算
        sin_angle = (wavelength * grating_lines_per_meter /
                     (2 * math.cos(blaze_angle / 2)))
        
        if abs(sin_angle) > 1:
            raise ValueError(f"Invalid grating angle: sin(angle) = {sin_angle}")
        
        grating_angle = math.asin(sin_angle)
        grating_angle += angle_correction
        
        # 補正後の波長を計算
        corrected_wavelength = (math.sin(grating_angle) /
                                grating_lines_per_meter *
                                2 * math.cos(blaze_angle / 2))
        
        return corrected_wavelength
        
    except (ValueError, OverflowError) as e:
        raise ValueError(f"Error in wavelength correction: {e}")


def calculate_wavelength_array(reference_wavelength: float,
                               pixel_index: int,
                               focal_length: float,
                               blaze_angle: float = 21.56 * math.pi / 180,
                               grating_lines_per_meter: float = 300e-6,
                               pixel_size: float = 50e-6,
                               num_pixels: int = 256,
                               reverse_order: bool = True
                               ) -> npt.NDArray[np.float64]:
    """
    CCDピクセルインデックスから波長配列を計算
    
    Args:
        reference_wavelength: 基準波長 [m]
        pixel_index: 基準ピクセルのインデックス
        focal_length: 焦点距離 [m]
        blaze_angle: ブレーズ角 [rad]
        grating_lines_per_meter: 回折格子の溝密度 [lines/m]
        pixel_size: ピクセルサイズ [m]
        num_pixels: ピクセル数
        reverse_order: 配列を反転するかどうか
        
    Returns:
        波長配列 [m]
        
    Raises:
        ValueError: 入力パラメータが不正な場合
    """
    if reference_wavelength <= 0:
        raise ValueError("Reference wavelength must be positive")
    
    if focal_length <= 0:
        raise ValueError("Focal length must be positive")
    
    if pixel_size <= 0:
        raise ValueError("Pixel size must be positive")
    
    if num_pixels <= 0:
        raise ValueError("Number of pixels must be positive")
    
    if grating_lines_per_meter <= 0:
        raise ValueError("Grating lines per meter must be positive")
    
    try:
        # ピクセルインデックス配列を作成
        if reverse_order:
            pixel_indices = pixel_index + np.arange(0, -num_pixels, -1)
        else:
            pixel_indices = np.arange(0, num_pixels) - pixel_index
        
        # 各ピクセルの角度配列を計算
        angle_array = np.arctan(pixel_indices * pixel_size / focal_length)
        
        # 基準角度を計算
        sin_reference_angle = (grating_lines_per_meter * reference_wavelength /
                               (2 * math.cos(blaze_angle / 2)))
        
        if abs(sin_reference_angle) > 1:
            raise ValueError(
                f"Invalid reference angle: sin(angle) = {sin_reference_angle}")
        
        reference_angle = math.asin(sin_reference_angle)
        
        # 波長配列を計算
        wavelength_array = (2 / grating_lines_per_meter *
                            np.cos(blaze_angle / 2 + angle_array / 2) *
                            np.sin(reference_angle + angle_array / 2))
        
        if reverse_order:
            wavelength_array = np.flip(wavelength_array)
        
        return wavelength_array
        
    except (ValueError, OverflowError) as e:
        raise ValueError(f"Error in wavelength array calculation: {e}")


# 後方互換性のための旧関数名（非推奨）
def find_transition(f: float, J_max: int, v_max: int
                    ) -> Tuple[str, int, float, int, int]:
    """非推奨: find_best_transition を使用してください"""
    import warnings
    warnings.warn(
        "find_transition is deprecated. Use find_best_transition instead.",
        DeprecationWarning, stacklevel=2)
    
    result = find_best_transition(f, J_max, v_max)
    return (result.color, result.vibrational_level, result.wavenumber,
            result.rotational_level, result.delta_j)


def find_transition_w_range(f: float, w: float, J_max: int, v_max: int
                            ) -> List[List[float]]:
    """非推奨: find_transitions_in_range を使用してください"""
    import warnings
    warnings.warn(
        "find_transition_w_range is deprecated. " +
        "Use find_transitions_in_range instead.",
        DeprecationWarning, stacklevel=2)
    
    results = find_transitions_in_range(f, w, J_max, v_max)
    return [[r.wavenumber, r.vibrational_level, r.rotational_level,
             r.delta_j, r.difference] for r in results]


def transition_wavenumber(v_l: int, j_l: int, dj: int) -> float:
    """非推奨: calculate_transition_wavenumber を使用してください"""
    import warnings
    warnings.warn(
        "transition_wavenumber is deprecated. " +
        "Use calculate_transition_wavenumber instead.",
        DeprecationWarning, stacklevel=2)
    
    return calculate_transition_wavenumber(v_l, j_l, dj)


def wl_correct_iHR320(wl: float, **kwargs) -> float:
    """非推奨: correct_wavelength_ihr320 を使用してください"""
    import warnings
    warnings.warn(
        "wl_correct_iHR320 is deprecated. " +
        "Use correct_wavelength_ihr320 instead.",
        DeprecationWarning, stacklevel=2)
    
    return correct_wavelength_ihr320(wl, **kwargs)


def wl_array_board(lambda_ref: float, ind_pixel: int, f: float, **kwargs
                   ) -> npt.NDArray[np.float64]:
    """非推奨: calculate_wavelength_array を使用してください"""
    import warnings
    warnings.warn(
        "wl_array_board is deprecated. " +
        "Use calculate_wavelength_array instead.",
        DeprecationWarning, stacklevel=2)
    
    return calculate_wavelength_array(lambda_ref, ind_pixel, f, **kwargs)
