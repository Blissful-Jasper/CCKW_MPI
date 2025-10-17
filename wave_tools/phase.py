#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 15 20:55:03 2025

@author: xpji
"""
import xarray as xr
import numpy as np
from joblib import Parallel, delayed
from typing import Optional, List, Tuple
from scipy.signal import butter, filtfilt
import time


def butter_lowpass_filter(data: np.ndarray, cutoff_freq: float, fs: float, order: int = 4) -> np.ndarray:
    """
    应用 Butterworth 低通滤波器，滤除高于截止频率的信号成分。
    
    Parameters
    ----------
    data : np.ndarray
        一维时间序列数据。
    cutoff_freq : float
        截止频率 (Hz)。例如 1/10 表示滤除周期小于 10 天的信号。
    fs : float
        采样频率 (Hz)。例如逐日数据为 1。
    order : int, optional
        滤波器阶数，默认值为 4。
    
    Returns
    -------
    np.ndarray
        滤波后的时间序列。
    """
    nyq = 0.5 * fs  # 奈奎斯特频率
    normalized_cut = cutoff_freq / nyq
    b, a = butter(order, normalized_cut, btype='lowpass')
    return filtfilt(b, a, data)


def remove_10d_from_daily_data(data: np.ndarray, fs: float = 1.0, cutoff: float = 1/10, order: int = 4) -> np.ndarray:
    """
    对三维数组 (time, lat, lon) 应用低通滤波，滤除周期小于10天的信号。
    
    Parameters
    ----------
    data : np.ndarray
        三维数组 (time, lat, lon)。
    fs : float
        采样频率 (默认 1/day)。
    cutoff : float
        截止频率，默认 1/10。
    order : int
        Butterworth滤波器阶数。
    
    Returns
    -------
    np.ndarray
        滤波后的三维数组。
    """
    nt, nlat, nlon = data.shape
    filtered = np.empty_like(data)
    for i in range(nlat):
        for j in range(nlon):
            filtered[:, i, j] = butter_lowpass_filter(data[:, i, j], cutoff_freq=cutoff, fs=fs, order=order)
    return filtered


def remove_clm(data: xr.DataArray, fs: float = 1.0, cutoff: float = 1/10, order: int = 4) -> xr.DataArray:
    """
    去除季节循环气候态并对异常进行低通滤波，提取长于10天的信号分量。

    Parameters
    ----------
    data : xr.DataArray
        输入的三维逐日数据，需包含 'time', 'lat', 'lon' 三个维度。
    fs : float
        采样频率，默认每日为 1。
    cutoff : float
        截止频率，默认周期为 10 天。
    order : int
        Butterworth 滤波器阶数。

    Returns
    -------
    xr.DataArray
        滤波后的异常数据（时间序列滤除了 10 天以下信号）。
    """
    if not isinstance(data, xr.DataArray):
        raise TypeError("输入数据必须是 xarray.DataArray 类型")
    if 'time' not in data.dims:
        raise ValueError("输入数据必须包含 'time' 维度")

    # 去除气候态（季节循环）
    climatology = data.mean('time')
    anomalies = data - climatology

    # 低通滤波（保留10天以上周期）
    # filtered_data = remove_10d_from_daily_data(anomalies.values, fs=fs, cutoff=cutoff, order=order)

    return xr.DataArray(
        data=anomalies,
        dims=data.dims,
        coords=data.coords,
        attrs=data.attrs,
        name='lowpass_filtered_anomaly'
    )


def find_local_extrema(V: np.ndarray) -> np.ndarray:
    """
    查找二维时序数据中的局部最大值和最小值点。
    
    参数:
    ----------
    V : np.ndarray
        输入二维数组，形状为 [time, lon]

    返回:
    ----------
    np.ndarray
        与输入同形状的数组，其中：
        - 1 表示局部最小值
        - -1 表示局部最大值
        - np.nan 表示其他位置
    """
    nt, nlon = V.shape
    local_min_max_id = np.full((nt, nlon), np.nan, dtype=np.float32)

    prev = V[:-2, :]
    curr = V[1:-1, :]
    next_ = V[2:, :]

    is_local_min = (curr <= prev) & (curr <= next_)
    is_local_max = (curr >= prev) & (curr >= next_)

    local_min_max_id[1:-1, :] = np.where(is_local_min, 1, local_min_max_id[1:-1, :])
    local_min_max_id[1:-1, :] = np.where(is_local_max, -1, local_min_max_id[1:-1, :])

    return local_min_max_id


def find_peak_influence_range(
    peak_idx: int,
    peak_value: float,
    zero_idx: np.ndarray,
    peak_indices: np.ndarray,
    V: np.ndarray,
    V_std: float,
    Nstd: float
) -> Tuple[int, int]:
    """
    确定某个峰值在时间序列中的影响范围。

    参数:
    ----------
    peak_idx : int
        峰值的时间索引。
    peak_value : float
        峰值对应的值。
    zero_idx : np.ndarray
        零点索引数组。
    peak_indices : np.ndarray
        所有峰值索引。
    V : np.ndarray
        当前经度对应的原始数据 (1D 时间序列)。
    V_std : float
        标准差阈值。
    Nstd : float
        判断显著峰值的标准差倍数。

    返回:
    ----------
    Tuple[int, int]
        峰值影响范围的左右时间索引边界 (inclusive)。
    """
    if np.abs(peak_value) < V_std * Nstd:
        return peak_idx, peak_idx

    # 与零点距离（正数 = 零点在右边，负数 = 零点在左边）
    dpeak_zero = zero_idx - peak_idx

    # 右边界
    pos_dist = np.where(dpeak_zero >= 0, dpeak_zero, np.inf)
    if np.all(np.isinf(pos_dist)):
        id_r = peak_idx
    else:
        id_r = zero_idx[np.argmin(pos_dist)]
        next_peaks = peak_indices[peak_indices > peak_idx]
        if next_peaks.size > 0:
            id_r = min(id_r, next_peaks.min() - 1)

    # 左边界
    neg_dist = np.where(dpeak_zero < 0, dpeak_zero, -np.inf)
    if np.all(neg_dist == -np.inf):
        id_l = peak_idx
    else:
        id_l = zero_idx[np.argmax(neg_dist)] + 1
        prev_peaks = peak_indices[peak_indices < peak_idx]
        if prev_peaks.size > 0:
            id_l = max(id_l, prev_peaks.max() + 1)

    return int(id_l), int(id_r)


   
def process_single_longitude(
    ilon: int,
    V: np.ndarray,
    local_min_max_id: np.ndarray,
    V_std: float,
    Nstd: float
) -> np.ndarray:
    """
    处理某个经度的时间序列，提取显著峰值并标记其影响范围。

    参数:
    ----------
    ilon : int
        经度索引。
    V : np.ndarray
        原始二维数据数组 [time, lon]。
    local_min_max_id : np.ndarray
        局部极值标记数组 [time, lon]，值为 1 (min), -1 (max), or NaN。
    V_std : float
        当前数据的标准差。
    Nstd : float
        显著性阈值，标准差倍数。

    返回:
    ----------
    np.ndarray
        一维数组 [time]，表示该经度上每个时刻对应的显著峰值或 NaN。
    """
    nt = V.shape[0]
    V_peak_lon = np.full(nt, np.nan, dtype=np.float32)

    # 找零交叉点（符号变化），时间维度上
    zero_idx = np.where(V[:-1, ilon] * V[1:, ilon] <= 0)[0]

    # 查找局部极值点索引
    peak_idx = np.where(np.abs(local_min_max_id[:, ilon]) == 1)[0]
    if len(peak_idx) == 0:
        return V_peak_lon

    # 遍历所有峰值点，判断是否显著，确定左右影响范围，并赋值
    for idx in peak_idx:
        peak_value = V[idx, ilon]
        if np.abs(peak_value) < V_std * Nstd:
            continue

        id_l, id_r = find_peak_influence_range(
            peak_idx=idx,
            peak_value=peak_value,
            zero_idx=zero_idx,
            peak_indices=peak_idx,
            V=V[:, ilon],
            V_std=V_std,
            Nstd=Nstd
        )

        V_peak_lon[id_l:id_r + 1] = peak_value

    return V_peak_lon    


def optimize_peak_detection(
    V: np.ndarray,
    kelvin_ref: xr.DataArray,
    V_std: float,
    Nstd: float = 1.0,
    use_parallel: bool = True,
    n_jobs: int = -1
) -> Tuple[xr.DataArray, xr.DataArray]:
    """
    主函数：执行局部极值查找与显著峰值影响范围分配

    参数:
        V (np.ndarray): 输入二维数组，形状为 [时间, 经度]，单位任意
        kelvin_ref (xr.DataArray): xarray 结构的参考数据，提供坐标与维度信息
        V_std (float): 输入数据的标准差（单位同V）
        Nstd (float, 可选): 判定显著性峰值的标准差倍数阈值，默认1
        use_parallel (bool, 可选): 是否使用并行处理（基于joblib），默认True
        n_jobs (int, 可选): 并行核数，-1为全部核心

    返回:
        V_peak_da (xr.DataArray): 峰值影响范围值分配结果，单位同V
        local_min_max_da (xr.DataArray): 局部极值标识结果，1为最小值，-1为最大值，非极值为NaN
    """
    nt, nlon = V.shape

    # Step 1: 查找局部极值
    t0 = time.time()
    local_min_max_id = find_local_extrema(V)
    print(f"[√] 局部极值检测耗时: {time.time() - t0:.2f} 秒")

    # 包装为 DataArray
    local_min_max_da = xr.DataArray(
        data=local_min_max_id,
        dims=kelvin_ref.dims,
        coords=kelvin_ref.coords,
        name='local_extrema'
    )

    # Step 2: 每个经度上进行峰值影响分配
    t1 = time.time()

    if use_parallel:
        # 并行处理每个经度
        V_peak_cols = Parallel(n_jobs=n_jobs)(
            delayed(process_single_longitude)(ilon, V, local_min_max_id, V_std, Nstd)
            for ilon in range(nlon)
        )
        V_peak = np.column_stack(V_peak_cols)
    else:
        # 顺序处理
        V_peak = np.full((nt, nlon), np.nan, dtype=np.float32)
        for ilon in range(nlon):
            V_peak[:, ilon] = process_single_longitude(ilon, V, local_min_max_id, V_std, Nstd)

    print(f"[√] 峰值影响范围赋值耗时: {time.time() - t1:.2f} 秒")

    # 包装为 DataArray
    V_peak_da = xr.DataArray(
        data=V_peak,
        dims=kelvin_ref.dims,
        coords=kelvin_ref.coords,
        name='peak_influence'
    )

    return V_peak_da, local_min_max_da






























