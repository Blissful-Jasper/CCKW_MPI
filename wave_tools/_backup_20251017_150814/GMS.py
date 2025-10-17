#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 15 10:44:51 2025

@author: xpji
"""


from typing import Tuple, Dict, Union, Optional, List, Any
from metpy.calc import mixing_ratio_from_specific_humidity
from metpy.units import units
import time
import os
import glob
from DSE import *
from geocat.comp import delta_pressure
import xarray as xr

def extract_model_name(path: str) -> str:
    return os.path.basename(path).split("_")[2]

def get_unique_file_or_raise(path_pattern: str, varname: str, model: str) -> str:
    files = glob.glob(path_pattern)
    if not files:
        raise FileNotFoundError(f"❌ 未找到 {varname} 文件，模型：{model}")
    if len(files) > 1:
        print(f"⚠️ 警告：找到多个 {varname} 文件，使用第一个：{files[0]}")
    return files[0]

def vertically_integrated_moist_flux_divergence(
    hus: xr.DataArray, ua: xr.DataArray, va: xr.DataArray,
    lat: xr.DataArray, lon: xr.DataArray
) -> xr.DataArray:
    """
    返回单位为 W/m²（能量通量散度，L × ∇·(rV)）
    """
    r = mixing_ratio_from_specific_humidity(hus)
    r = xr.DataArray(r, coords=hus.coords, dims=hus.dims)

    r_u = r * ua
    r_v = r * va

    dx, dy = compute_dx_dy(lat, lon)

    dr_u_dx = xr.DataArray(np.gradient(r_u, axis=-1) / dx, coords=r_u.coords, dims=r_u.dims)
    dr_v_dy = xr.DataArray(np.gradient(r_v, axis=-2) / dy, coords=r_v.coords, dims=r_v.dims)

    div_rV = dr_u_dx + dr_v_dy

    # dp = np.gradient(hus.plev.values * 100)  # Pa
    # weights = (np.abs(dp) / g).reshape(1, -1, 1, 1)

    return L *( div_rV.integrate("plev"))

def calc_horizontal_GMS(
    ta: xr.DataArray, zg: xr.DataArray, plev: np.ndarray,
    lon: xr.DataArray, lat: xr.DataArray,
    ua: xr.DataArray, va: xr.DataArray, hus: xr.DataArray
) -> xr.DataArray:
    """
    输出：水平 GMS，单位：无量纲（单位化的能量通量比）
    """
    dse = calc_dse(ta, zg, plev)  # J/kg

    dx, dy = compute_dx_dy(lat, lon)
    ds_dx = xr.DataArray(np.gradient(dse, axis=-1) / dx, coords=dse.coords, dims=dse.dims)
    ds_dy = xr.DataArray(np.gradient(dse, axis=-2) / dy, coords=dse.coords, dims=dse.dims)

    v_grad_s = ua * ds_dx + va * ds_dy

    iv_dot_grad_s = (v_grad_s).integrate("plev")
    i_r = vertically_integrated_moist_flux_divergence(hus, ua, va, lat, lon)

    return -T_ref * iv_dot_grad_s / i_r

def calc_vertical_GMS(
    ta: xr.DataArray, zg: xr.DataArray, plev: np.ndarray,
    wa: xr.DataArray, hus: xr.DataArray, ua: xr.DataArray, va: xr.DataArray,
    lat: xr.DataArray, lon: xr.DataArray
) -> xr.DataArray:
    """
    输出：垂直 GMS，单位：无量纲
    """
    dse = calc_dse(ta, zg, plev)  # J/kg
    ds_dp = dse.differentiate("plev")  # J/kg/Pa，注意 dse 必须按 plev 定义
    ids_dp = (wa * ds_dp).integrate('plev')  # J/m²/s
    i_r = vertically_integrated_moist_flux_divergence(hus, ua, va, lat, lon)

    return -T_ref * ids_dp / i_r
    
def gross_moist_stability(
    ta_path: str, zg_path: str, ua_path: str, va_path: str, wa_path: str, hus_path: str
) -> tuple[xr.DataArray, xr.DataArray]:
    ta, lon, lat = load_data(ta_path, 'ta')
    zg, _, _     = load_data(zg_path, 'zg')
    ua, _, _     = load_data(ua_path, 'ua')
    va, _, _     = load_data(va_path, 'va')
    wa, _, _     = load_data(wa_path, 'wap')
    hus,_,_      = load_data(hus_path, 'hus')
    plev = ta.plev.values


    h_gms = calc_horizontal_GMS(ta, zg, plev, lon, lat, ua, va, hus)
    v_gms = calc_vertical_GMS(ta, zg, plev, wa, hus, ua, va, lat, lon)

    return h_gms, v_gms
