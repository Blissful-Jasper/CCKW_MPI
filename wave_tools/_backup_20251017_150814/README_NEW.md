# Wave Tools - 热带大气波动分析工具包

> **简洁、清晰、实用的Python工具包，专注于热带大气波动分析**

**作者**: Jianpu | **邮箱**: xianpuji@hhu.edu.cn | **机构**: Hohai University

---

## 📦 安装

```bash
# 安装依赖
pip install -r requirements.txt

# 在Python中导入
from wave_tools import *
```

---

## 🎯 核心功能

本工具包提供**热带大气波动分析**的完整解决方案：

1. **Matsuno理论模态计算** - 计算赤道浅水方程的理论色散关系
2. **Wheeler-Kiladis频谱分析** - 波数-频率功率谱诊断
3. **波动滤波与提取** - 时空滤波提取Kelvin、ER、MRG等波动
4. **波动相位分析** - 峰值检测和相位合成
5. **诊断工具** - 湿稳定度(GMS)、Radon变换等
6. **专业绘图** - WK频谱图、地图、泰勒图等

---

## 📁 模块结构

```
wave_tools/
├── matsuno.py       # Matsuno理论模态
├── spectral.py      # 频谱分析
├── filters.py       # 波动滤波
├── phase.py         # 相位分析
├── diagnostics.py   # 诊断工具
├── plotting.py      # 所有绘图功能
└── utils.py         # 工具函数
```

---

## 🚀 快速开始

### 示例1: Wheeler-Kiladis频谱分析

```python
from wave_tools.spectral import calculate_wk_spectrum
import xarray as xr

# 加载数据
data = xr.open_dataset('olr.nc')['olr']

# 计算频谱（一行代码）
power_sym, power_asym, background = calculate_wk_spectrum(
    data, 
    window_days=96, 
    skip_days=30,
    output_path='spectrum.nc'
)
```

### 示例2: 提取Kelvin波

```python
from wave_tools.filters import WaveFilter

wf = WaveFilter()
kelvin = wf.extract_wave_signal(data, wave_name='kelvin', use_parallel=True)
kelvin.to_netcdf('kelvin_filtered.nc')
```

### 示例3: 绘制频谱图

```python
from wave_tools.plotting import plot_wk_spectrum

plot_wk_spectrum(
    power_sym, power_asym, background,
    wavenumber, frequency,
    add_matsuno_lines=True,
    save_path='wk_spectrum.png'
)
```

---

## 📚 完整函数索引

### 1. matsuno.py - Matsuno理论模态

#### 主要函数

| 函数 | 功能 | 参数 |
|------|------|------|
| `kelvin_mode(he, latitude, max_wn, n_wn)` | Kelvin波色散关系 | he: 等效深度(m)<br/>max_wn: 最大波数<br/>n_wn: 波数点数 |
| `er_n(he, n, latitude, max_wn, n_wn)` | 赤道Rossby波色散关系 | n: 经向模态数 |
| `mrg_mode(he, latitude, max_wn, n_wn)` | 混合Rossby重力波 | 同上 |
| `eig_n(he, n, latitude, max_wn, n_wn)` | 东传惯性重力波 | 同上 |
| `wig_n(he, n, latitude, max_wn, n_wn)` | 西传惯性重力波 | 同上 |
| `matsuno_modes_wk(he, n, max_wn)` | **批量计算所有模态** | he: 等效深度列表<br/>n: 经向模态数列表 |

**返回**: DataFrame，index为波数，columns为各模态的频率

**应用**: 生成理论框架，用于对比观测/模式的波动功率谱

---

### 2. spectral.py - Wheeler-Kiladis频谱分析

#### 主要类

**`WKSpectralAnalysis`** - 频谱分析主类

| 方法 | 功能 | 参数说明 |
|------|------|----------|
| `load_data(data, data_path, variable, lat_range, time_range)` | 加载数据 | data: xr.DataArray（可选）<br/>data_path: NetCDF路径（可选）<br/>lat_range: 纬度范围，如(-15,15) |
| `preprocess()` | 预处理 | 去趋势、去年循环、对称/反对称分解 |
| `compute_spectrum()` | 计算功率谱 | 波数-频率2D FFT |
| `smooth_background(wave_limit)` | 平滑背景谱 | wave_limit: 平滑的波数限制 |
| `save(output_path)` | 保存到NetCDF | output_path: 输出文件路径 |

**配置类**: `SpectralConfig`

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `WINDOW_SIZE_DAYS` | 96 | 时间窗口大小（天） |
| `WINDOW_SKIP_DAYS` | 30 | 窗口跳跃间隔（天） |
| `SAMPLES_PER_DAY` | 1 | 每天采样次数 |

#### 便捷函数

| 函数 | 功能 | 返回 |
|------|------|------|
| `calculate_wk_spectrum(data, window_days, skip_days, output_path)` | **一步完成频谱计算** | power_symmetric, power_antisymmetric, background |

---

### 3. filters.py - 波动滤波

#### 主要类

**`WaveFilter`** - 波动滤波器

| 方法 | 功能 | 参数说明 |
|------|------|----------|
| `extract_wave_signal(ds, wave_name, obs_per_day, use_parallel, n_jobs)` | **提取波动信号** | ds: 输入数据<br/>wave_name: 波动类型<br/>use_parallel: 是否并行 |
| `add_wave_param(wave_name, freq_range, wnum_range, equiv_depth)` | 添加自定义波动 | freq_range: 周期范围(天)<br/>wnum_range: 波数范围<br/>equiv_depth: 等效深度(m) |
| `get_available_waves()` | 获取可用波动列表 | 返回: ['kelvin', 'er', 'mrg', ...] |
| `get_wave_params(wave_name)` | 查看波动参数 | 返回: dict |

#### 预定义波动类型

| 波动 | 周期(天) | 波数 | 等效深度(m) | 特征 |
|------|---------|------|------------|------|
| `'kelvin'` | 3-20 | 2-14 | 8-90 | 东传，赤道对称 |
| `'er'` | 9-72 | -10~-1 | 8-90 | 西传，低频 |
| `'mrg'` | 3-10 | -10~-1 | 8-90 | 西传，高频 |
| `'ig'` | 1-14 | 1-5 | 8-90 | 东传惯性重力波 |
| `'mjo'` | 20-100 | 1-5 | NaN | Madden-Julian振荡 |
| `'td'` | 2.5-5 | -20~-6 | NaN | 热带低压型扰动 |

---

### 4. phase.py - 波动相位分析

#### 主要函数

| 函数 | 功能 | 参数 | 返回 |
|------|------|------|------|
| `optimize_peak_detection(V, kelvin_ref, V_std, Nstd, use_parallel)` | **检测显著峰值** | V: 数据数组<br/>V_std: 标准差<br/>Nstd: 阈值倍数 | V_peak: 峰值影响范围<br/>local_extrema: 极值标识 |
| `remove_clm(data, fs, cutoff, order)` | 去除季节循环 | cutoff: 截止频率<br/>order: Butterworth滤波器阶数 | 滤波后的异常场 |
| `find_local_extrema(V)` | 查找局部极值 | V: [time, lon]数组 | 1=极小值，-1=极大值 |
| `butter_lowpass_filter(data, cutoff_freq, fs, order)` | Butterworth低通滤波 | fs: 采样频率 | 滤波后的序列 |

**应用**: 相位合成分析、传播速度计算、生命周期统计

---

### 5. diagnostics.py - 诊断工具

#### 主要函数

| 函数 | 功能 | 参数 | 返回 |
|------|------|------|------|
| `calc_horizontal_GMS(ta, zg, plev, lon, lat, ua, va, hus)` | **水平湿稳定度** | ta: 温度<br/>zg: 位势高度<br/>ua,va: 风场<br/>hus: 比湿 | GMS (无量纲) |
| `calc_vertical_GMS(ta, zg, plev, wa, hus, ua, va, lat, lon)` | **垂直湿稳定度** | wa: 垂直速度 | GMS (无量纲) |
| `gross_moist_stability(ta_path, zg_path, ua_path, va_path, wa_path, hus_path)` | **完整GMS计算** | 各变量文件路径 | h_gms, v_gms |

**物理意义**: 
- GMS > 0: 负反馈，抑制对流
- GMS < 0: 正反馈，促进对流组织化

---

### 6. plotting.py - 绘图功能

#### CCKW包络图

| 函数 | 功能 | 参数 |
|------|------|------|
| `get_cckw_envelope_curve(he, fmax)` | 计算包络曲线坐标 | he: 等效深度列表<br/>fmax: 最大频率列表 |
| `plot_cckw_envelope(he, fmax, save_path, dpi)` | 绘制CCKW包络图 | save_path: 保存路径 |

#### Wheeler-Kiladis频谱图

| 函数 | 功能 | 关键参数 |
|------|------|----------|
| `plot_wk_spectrum(power_symmetric, power_antisymmetric, background, wavenumber, frequency, ...)` | **绘制WK频谱** | max_wn: 最大波数<br/>max_freq: 最大频率<br/>add_matsuno_lines: 是否添加理论曲线<br/>he: 等效深度列表<br/>cpd_lines: 标注周期线 |

#### 地图绘图

| 函数 | 功能 | 参数 |
|------|------|------|
| `setup_map_axes(ax, title, box)` | 设置地图坐标轴 | box: [lon_min, lon_max, lat_min, lat_max] |
| `plot_spatial_field(data, ax, cmap, title, box, levels)` | 绘制空间场 | data: xr.DataArray<br/>cmap: 色标 |

#### 泰勒图

**类**: `TaylorDiagram(refstd, fig, rect, srange, extend)`

| 方法 | 功能 | 参数 |
|------|------|------|
| `add_sample(stddev, corrcoef, *args, **kwargs)` | 添加模式样本点 | stddev: 标准差<br/>corrcoef: 相关系数 |
| `add_contours(levels, **kwargs)` | 添加RMS等值线 | levels: 等值线水平数 |

#### 通用工具

| 函数 | 功能 | 参数 |
|------|------|------|
| `save_figure(fig, filename, folder, fmt, dpi)` | 保存图形 | fmt: 'pdf'/'png'/'jpg'<br/>dpi: 分辨率 |

---

### 7. utils.py - 工具函数

#### 数据处理

| 函数 | 功能 | 参数 | 返回 |
|------|------|------|------|
| `load_data(path, var, lat_range)` | 加载NetCDF数据 | path: 文件路径<br/>var: 变量名<br/>lat_range: 纬度范围 | data, lon, lat |
| `filter_series(series, min_wn, max_wn)` | 波数范围过滤 | min_wn, max_wn: 波数边界 | 过滤后的序列 |
| `filter_paths_by_models(paths, model_names, loc, sep)` | 按模式名过滤文件 | paths: 文件路径列表<br/>model_names: 模式名列表 | 匹配的路径列表 |
| `extract_model_name(path, loc, sep)` | 从路径提取模式名 | loc: 分隔后的位置<br/>sep: 分隔符 | 模式名 |

#### HEALPix网格转换

| 函数 | 功能 | 参数 | 返回 |
|------|------|------|------|
| `dataarray_healpix_to_equatorial_latlon(healpix_dataarray, nside, nest, minmax_lat)` | **HEALPix→等经纬度** | nside: HEALPix参数<br/>nest: 是否nested<br/>minmax_lat: 纬度范围 | 等经纬度DataArray |
| `get_region_healpix_(zoom, extent, nest)` | 获取区域网格索引 | extent: [lon_min, lon_max, lat_min, lat_max] | 网格索引数组 |

**优化**: 使用Numba加速，速度提升10倍

#### Radon变换

| 函数 | 功能 | 参数 | 返回 |
|------|------|------|------|
| `calc_radon_angle(field, theta_range)` | 计算传播角度 | field: 2D数组 | theta, intensity, theta_max |
| `calc_c_from_theta(theta_deg, dx_deg, dt_sec, lat)` | 计算相速度 | theta_deg: 角度(度)<br/>dx_deg: 经度间隔<br/>dt_sec: 时间间隔 | 速度(m/s) |
| `plot_radon_energy_distribution(theta, energy, title, color, ax)` | 绘制能量分布 | energy: Radon能量 | theta_max, theta_ci |

#### 其他工具

| 函数 | 功能 | 参数 |
|------|------|------|
| `create_cmap_from_string(color_string)` | 创建自定义色标 | color_string: 颜色列表字符串 |
| `set_matplotlib_font(font_dir, arial_font_path)` | 设置Matplotlib字体 | font_dir: 字体目录 |

---

## 📖 使用示例

### 完整工作流程

```python
import xarray as xr
from wave_tools.spectral import WKSpectralAnalysis, SpectralConfig
from wave_tools.filters import WaveFilter
from wave_tools.phase import optimize_peak_detection, remove_clm
from wave_tools.plotting import plot_wk_spectrum
from wave_tools.diagnostics import gross_moist_stability

# ===== 1. 频谱诊断 =====
data = xr.open_dataset('olr.nc')['olr'].sel(lat=slice(-15, 15))

config = SpectralConfig()
config.WINDOW_SIZE_DAYS = 96

analysis = WKSpectralAnalysis(config)
analysis.load_data(data=data)
analysis.preprocess()
analysis.compute_spectrum()
analysis.smooth_background()
analysis.save('spectrum.nc')

# 绘制频谱
plot_wk_spectrum(
    analysis.power_symmetric,
    analysis.power_antisymmetric,
    analysis.background,
    analysis.wavenumber.values,
    analysis.frequency.values,
    add_matsuno_lines=True,
    save_path='wk_spectrum.png'
)

# ===== 2. 波动提取 =====
wf = WaveFilter()
kelvin = wf.extract_wave_signal(data, wave_name='kelvin', use_parallel=True)
er = wf.extract_wave_signal(data, wave_name='er', use_parallel=True)

# 保存结果
kelvin.to_netcdf('kelvin_filtered.nc')
er.to_netcdf('er_filtered.nc')

# ===== 3. 相位分析 =====
kelvin_lowpass = remove_clm(kelvin, cutoff=1/10, order=4)
kelvin_zonal = kelvin_lowpass.mean('lat')

V_peak, local_extrema = optimize_peak_detection(
    V=kelvin_zonal.values,
    kelvin_ref=kelvin_zonal,
    V_std=kelvin_zonal.std().values,
    Nstd=1.0,
    use_parallel=True
)

# ===== 4. 诊断分析 =====
h_gms, v_gms = gross_moist_stability(
    'ta.nc', 'zg.nc', 'ua.nc', 'va.nc', 'wa.nc', 'hus.nc'
)
```

---

## 🔧 依赖库

### 核心依赖
```
numpy >= 1.19.0
xarray >= 0.16.0
scipy >= 1.5.0
matplotlib >= 3.3.0
pandas >= 1.1.0
```

### 可选依赖
```
cmaps          # 气象色标
cartopy        # 地图投影
healpy         # HEALPix
numba          # 加速计算
joblib         # 并行计算
metpy          # 气象工具
```

---

## 📝 引用

如果本工具包对您的研究有帮助，请引用：

```
Jianpu. (2025). Wave Tools: A Python package for tropical atmospheric wave analysis. 
Hohai University. Email: xianpuji@hhu.edu.cn
```

---

## 📧 联系方式

**问题与建议**: xianpuji@hhu.edu.cn

**参考文献**:
- Wheeler & Kiladis (1999). *J. Atmos. Sci.*, 56(3), 374-399.
- Matsuno (1966). *J. Meteor. Soc. Japan*, 44(1), 25-43.
- Kiladis et al. (2009). *Rev. Geophys.*, 47(2).

---

**版本**: 1.0.0 | **最后更新**: 2025-10-17
