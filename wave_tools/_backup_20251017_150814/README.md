# Wave Tools - 热带大气波动分析工具包

[![Python Version](https://img.shields.io/badge/python-3.7%2B-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

> **作者**: Jianpu  
> **联系**: xianpuji@hhu.edu.cn  
> **机构**: Hohai University

---

## 📖 目录

- [简介](#简介)
- [安装要求](#安装要求)
- [模块结构](#模块结构)
- [核心功能](#核心功能)
- [使用示例](#使用示例)
- [API 文档](#api-文档)
- [更新日志](#更新日志)

---

## 🌊 简介

`wave_tools` 是一个专门用于热带大气波动分析的 Python 工具包，提供了从数据预处理、波动提取、频谱分析到结果可视化的完整工作流程。该工具包特别适用于：

- 🌀 **热带波动分析**: Kelvin波、赤道Rossby波(ER)、混合Rossby重力波(MRG)、惯性重力波(IG)
- 📊 **频谱分析**: Wheeler-Kiladis功率谱分析
- 🎯 **波动提取**: 基于时空滤波的波动信号提取
- 📈 **可视化**: 专业的科学绘图和泰勒图
- 🔬 **诊断分析**: 湿稳定度(GMS)计算、相位分析

---

## 💻 安装要求

### 必需依赖
```bash
numpy >= 1.19.0
xarray >= 0.16.0
scipy >= 1.5.0
matplotlib >= 3.3.0
pandas >= 1.1.0
```

### 可选依赖
```bash
cmaps          # 气象专用色标
cartopy        # 地图投影
healpy         # HEALPix数据处理
numba          # 加速计算
joblib         # 并行计算
metpy          # 气象计算工具
```

### 安装方法
```bash
# 克隆或下载本仓库
git clone <repository_url>
cd wave_tools

# 安装依赖
pip install -r requirements.txt
```

---

## 📁 模块结构

```
wave_tools/
│
├── __init__.py                    # 包初始化文件
│
├── 核心模块
│   ├── core.py                    # CCKW包络绘图工具
│   ├── ma.py                      # Matsuno模态计算（色散关系）
│   └── utils.py                   # 通用工具函数集合
│
├── 频谱分析模块
│   ├── Spectral_Analysis.py       # Wheeler-Kiladis频谱分析（规则网格）
│   ├── spectral_analysis_icon.py  # 频谱分析（非结构网格/ICON模式）
│   ├── SpectrumPlotter.py         # 频谱图绘制工具
│   └── functions.py               # 频谱相关辅助函数
│
├── 波动分析模块
│   ├── Wave_Filter.py             # 波动滤波器（WK99方法）
│   └── Wave_Phase.py              # 波动相位和峰值检测
│
├── 可视化模块
│   ├── plot.py                    # 时空场绘图、趋势分析绘图
│   └── TaylorDiagram.py           # 泰勒图（模式评估）
│
└── 诊断分析模块
    └── GMS.py                     # 总体湿稳定度(Gross Moist Stability)计算
```

---

## 🔧 核心功能

### 1. **Matsuno模态与色散关系** (`ma.py`)

计算热带大气波动的理论色散关系，包括：
- Kelvin波
- 赤道Rossby波 (ER)
- 混合Rossby重力波 (MRG)
- 惯性重力波 (IG, WIG)

**主要函数**:
- `kelvin_mode(he, latitude, max_wn, n_wn)` - 计算Kelvin波色散曲线
- `er_n(he, n, latitude, max_wn, n_wn)` - 计算ER波色散曲线
- `matsuno_modes_wk(he, n, max_wn)` - 批量计算多个等效深度的Matsuno模态

**使用场景**: 生成理论波动能谱框架，用于对比观测/模式的波动功率谱。

---

### 2. **Wheeler-Kiladis频谱分析** (`Spectral_Analysis.py`)

实现经典的Wheeler-Kiladis (1999) 波数-频率谱分析方法。

**核心类**: `SpectralAnalysis`

**主要方法**:
```python
# 数据加载
.load_data(data_path, variable, lat_range, time_range)

# 预处理（去趋势、去年循环、对称/反对称分解）
.preprocess_data()

# 计算功率谱
.compute_spectra()

# 平滑背景谱
.smooth_background()

# 绘制谱图
.plot_spectra(filename, output_path)

# 保存结果
.save_spectra(output_path)
```

**典型工作流程**:
1. 数据加载 → 2. 去气候态 → 3. 对称/反对称分解 → 4. 2D FFT → 5. 功率谱计算 → 6. 背景谱平滑 → 7. 归一化绘图

**输出**: 对称分量功率谱、反对称分量功率谱、背景谱

---

### 3. **波动滤波与提取** (`Wave_Filter.py`)

基于WK99时空滤波方法提取特定波动信号。

**核心类**: `WaveFilter`

**预定义波动类型**:
| 波动类型 | 周期范围(天) | 波数范围 | 等效深度(m) |
|---------|------------|---------|-----------|
| Kelvin  | 3-20       | 2-14    | 8-90      |
| ER      | 9-72       | -10 - -1| 8-90      |
| MRG     | 3-10       | -10 - -1| 8-90      |
| IG      | 1-14       | 1-5     | 8-90      |
| MJO     | 20-100     | 1-5     | NaN       |
| TD      | 2.5-5      | -20 - -6| NaN       |

**主要方法**:
```python
wf = WaveFilter()

# 提取特定波动
kelvin_signal = wf.extract_wave_signal(
    ds=data,
    wave_name='kelvin',
    obs_per_day=1,
    use_parallel=True,
    n_jobs=-1
)

# 添加自定义波动参数
wf.add_wave_param(
    wave_name='custom_wave',
    freq_range=(5, 15),
    wnum_range=(3, 10),
    equiv_depth=(12, 50)
)
```

**关键特性**:
- 支持并行计算（joblib）
- 自动去除年循环
- 色散关系约束滤波
- 与NCL结果对比验证

---

### 4. **波动相位与峰值检测** (`Wave_Phase.py`)

识别波动的局部极值和相位信息。

**主要函数**:
```python
# 优化的峰值检测（支持并行）
V_peak, local_extrema = optimize_peak_detection(
    V=data.values,
    kelvin_ref=data,
    V_std=data.std().values,
    Nstd=1.0,
    use_parallel=True
)

# 去除季节循环并低通滤波
filtered_anomaly = remove_clm(data, cutoff=1/10, order=4)
```

**应用**: 
- 波动相位合成分析
- 波动传播速度计算
- 波动生命周期统计

---

### 5. **频谱可视化** (`SpectrumPlotter.py`)

专业的Wheeler-Kiladis频谱图绘制工具。

**核心类**: `SpectrumPlotter`

**特色功能**:
- 自动叠加Matsuno模态理论曲线
- 标注Kelvin、ER、IG等波动区域
- 自定义色标和等值线范围
- 支持模式间对比（CMIP6等）

**示例**:
```python
plotter = SpectrumPlotter(cpd_lines=[3, 6, 30])

plotter.plot_cmip_future(
    cmip_f=[symmetric_spectrum],
    wn=wavenumbers,
    fq=frequencies,
    max_wn_plot=15,
    max_freq_plot=0.5,
    matsuno_lines=True,
    he=[8, 25, 90],
    ofil='spectrum_plot'
)
```

---

### 6. **通用工具函数** (`utils.py`)

提供基础设施支持的工具函数集合。

**分类功能**:

#### 数据处理
- `load_data(path, var, lat_range)` - 标准化数据加载
- `filter_series(series, min_wn, max_wn)` - 波数范围过滤
- `filter_paths_by_models(paths, model_names)` - 按模式名称过滤文件路径

#### HEALPix网格转换（优化版）
- `dataarray_healpix_to_equatorial_latlon()` - HEALPix → 等经纬度网格
- 使用Numba加速的周期性插值
- 支持向量化批处理

#### Radon变换分析
- `calc_radon_angle(field)` - 计算波动传播角度
- `calc_c_from_theta(theta, dx, dt, lat)` - 计算相速度
- `plot_radon_energy_distribution()` - 绘制Radon能量分布

#### 图形工具
- `save_figure(fig, filename, folder, fmt, dpi)` - 标准化图像保存
- `create_cmap_from_string(color_string)` - 自定义色标创建
- `get_curve(he, fmax)` - CCKW包络曲线坐标

---

### 7. **地图绘图工具** (`plot.py`)

基于Cartopy的地理场绘图和趋势分析。

**主要函数**:
```python
# 设置地图子图
make_space_fig(ax, title='Map Title', box=[0, 360, -30, 30])

# 绘制空间数据
plot_space_data(
    data=pr_data,
    ax=ax,
    cmap='RdBu_r',
    fmt='mm/day',
    title='Precipitation',
    levels=np.linspace(0, 10, 21)
)

# 多波动年际趋势对比
plot_multiple_wave_trends(
    wave_filters=[kelvin, er, mrg, td],
    wave_names=['Kelvin', 'ER', 'MRG', 'TD'],
    regions=['Global', 'Indian-Pacific', 'Indian'],
    save_path='wave_trends.png'
)
```

---

### 8. **泰勒图** (`TaylorDiagram.py`)

模式与观测的统计对比工具。

**核心类**: `TaylorDiagram`

**使用示例**:
```python
import numpy as np
from wave_tools.TaylorDiagram import TaylorDiagram

# 参考标准差
refstd = obs_data.std()

# 创建泰勒图
dia = TaylorDiagram(refstd, fig=fig, rect=111)

# 添加模式样本点
for model in models:
    stddev = model_data[model].std()
    corrcoef = np.corrcoef(obs_data.flatten(), model_data[model].flatten())[0, 1]
    dia.add_sample(stddev, corrcoef, marker='o', label=model)

# 添加格网和RMS等值线
dia.add_grid()
dia.add_contours(levels=5)
```

---

### 9. **总体湿稳定度分析** (`GMS.py`)

计算Gross Moist Stability (GMS)，评估对流-大尺度环流耦合强度。

**主要函数**:
```python
# 水平GMS
h_gms = calc_horizontal_GMS(ta, zg, plev, lon, lat, ua, va, hus)

# 垂直GMS
v_gms = calc_vertical_GMS(ta, zg, plev, wa, hus, ua, va, lat, lon)

# 完整GMS计算
h_gms, v_gms = gross_moist_stability(
    ta_path, zg_path, ua_path, va_path, wa_path, hus_path
)
```

**物理意义**: 
- GMS > 0: 负反馈，抑制对流
- GMS < 0: 正反馈，促进对流组织化

---

## 📚 使用示例

### 示例 1: Wheeler-Kiladis频谱分析

```python
from wave_tools.Spectral_Analysis import SpectralAnalysis, Config
import xarray as xr

# 1. 配置参数
config = Config()
config.WINDOW_SIZE_DAYS = 96
config.WINDOW_SKIP_DAYS = 30

# 2. 创建分析对象
analysis = SpectralAnalysis(config)

# 3. 加载数据
analysis.load_data(
    data_path="path/to/olr_data.nc",
    variable="olr",
    lat_range=[-15, 15],
    time_range=('1997', '2014')
)

# 4. 预处理（去气候态、对称分解）
analysis.preprocess_data()

# 5. 计算功率谱
analysis.compute_spectra()

# 6. 平滑背景谱
analysis.smooth_background()

# 7. 绘制结果
analysis.plot_spectra(filename='wk_spectrum', output_path='./figures')

# 8. 保存频谱数据
analysis.save_spectra(output_path='./output/spectra.nc')
```

---

### 示例 2: 提取Kelvin波信号

```python
from wave_tools.Wave_Filter import WaveFilter
import xarray as xr

# 1. 加载数据
ds = xr.open_dataset('olr_daily.nc')['olr']
ds = ds.sel(lat=slice(-15, 15), time=slice('1997', '2014'))

# 2. 创建滤波器
wf = WaveFilter()

# 3. 提取Kelvin波（并行计算）
kelvin_signal = wf.extract_wave_signal(
    ds=ds,
    wave_name='kelvin',
    obs_per_day=1,
    use_parallel=True,
    n_jobs=-1,
    n_harm=3
)

# 4. 保存结果
kelvin_signal.to_netcdf('kelvin_filtered.nc')

# 5. 快速可视化
kelvin_signal.isel(time=100).plot(cmap='RdBu_r')
```

---

### 示例 3: 绘制带Matsuno模态的频谱

```python
from wave_tools.SpectrumPlotter import SpectrumPlotter
import xarray as xr
import numpy as np

# 加载频谱数据
ds = xr.open_dataset('spectrum_output.nc')
sym_spec = ds['power_symmetric'] / ds['background']

# 创建绘图对象
plotter = SpectrumPlotter(cpd_lines=[3, 6, 30])

# 绘图（自动叠加理论曲线）
plotter.plot_cmip_future(
    cmip_f=[sym_spec],
    wn=sym_spec.wavenumber.values,
    fq=sym_spec.frequency.values,
    figsize=(8, 6),
    max_wn_plot=15,
    max_freq_plot=0.5,
    matsuno_lines=True,
    he=[8, 25, 90],
    meridional_modes=[1],
    cmap='RdBu_r',
    contour_range=[1, 2, 0.1],
    labels=True,
    ofil='spectrum_with_matsuno'
)
```

---

### 示例 4: 计算波动相位

```python
from wave_tools.Wave_Phase import optimize_peak_detection, remove_clm
import xarray as xr

# 1. 加载滤波后的波动信号
kelvin = xr.open_dataset('kelvin_filtered.nc')['olr']

# 2. 去季节循环并低通滤波
kelvin_lowpass = remove_clm(kelvin, cutoff=1/10, order=4)

# 3. 纬度平均
kelvin_zonal = kelvin_lowpass.mean('lat')

# 4. 检测峰值和相位
V_peak, local_extrema = optimize_peak_detection(
    V=kelvin_zonal.values,
    kelvin_ref=kelvin_zonal,
    V_std=kelvin_zonal.std().values,
    Nstd=1.0,
    use_parallel=True,
    n_jobs=-1
)

# 5. 保存结果
V_peak.to_netcdf('kelvin_phase.nc')
```

---

## 🔬 API 文档

### `SpectralAnalysis` 类

**初始化**
```python
SpectralAnalysis(config=None)
```
- `config`: `Config` 对象，包含分析参数

**方法**

| 方法 | 说明 | 返回值 |
|-----|------|--------|
| `load_data()` | 加载netCDF数据 | self |
| `preprocess_data()` | 预处理（去趋势、去年循环、分解） | self |
| `compute_spectra()` | 计算波数-频率功率谱 | self |
| `smooth_background()` | 平滑背景谱 | self |
| `plot_spectra()` | 绘制频谱图 | None |
| `save_spectra()` | 保存频谱到netCDF | self |

---

### `WaveFilter` 类

**初始化**
```python
WaveFilter()
```

**方法**

| 方法 | 说明 | 返回值 |
|-----|------|--------|
| `extract_wave_signal()` | 提取特定波动信号 | xr.DataArray |
| `add_wave_param()` | 添加自定义波动参数 | None |
| `get_available_waves()` | 获取可用波动类型列表 | list |
| `get_wave_params()` | 获取波动参数 | dict |
| `check_filter_wave()` | 与NCL结果对比验证 | (fig, axes) |

---

### `SpectrumPlotter` 类

**初始化**
```python
SpectrumPlotter(cpd_lines, norm=None)
```
- `cpd_lines`: 要标注的周期线（天）
- `norm`: matplotlib的Normalize对象

**方法**

| 方法 | 说明 | 返回值 |
|-----|------|--------|
| `plot_cmip_future()` | 绘制频谱（支持多模式对比） | None |
| `set_ax()` | 设置坐标轴样式 | None |

---

## 📝 更新日志

### v1.0.0 (2025-10-17)
- ✨ 初始版本发布
- ✅ 完整的Wheeler-Kiladis频谱分析
- ✅ WK99波动滤波器实现
- ✅ Matsuno模态计算
- ✅ HEALPix网格转换（Numba加速）
- ✅ 波动相位检测
- ✅ 泰勒图工具
- ✅ GMS诊断分析

### 计划特性
- [ ] 自动化批处理脚本
- [ ] Jupyter Notebook教程
- [ ] 单元测试覆盖
- [ ] 完整文档（Read the Docs）
- [ ] MJO相位图绘制工具
- [ ] 波动传播动画生成

---

## 📧 联系方式

**作者**: Jianpu  
**Email**: xianpuji@hhu.edu.cn  
**机构**: Hohai University

**反馈与建议**: 欢迎通过Issue或Email联系作者

---

## 📄 许可证

MIT License

---

## 🙏 致谢

本工具包的开发参考了以下经典文献和工具：

- Wheeler, M., & Kiladis, G. N. (1999). Convectively coupled equatorial waves. *Journal of the Atmospheric Sciences*.
- Matsuno, T. (1966). Quasi-geostrophic motions in the equatorial area. *Journal of the Meteorological Society of Japan*.
- NCL (NCAR Command Language) - 波动滤波参考实现

---

**Happy Wave Hunting! 🌊**
