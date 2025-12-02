# Wave Tools - 热带大气波动分析工具包

> **简洁、清晰、实用的Python工具包，专注于热带大气波动分析**

**作者**: Jianpu | **邮箱**: xianpuji@hhu.edu.cn | **机构**: Hohai University

---

## 📦 安装

```bash
# 安装依赖
pip install -r wave_tools/requirements.txt

# 在Python中导入
from wave_tools import *
```

---

## 🎯 核心功能

本工具包提供**热带大气波动分析**的完整解决方案：

1. **Matsuno理论模态计算** - 计算赤道浅水方程的理论色散关系
2. **Wheeler-Kiladis频谱分析** - 波数-频率功率谱诊断
3. **波动滤波与提取** - 时空滤波提取Kelvin、ER、MRG等波动
4. **波动相位分析** - 峰值检测和相位合成分析
5. **EOF分析** - 支持SVD和xeofs两种方法的经验正交函数分解
6. **诊断工具** - 总体湿稳定度(GMS)、Radon变换、静力能计算等
7. **专业绘图** - WK频谱图、地图、泰勒图、CCKW包络图等

---

## 📁 模块结构

```
wave_tools/
├── matsuno.py       # Matsuno理论模态色散关系
├── spectral.py      # Wheeler-Kiladis频谱分析
├── filters.py       # 波动滤波与信号提取
├── phase.py         # 相位分析与峰值检测
├── eof.py           # EOF/PCA经验正交函数分解
├── diagnostics.py   # 大气动力学和热力学诊断
├── plotting.py      # 所有绘图功能
└── utils.py         # 工具函数（HEALPix转换、Radon变换等）
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

### 示例2: 提取Kelvin波（方法1：传统WaveFilter）

```python
from wave_tools.filters import WaveFilter

wf = WaveFilter()
kelvin = wf.extract_wave_signal(data, wave_name='kelvin', use_parallel=True)
kelvin.to_netcdf('kelvin_filtered.nc')
```

### 示例2.5: 提取Kelvin波（方法2：CCKWFilter - 推荐⭐）

```python
from wave_tools.filters import CCKWFilter
import xarray as xr

# 读取数据
pr_data = xr.open_dataarray('pr_data.nc')

# 初始化滤波器
wave_filter = CCKWFilter(
    ds=pr_data,
    sel_dict={'time': slice('1980-01-01', '1993-12-31'), 'lat': slice(-15, 15)},
    wave_name='kelvin',
    units='mm/day',
    spd=1,
    n_workers=4
)

# 方式1：逐步执行
wave_filter.load_data()
wave_filter.detrend_data()
wave_filter.fft_transform()
wave_filter.apply_filter()
wave_filter.inverse_fft()
filtered_data = wave_filter.create_output()

# 方式2：一键执行（推荐）
filtered_data = wave_filter.process()

# 计算标准差
std_data = filtered_data.std(dim='time')
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

**`WaveFilter`** - 波动滤波器（谐波分析方法）

| 方法 | 功能 | 参数说明 |
|------|------|----------|
| `extract_wave_signal(ds, wave_name, obs_per_day, use_parallel, n_jobs)` | **提取波动信号** | ds: 输入数据<br/>wave_name: 波动类型<br/>use_parallel: 是否并行 |
| `add_wave_param(wave_name, freq_range, wnum_range, equiv_depth)` | 添加自定义波动 | freq_range: 周期范围(天)<br/>wnum_range: 波数范围<br/>equiv_depth: 等效深度(m) |
| `get_available_waves()` | 获取可用波动列表 | 返回: ['kelvin', 'er', 'mrg', ...] |
| `get_wave_params(wave_name)` | 查看波动参数 | 返回: dict |

**`CCKWFilter`** - 对流耦合波动滤波器（频率-波数空间滤波，使用Dask并行）⭐

| 方法 | 功能 | 参数说明 |
|------|------|----------|
| `__init__(ds, var, sel_dict, wave_name, units, spd, n_workers, verbose)` | 初始化滤波器 | ds: 数据源<br/>var: 变量名（Dataset时需要）<br/>sel_dict: 选择字典<br/>wave_name: 'kelvin'或'er'<br/>spd: 每天采样次数<br/>n_workers: 并行工作进程数 |
| `load_data()` | 加载并预处理数据 | 自动应用sel_dict筛选 |
| `detrend_data()` | 数据去趋势 | 移除年际变化和线性趋势，应用Tukey窗口 |
| `fft_transform()` | 2D FFT变换 | 计算波数-频率网格 |
| `apply_filter()` | 应用滤波器 | 根据波动类型设置mask |
| `inverse_fft()` | 逆FFT变换 | 获取滤波后的实空间数据 |
| `create_output()` | 创建输出 | 返回带属性的xr.DataArray |
| `process()` | **一键执行完整流程** | 返回滤波后数据 |

**CCKWFilter特点**：
- ✅ 基于Wheeler-Kiladis频率-波数滤波方法
- ✅ 使用Dask进行大规模数据并行处理
- ✅ 自动应用浅水波色散关系约束
- ✅ 支持Kelvin波和ER波
- ✅ 提供详细的处理信息和诊断输出

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
| `optimize_peak_detection(V, kelvin_ref, V_std, Nstd, use_parallel, n_jobs)` | **检测显著峰值及影响范围** | V: 数据数组[time, lon]<br/>kelvin_ref: 参考数据<br/>V_std: 标准差<br/>Nstd: 显著性阈值倍数<br/>use_parallel: 是否并行<br/>n_jobs: 并行核数 | V_peak: 峰值影响范围<br/>local_extrema: 极值标识(1=最小值，-1=最大值) |
| `remove_clm(data, fs, cutoff, order)` | 去除季节循环并低通滤波 | data: xr.DataArray<br/>fs: 采样频率(Hz)<br/>cutoff: 截止频率<br/>order: Butterworth滤波器阶数 | 滤波后的异常场 |
| `find_local_extrema(V)` | 查找二维场的局部极值 | V: [time, lon]数组 | 标识数组(1=局部最小值，-1=局部最大值，NaN=其他) |
| `butter_lowpass_filter(data, cutoff_freq, fs, order)` | Butterworth低通滤波 | data: 1D时间序列<br/>cutoff_freq: 截止频率(Hz)<br/>fs: 采样频率<br/>order: 滤波器阶数 | 滤波后的序列 |
| `find_peak_influence_range(peak_idx, peak_value, zero_idx, peak_indices, V, V_std, Nstd)` | 确定单个峰值的时间影响范围 | peak_idx: 峰值索引<br/>zero_idx: 零点索引数组 | (left_idx, right_idx): 左右边界 |
| `process_single_longitude(ilon, V, local_min_max_id, V_std, Nstd)` | 处理单个经度的峰值检测 | ilon: 经度索引<br/>V: 原始数据<br/>local_min_max_id: 极值标记 | 1D峰值数组[time] |

**核心算法**: 基于零交叉点和局部极值的峰值影响范围识别，支持并行加速

**应用**: 波动相位合成分析、传播速度计算、生命周期统计、事件检测

---

### 5. eof.py - 经验正交函数分析

#### 主要类

**`EOFAnalyzer`** - EOF分析主类

| 参数 | 说明 |
|------|------|
| `method` | 分解方法: 'svd'(numpy实现) 或 'xeofs'(xeofs库) |
| `apply_land_mask` | 是否应用陆地/海洋掩膜 |
| `ocean_only` | True保留海洋点，False保留陆地点 |
| `mask_resolution` | 掩膜分辨率: 'c'(粗~50km), 'l'(低~10km), 'i'(中~2km), 'h'(高~400m), 'f'(全~100m) |

| 方法 | 功能 | 参数说明 |
|------|------|----------|
| `fit(data, n_modes, dim_names)` | 执行EOF分解 | data: xr.DataArray<br/>n_modes: 模态数<br/>dim_names: 维度名称dict |
| `plot_vertical_profiles(n_modes, figsize, save_path)` | 绘制垂直剖面 | n_modes: 绘制模态数 |
| `save_results(filename)` | 保存分析结果 | filename: 输出文件路径 |
| `load_results(filename)` | 加载已保存结果 | filename: 输入文件路径 |

**支持特性**:
- 自动维度检测和转置
- 气候态去除（谐波滤波）
- 陆地/海洋掩膜功能
- 解释方差计算
- 结果序列化存储

**示例**:
```python
analyzer = EOFAnalyzer(method='svd', apply_land_mask=True, ocean_only=True)
results = analyzer.fit(data, n_modes=4)
fig = analyzer.plot_vertical_profiles(n_modes=4)
```

---

### 6. diagnostics.py - 诊断工具

#### 主要函数

| 函数 | 功能 | 参数 | 返回 |
|------|------|------|------|
| `calc_horizontal_GMS(ta, zg, plev, lon, lat, ua, va, hus)` | **水平总体湿稳定度** | ta: 温度(K)<br/>zg: 位势高度(m)<br/>plev: 压力层(Pa)<br/>ua,va: 风场(m/s)<br/>hus: 比湿(kg/kg) | GMS (无量纲) |
| `calc_vertical_GMS(ta, zg, plev, wa, hus, ua, va, lat, lon)` | **垂直总体湿稳定度** | wa: 垂直速度(Pa/s) | GMS (无量纲) |
| `gross_moist_stability(ta_path, zg_path, ua_path, va_path, wa_path, hus_path)` | **完整GMS计算** | 各变量NetCDF文件路径 | h_gms, v_gms |
| `calc_dse(ta, zg, plev)` | 计算干静力能 | ta: 温度<br/>zg: 位势高度 | DSE (J/kg) |
| `compute_dx_dy(lat, lon)` | 计算网格间距 | lat, lon: 纬度经度数组 | dx, dy (米) |
| `vertically_integrated_moist_flux_divergence(hus, ua, va, lat, lon)` | 垂直积分的湿通量散度 | hus: 比湿<br/>ua, va: 风场 | 能量通量散度 (W/m²) |

**物理常数**:
- Cp = 1004.0 J/(kg·K) - 定压比热
- g = 9.81 m/s² - 重力加速度
- L = 2.5×10⁶ J/kg - 蒸发潜热
- T_ref = 300.0 K - 参考温度

**计算公式**:
- DSE = Cp·T + g·Z
- 水平GMS = -T_ref · ∫(V·∇s)dp / ∫L·∇·(rV)dp
- 垂直GMS = -T_ref · ∫(ω·∂s/∂p)dp / ∫L·∇·(rV)dp

**物理意义**: 
- GMS > 0: 负反馈，抑制对流发展
- GMS < 0: 正反馈，促进对流组织化

**依赖库**: metpy, geocat-comp (可选，提供简化实现)

---

### 7. plotting.py - 绘图功能

#### CCKW包络图

| 函数 | 功能 | 参数 |
|------|------|------|
| `get_cckw_envelope_curve(he, fmax)` | 计算CCKW包络曲线坐标 | he: 等效深度列表(m)<br/>fmax: 最大频率列表 |
| `plot_cckw_envelope(he, fmax, save_path, dpi)` | 绘制CCKW包络示意图 | save_path: 保存路径<br/>dpi: 分辨率(默认200) |

#### Wheeler-Kiladis频谱图

| 函数 | 功能 | 关键参数 |
|------|------|----------|
| `plot_wk_spectrum(power_symmetric, power_antisymmetric, background, wavenumber, frequency, ...)` | **绘制WK频谱** | max_wn: 最大波数<br/>max_freq: 最大频率<br/>add_matsuno_lines: 是否添加理论曲线<br/>he: 等效深度列表<br/>cpd_lines: 标注周期线<br/>cmap: 色标<br/>levels: 等值线水平 |
| `set_axis_for_wave(ax, text_size, freq_lines, depth, cpd_lines, max_wn_plot, max_freq_plot)` | 设置波动分析坐标轴 | freq_lines: 是否显示频率线<br/>depth: 是否显示等效深度线<br/>cpd_lines: 周期线列表(天) |

**默认设置**:
- 频谱等值线: [1, 1.2, 1.4, 1.6, 1.8, 2.0]
- 周期标注: [3, 6, 30]天
- 等效深度: [8, 25, 90]米

#### 地图绘图

| 函数 | 功能 | 参数 |
|------|------|------|
| `setup_map_axes(ax, title, box)` | 设置地图坐标轴 | box: [lon_min, lon_max, lat_min, lat_max] |
| `plot_spatial_field(data, ax, cmap, title, box, levels, add_coastlines)` | 绘制空间场 | data: xr.DataArray<br/>cmap: 色标<br/>add_coastlines: 是否添加海岸线 |

**依赖**: cartopy (可选，用于地图投影)

#### 泰勒图

**类**: `TaylorDiagram(refstd, fig, rect, srange, extend)`

| 方法 | 功能 | 参数 |
|------|------|------|
| `add_sample(stddev, corrcoef, *args, **kwargs)` | 添加模式样本点 | stddev: 标准差<br/>corrcoef: 相关系数 |
| `add_contours(levels, **kwargs)` | 添加RMS等值线 | levels: 等值线水平数 |

#### 通用工具

| 函数 | 功能 | 参数 |
|------|------|------|
| `save_figure(fig, filename, folder, fmt, dpi)` | 保存图形 | fmt: 'pdf'/'png'/'jpg'<br/>dpi: 分辨率(默认600) |

**色标支持**: 默认使用cmaps.amwg_blueyellowred，回退到'RdBu_r'

---

### 8. utils.py - 工具函数

#### 数据处理

| 函数 | 功能 | 参数 | 返回 |
|------|------|------|------|
| `load_data(path, var, lat_range)` | 加载NetCDF数据 | path: 文件路径<br/>var: 变量名<br/>lat_range: 纬度范围(默认-15,15) | data, lon, lat |
| `filter_series(series, min_wn, max_wn)` | 波数范围过滤 | series: 时间序列<br/>min_wn, max_wn: 波数边界 | 过滤后的序列 |
| `filter_paths_by_models(paths, model_names, loc, sep, case_sensitive, strict, missing_ok)` | **按模式名过滤文件路径** | paths: 文件路径列表<br/>model_names: 模式名列表<br/>loc: 分隔后的位置<br/>sep: 分隔符(默认'_')<br/>case_sensitive: 大小写敏感<br/>strict: 严格模式 | 匹配的路径列表 |
| `extract_model_name(path, loc, sep)` | 从路径提取模式名 | path: 文件路径<br/>loc: 分隔后的位置(默认1)<br/>sep: 分隔符(默认'_') | 模式名字符串 |

#### HEALPix网格转换

| 函数 | 功能 | 参数 | 返回 |
|------|------|------|------|
| `dataarray_healpix_to_equatorial_latlon(healpix_dataarray, nside, nest, minmax_lat)` | **HEALPix→等经纬度** | nside: HEALPix参数<br/>nest: 是否nested排序<br/>minmax_lat: 纬度范围 | 等经纬度DataArray |
| `get_region_healpix_(zoom, extent, nest)` | 获取区域网格索引 | extent: [lon_min, lon_max, lat_min, lat_max] | 网格索引数组 |

**优化**: 使用Numba @jit加速，速度提升10倍+

**依赖**: healpy (可选)

#### Radon变换

| 函数 | 功能 | 参数 | 返回 |
|------|------|------|------|
| `calc_radon_angle(field, theta_range)` | **计算传播角度** | field: 2D时空数组<br/>theta_range: 角度范围(默认0-180度) | theta, intensity, theta_max |
| `calc_c_from_theta(theta_deg, dx_deg, dt_sec, lat)` | 从角度计算相速度 | theta_deg: 角度(度)<br/>dx_deg: 经度间隔<br/>dt_sec: 时间间隔<br/>lat: 纬度 | 速度(m/s) |
| `plot_radon_energy_distribution(theta, energy, title, color, ax)` | 绘制Radon能量分布 | energy: Radon能量<br/>ax: matplotlib轴(可选) | theta_max, theta_ci (95%置信区间) |

**应用**: 
- 波动传播方向诊断
- 相速度估算
- 东传/西传波动识别

**物理原理**: Radon变换将时空场投影到不同角度，能量最大处对应主导传播方向

#### 色标和字体工具

| 函数 | 功能 | 参数 |
|------|------|------|
| `create_cmap_from_string(color_string)` | 从字符串创建色标 | color_string: 颜色列表(十六进制或RGB) | ListedColormap对象 |
| `set_matplotlib_font(font_dir, arial_font_path)` | 设置Matplotlib字体 | font_dir: 字体目录<br/>arial_font_path: Arial字体路径 | None |
| `save_figure(fig, filename, folder, fmt, dpi)` | 保存图形 | fig: matplotlib图形对象<br/>folder: 保存目录<br/>fmt: 格式<br/>dpi: 分辨率 | None |

#### 其他工具

| 函数 | 功能 |
|------|------|
| `get_curve(he, fmax)` | 获取Kelvin波包络曲线（与plotting.py中的get_cckw_envelope_curve相同） |

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
from wave_tools.eof import EOFAnalyzer

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
kelvin = wf.extract_wave_signal(data, wave_name='kelvin', use_parallel=True, n_harm=3)
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

# ===== 4. EOF分析 =====
analyzer = EOFAnalyzer(method='svd', apply_land_mask=True, ocean_only=True)
results = analyzer.fit(data, n_modes=4)
fig = analyzer.plot_vertical_profiles(n_modes=4, save_path='eof_profiles.png')
analyzer.save_results('eof_results.pkl')

# ===== 5. 诊断分析 =====
h_gms, v_gms = gross_moist_stability(
    'ta.nc', 'zg.nc', 'ua.nc', 'va.nc', 'wa.nc', 'hus.nc'
)

print(f"水平GMS: {h_gms.mean().values:.3f}")
print(f"垂直GMS: {v_gms.mean().values:.3f}")
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
joblib >= 1.0.0       # 并行计算
```

### 可选依赖
```
cmaps                 # NCL气象色标
cartopy              # 地图投影
healpy               # HEALPix球面网格处理
numba                # JIT加速计算
metpy                # 气象计算工具
geocat-comp          # NCAR地球科学计算
xeofs                # EOF分析（可选方法）
global-land-mask     # 陆地/海洋掩膜
scikit-image         # Radon变换
```

**安装全部依赖**:
```bash
pip install numpy xarray scipy matplotlib pandas joblib
pip install cmaps cartopy healpy numba metpy geocat-comp xeofs global-land-mask scikit-image
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

**版本**: 1.1.0 | **最后更新**: 2025-11-04
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
joblib >= 1.0.0       # 并行计算
```

### 可选依赖
```
cmaps                 # NCL气象色标
cartopy              # 地图投影
healpy               # HEALPix球面网格处理
numba                # JIT加速计算
metpy                # 气象计算工具
geocat-comp          # NCAR地球科学计算
xeofs                # EOF分析（可选方法）
global-land-mask     # 陆地/海洋掩膜
scikit-image         # Radon变换
```

**安装全部依赖**:
```bash
pip install numpy xarray scipy matplotlib pandas joblib
pip install cmaps cartopy healpy numba metpy geocat-comp xeofs global-land-mask scikit-image
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

**版本**: 1.1.0 | **最后更新**: 2025-11-04
