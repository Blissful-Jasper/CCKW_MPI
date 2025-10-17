# 快速开始指南

## 1. 安装

```bash
# 克隆项目
git clone <your-repo-url>
cd wave_tools

# 安装依赖
pip install -r requirements.txt
```

## 2. 五分钟入门

### 例子1: 计算Wheeler-Kiladis频谱

```python
from wave_tools import SpectralAnalysis, Config
import xarray as xr

# 加载数据
data = xr.open_dataset('your_data.nc')['olr']

# 运行分析
config = Config()
analysis = SpectralAnalysis(config)
analysis.load_data(data=data)
analysis.preprocess_data()
analysis.compute_spectra()
analysis.smooth_background()
analysis.plot_spectra()
```

### 例子2: 提取Kelvin波

```python
from wave_tools import WaveFilter
import xarray as xr

# 加载数据
data = xr.open_dataset('your_data.nc')['olr']

# 提取Kelvin波
wf = WaveFilter()
kelvin = wf.extract_wave_signal(data, wave_name='kelvin')

# 保存结果
kelvin.to_netcdf('kelvin_filtered.nc')
```

### 例子3: 绘制带理论曲线的频谱

```python
from wave_tools import SpectrumPlotter
import xarray as xr

# 加载频谱数据
ds = xr.open_dataset('spectrum.nc')
sym_spec = ds['power_symmetric'] / ds['background']

# 绘图
plotter = SpectrumPlotter(cpd_lines=[3, 6, 30])
plotter.plot_cmip_future(
    cmip_f=[sym_spec],
    wn=sym_spec.wavenumber.values,
    fq=sym_spec.frequency.values,
    max_wn_plot=15,
    max_freq_plot=0.5,
    matsuno_lines=True,
    ofil='spectrum_plot'
)
```

## 3. 数据格式要求

### 输入数据格式
- **文件类型**: netCDF (.nc)
- **必需维度**: time, lat, lon
- **维度顺序**: (time, lat, lon) 或 (time, lon, lat)
- **时间分辨率**: 支持日数据、6小时数据等
- **纬度范围**: 建议包含 ±15° 以覆盖热带区域

### 示例数据结构
```python
<xarray.DataArray 'olr' (time: 6575, lat: 31, lon: 144)>
Coordinates:
  * time     (time) datetime64[ns] 1997-01-01 ... 2014-12-31
  * lat      (lat) float32 -15.0 -14.0 -13.0 ... 13.0 14.0 15.0
  * lon      (lon) float32 0.0 2.5 5.0 ... 352.5 355.0 357.5
Attributes:
    units:      W/m^2
    long_name:  Outgoing Longwave Radiation
```

## 4. 常见问题

### Q1: 如何处理缺失值？
```python
# 在加载数据后
data = data.interpolate_na(dim='time', method='linear')
```

### Q2: 内存不足怎么办？
```python
# 使用分块加载
data = xr.open_dataset('large_file.nc', chunks={'time': 100})

# 或者减小分析窗口
config.WINDOW_SIZE_DAYS = 64  # 默认96天
```

### Q3: 如何加速计算？
```python
# 使用并行计算
wf = WaveFilter()
kelvin = wf.extract_wave_signal(
    data, 
    wave_name='kelvin',
    use_parallel=True,
    n_jobs=-1  # 使用所有CPU核心
)
```

## 5. 推荐工作流程

```
1. 数据准备
   ├── 下载原始数据 (OLR/降水等)
   ├── 格式转换为netCDF
   └── 质量控制（缺失值、异常值）

2. 频谱诊断
   ├── 计算Wheeler-Kiladis频谱
   ├── 识别主要波动信号
   └── 与理论对比

3. 波动提取
   ├── 基于频谱结果选择波动类型
   ├── 应用时空滤波
   └── 验证滤波效果

4. 深入分析
   ├── 相位合成分析
   ├── 传播特征
   ├── 波动强度趋势
   └── 湿稳定度诊断

5. 可视化与发表
   ├── 生成高质量图件
   ├── 数据归档
   └── 撰写论文/报告
```

## 6. 更多资源

- 📚 完整文档: 见 `README.md`
- 📧 技术支持: xianpuji@hhu.edu.cn
- 🐛 问题反馈: 提交 GitHub Issue
- 📖 参考文献: 见 `REFERENCES.md`

---

祝研究顺利！🎉
