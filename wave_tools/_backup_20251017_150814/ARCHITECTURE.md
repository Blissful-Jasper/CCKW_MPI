# Wave Tools 模块架构图

```
wave_tools/
│
├── 📦 核心分析层 (Core Analysis Layer)
│   │
│   ├── 🌊 频谱分析 (Spectral Analysis)
│   │   ├── Spectral_Analysis.py         # WK谱计算（规则网格）
│   │   ├── spectral_analysis_icon.py    # WK谱计算（非结构网格）
│   │   └── SpectrumPlotter.py           # 频谱可视化
│   │
│   ├── 🔍 波动提取 (Wave Extraction)
│   │   ├── Wave_Filter.py               # 时空滤波器
│   │   └── Wave_Phase.py                # 相位检测
│   │
│   └── 📐 理论模态 (Theoretical Modes)
│       ├── ma.py                        # Matsuno模态计算
│       └── core.py                      # CCKW包络工具
│
├── 🛠️ 工具层 (Utilities Layer)
│   │
│   ├── 📊 数据处理 (Data Processing)
│   │   └── utils.py
│   │       ├── load_data()              # 数据加载
│   │       ├── filter_series()          # 序列过滤
│   │       └── filter_paths_by_models() # 文件过滤
│   │
│   ├── 🗺️ 网格转换 (Grid Conversion)
│   │   └── utils.py
│   │       ├── dataarray_healpix_to_equatorial_latlon()  # HEALPix转换
│   │       ├── get_region_healpix_()                     # 区域提取
│   │       └── _fast_periodic_interp()                   # 周期性插值(Numba加速)
│   │
│   ├── 📏 诊断工具 (Diagnostic Tools)
│   │   └── utils.py
│   │       ├── calc_radon_angle()       # Radon变换
│   │       ├── calc_c_from_theta()      # 相速度计算
│   │       └── plot_radon_energy_distribution()  # 能量分布
│   │
│   └── 🎨 图形工具 (Graphics Utilities)
│       └── utils.py
│           ├── save_figure()            # 图像保存
│           ├── create_cmap_from_string()  # 自定义色标
│           └── get_curve()              # CCKW曲线坐标
│
├── 📈 可视化层 (Visualization Layer)
│   │
│   ├── 🗺️ 地图绘图 (Map Plotting)
│   │   └── plot.py
│   │       ├── make_space_fig()         # 地图设置
│   │       ├── plot_space_data()        # 空间数据绘制
│   │       └── plot_multiple_wave_trends()  # 趋势分析
│   │
│   └── 📊 统计图表 (Statistical Plots)
│       └── TaylorDiagram.py             # 泰勒图
│           ├── add_sample()             # 添加样本点
│           ├── add_grid()               # 添加网格
│           └── add_contours()           # 添加RMS等值线
│
├── 🔬 诊断分析层 (Diagnostic Layer)
│   │
│   └── GMS.py                           # 总体湿稳定度
│       ├── calc_horizontal_GMS()        # 水平GMS
│       ├── calc_vertical_GMS()          # 垂直GMS
│       └── gross_moist_stability()      # 完整GMS计算
│
└── 🔧 辅助模块 (Auxiliary Modules)
    │
    └── functions.py                     # 频谱辅助函数
        └── Spectrum类                   # 频谱类（早期版本）


═══════════════════════════════════════════════════════════════════

📋 数据流图 (Data Flow)

原始数据 (NetCDF)
    │
    ├─→ [Spectral_Analysis] → 功率谱 → [SpectrumPlotter] → WK谱图
    │        │
    │        └─→ 对称/反对称分量
    │
    ├─→ [Wave_Filter] → 滤波信号 → [Wave_Phase] → 相位信息
    │        │
    │        └─→ Kelvin/ER/MRG/IG/MJO/TD
    │
    └─→ [GMS] → 湿稳定度 → 对流-环流耦合诊断


═══════════════════════════════════════════════════════════════════

🔗 模块依赖关系 (Dependencies)

核心依赖:
    numpy, xarray, scipy, matplotlib, pandas

频谱分析:
    cmaps (气象色标)
    scipy.fft (傅里叶变换)

波动滤波:
    joblib (并行计算)
    numba (加速计算)

可视化:
    cartopy (地图投影)
    TaylorDiagram (统计图)

诊断工具:
    metpy (气象计算)
    healpy (HEALPix处理)


═══════════════════════════════════════════════════════════════════

📊 典型使用场景 (Use Cases)

场景1: 波动谱诊断
    输入: OLR日数据 → Spectral_Analysis → 输出: WK频谱图
    
场景2: 波动信号提取
    输入: OLR日数据 → Wave_Filter → 输出: Kelvin波时间序列
    
场景3: 模式评估
    输入: 多模式数据 → Spectral_Analysis + TaylorDiagram → 输出: 评估报告
    
场景4: 相位合成
    输入: 滤波信号 → Wave_Phase → 输出: 波动生命周期演变
    
场景5: 湿稳定度诊断
    输入: 3D气象场 → GMS → 输出: 水平/垂直GMS


═══════════════════════════════════════════════════════════════════

🎯 设计原则 (Design Principles)

1. 模块化: 每个模块独立功能明确
2. 可扩展: 易于添加新波动类型和分析方法
3. 高性能: 支持并行计算和向量化操作
4. 用户友好: 清晰的API和丰富的文档
5. 科学严谨: 基于经典文献和验证方法
```

---

**说明**:
- 📦 = 主要模块分组
- 🌊 = 波动分析相关
- 🛠️ = 工具函数
- 📈 = 可视化功能
- 🔬 = 诊断分析
- 🔧 = 辅助功能
