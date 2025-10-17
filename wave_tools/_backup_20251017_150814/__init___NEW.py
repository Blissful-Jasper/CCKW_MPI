"""
Wave Tools - 热带大气波动分析工具包
====================================

简洁、清晰、实用的Python工具包，专注于热带大气波动分析

主要功能:
---------
1. Matsuno理论模态计算
2. Wheeler-Kiladis频谱分析
3. 波动滤波与提取
4. 波动相位分析
5. 诊断工具（GMS等）
6. 专业科学绘图

作者: Jianpu
邮箱: xianpuji@hhu.edu.cn
机构: Hohai University
"""

__version__ = "1.0.0"
__author__ = "Jianpu"
__email__ = "xianpuji@hhu.edu.cn"

# ===== 核心分析模块 =====

# Matsuno理论
from .matsuno import (
    kelvin_mode,
    er_n,
    mrg_mode,
    eig_n,
    wig_n,
    matsuno_modes_wk,
    matsuno_dataframe,
)

# 频谱分析
from .spectral import (
    WKSpectralAnalysis,
    SpectralConfig,
    calculate_wk_spectrum,
)

# 波动滤波
from .filters import WaveFilter

# 相位分析
from .phase import (
    optimize_peak_detection,
    remove_clm,
    find_local_extrema,
    butter_lowpass_filter,
)

# ===== 诊断工具 =====

from .diagnostics import (
    calc_horizontal_GMS,
    calc_vertical_GMS,
    gross_moist_stability,
)

# ===== 绘图功能 =====

from .plotting import (
    # CCKW包络
    get_cckw_envelope_curve,
    plot_cckw_envelope,
    
    # WK频谱
    plot_wk_spectrum,
    
    # 地图
    setup_map_axes,
    plot_spatial_field,
    
    # 泰勒图
    TaylorDiagram,
    
    # 工具
    save_figure,
)

# ===== 工具函数 =====

from .utils import (
    # 数据处理
    load_data,
    filter_series,
    filter_paths_by_models,
    extract_model_name,
    
    # HEALPix
    dataarray_healpix_to_equatorial_latlon,
    get_region_healpix_,
    
    # Radon
    calc_radon_angle,
    calc_c_from_theta,
    plot_radon_energy_distribution,
    
    # 其他
    create_cmap_from_string,
)


# ===== 导出列表 =====

__all__ = [
    # 核心类
    'WKSpectralAnalysis',
    'SpectralConfig',
    'WaveFilter',
    'TaylorDiagram',
    
    # 频谱分析
    'calculate_wk_spectrum',
    
    # Matsuno模态
    'kelvin_mode',
    'er_n',
    'mrg_mode',
    'matsuno_modes_wk',
    
    # 波动滤波（通过WaveFilter类）
    
    # 相位分析
    'optimize_peak_detection',
    'remove_clm',
    
    # 诊断
    'gross_moist_stability',
    'calc_horizontal_GMS',
    'calc_vertical_GMS',
    
    # 绘图
    'plot_wk_spectrum',
    'plot_cckw_envelope',
    'plot_spatial_field',
    'save_figure',
    
    # 工具
    'load_data',
    'filter_series',
    'calc_radon_angle',
]


# ===== 模块信息 =====

def get_version():
    """返回版本号"""
    return __version__

def list_available_waves():
    """列出可用的波动类型"""
    return ['kelvin', 'er', 'mrg', 'ig', 'mjo', 'td']

def print_info():
    """打印工具包信息"""
    info = f"""
    Wave Tools v{__version__}
    ========================
    热带大气波动分析工具包
    
    作者: {__author__}
    邮箱: {__email__}
    
    主要模块:
    - matsuno.py       Matsuno理论模态
    - spectral.py      频谱分析
    - filters.py       波动滤波
    - phase.py         相位分析
    - diagnostics.py   诊断工具
    - plotting.py      绘图功能
    - utils.py         工具函数
    
    快速开始:
    >>> from wave_tools import *
    >>> help(WKSpectralAnalysis)
    >>> help(WaveFilter)
    """
    print(info)
