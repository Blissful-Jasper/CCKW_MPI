"""
Wave Tools - 热带大气波动分析工具包
====================================

一个用于热带大气波动分析的综合Python工具包，支持：
- Wheeler-Kiladis频谱分析
- 波动滤波与提取（Kelvin, ER, MRG, IG, MJO, TD）
- Matsuno模态理论计算
- 波动相位分析
- 专业科学可视化

作者: Jianpu
联系: xianpuji@hhu.edu.cn
机构: Hohai University
"""

__version__ = "1.0.0"
__author__ = "Jianpu"
__email__ = "xianpuji@hhu.edu.cn"

# 核心模块
from .core import *
from .ma import *
from .utils import *

# 频谱分析
from .Spectral_Analysis import SpectralAnalysis, Config
from .SpectrumPlotter import SpectrumPlotter

# 波动分析
from .Wave_Filter import WaveFilter
from .Wave_Phase import optimize_peak_detection, remove_clm

# 可视化
from .plot import *
from .TaylorDiagram import TaylorDiagram

# 诊断工具
from .GMS import gross_moist_stability

__all__ = [
    # 核心功能
    'SpectralAnalysis',
    'Config',
    'WaveFilter',
    'SpectrumPlotter',
    
    # 工具函数
    'optimize_peak_detection',
    'remove_clm',
    'gross_moist_stability',
    
    # 可视化
    'TaylorDiagram',
]