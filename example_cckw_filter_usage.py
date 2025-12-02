#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CCKWFilter 使用示例
==================

演示如何使用 wave_tools.filters.CCKWFilter 进行Kelvin波滤波

作者: xpji
邮箱: xianpuji@hhu.edu.cn
日期: 2025-04-09
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os

# 导入CCKWFilter
from wave_tools.filters import CCKWFilter

# ============================================================================
# 示例1: 基本使用 - 处理单个数据集
# ============================================================================

def example1_basic_usage():
    """基本使用示例"""
    print("="*70)
    print("示例1: 基本使用")
    print("="*70)
    
    # 读取数据
    DATA_DIR = "./processed_data"
    pr_file = os.path.join(DATA_DIR, 'pr_cntl_2deg_interp.nc')
    pr_data = xr.open_dataarray(pr_file) * 86400  # 转换为 mm/day
    
    # 初始化滤波器
    wave_filter = CCKWFilter(
        ds=pr_data,
        sel_dict={
            'time': slice('1980-01-01', '1993-12-31'),
            'lat': slice(-15, 15)
        },
        wave_name='kelvin',
        units='mm/day',
        spd=1,
        n_workers=4,
        verbose=True
    )
    
    # 打印滤波器信息
    print(wave_filter)
    
    # 执行滤波（方式1：逐步执行）
    print("\n方式1: 逐步执行各个步骤")
    wave_filter.load_data()
    wave_filter.detrend_data()
    wave_filter.fft_transform()
    wave_filter.apply_filter()
    wave_filter.inverse_fft()
    filtered_data = wave_filter.create_output()
    
    print(f"\n✅ 滤波完成！")
    print(f"   输出数据形状: {filtered_data.shape}")
    print(f"   数值范围: [{filtered_data.min().values:.3f}, {filtered_data.max().values:.3f}]")
    
    return filtered_data


def example2_one_step():
    """一键执行示例"""
    print("\n" + "="*70)
    print("示例2: 一键执行")
    print("="*70)
    
    # 读取数据
    DATA_DIR = "./processed_data"
    pr_file = os.path.join(DATA_DIR, 'pr_cntl_2deg_interp.nc')
    pr_data = xr.open_dataarray(pr_file) * 86400
    
    # 初始化并一键执行
    wave_filter = CCKWFilter(
        ds=pr_data,
        sel_dict={
            'time': slice('1980-01-01', '1993-12-31'),
            'lat': slice(-15, 15)
        },
        wave_name='kelvin',
        units='mm/day',
        spd=1,
        n_workers=4
    )
    
    # 一键执行
    filtered_data = wave_filter.process()
    
    return filtered_data


# ============================================================================
# 示例3: 批量处理多个实验
# ============================================================================

def example3_batch_processing():
    """批量处理多个实验"""
    print("\n" + "="*70)
    print("示例3: 批量处理多个实验")
    print("="*70)
    
    DATA_DIR = "./processed_data"
    
    # 定义实验
    experiments = {
        'CNTL': 'pr_cntl_2deg_interp.nc',
        'P4K': 'pr_4k_2deg_interp.nc',
        '4CO2': 'pr_4co2_2deg_interp.nc'
    }
    
    # 存储结果
    kelvin_wave_data = {}
    kelvin_std = {}
    
    # 批量处理
    for exp_name, filename in experiments.items():
        print(f"\n📍 处理 {exp_name}...")
        
        # 读取数据
        pr_data = xr.open_dataarray(os.path.join(DATA_DIR, filename)) * 86400
        
        # 滤波
        wave_filter = CCKWFilter(
            ds=pr_data,
            sel_dict={
                'time': slice('1980-01-01', '1993-12-31'),
                'lat': slice(-15, 15)
            },
            wave_name='kelvin',
            units='mm/day',
            spd=1,
            n_workers=4,
            verbose=False  # 关闭详细输出
        )
        
        filtered_data = wave_filter.process()
        kelvin_wave_data[exp_name] = filtered_data
        
        # 计算标准差
        std_data = filtered_data.std(dim='time')
        kelvin_std[exp_name] = std_data
        
        print(f"   ✅ {exp_name} 完成！STD范围: "
              f"[{std_data.min().values:.3f}, {std_data.max().values:.3f}] mm/day")
    
    return kelvin_wave_data, kelvin_std


# ============================================================================
# 示例4: 提取ER波
# ============================================================================

def example4_er_wave():
    """提取ER波示例"""
    print("\n" + "="*70)
    print("示例4: 提取ER波")
    print("="*70)
    
    DATA_DIR = "./processed_data"
    pr_file = os.path.join(DATA_DIR, 'pr_cntl_2deg_interp.nc')
    pr_data = xr.open_dataarray(pr_file) * 86400
    
    # 提取ER波
    wave_filter = CCKWFilter(
        ds=pr_data,
        sel_dict={
            'time': slice('1980-01-01', '1993-12-31'),
            'lat': slice(-15, 15)
        },
        wave_name='er',  # 改为ER波
        units='mm/day',
        spd=1,
        n_workers=4
    )
    
    er_data = wave_filter.process()
    
    print(f"\n✅ ER波提取完成！")
    print(f"   ER波参数:")
    print(f"   - 周期范围: {wave_filter.tMin}-{wave_filter.tMax} 天")
    print(f"   - 波数范围: {wave_filter.kmin} 到 {wave_filter.kmax}")
    print(f"   - 等效深度: {wave_filter.hmin}-{wave_filter.hmax} 米")
    
    return er_data


# ============================================================================
# 示例5: 保存结果
# ============================================================================

def example5_save_results():
    """保存结果示例"""
    print("\n" + "="*70)
    print("示例5: 保存滤波结果")
    print("="*70)
    
    # 获取滤波数据
    kelvin_wave_data, kelvin_std = example3_batch_processing()
    
    # 保存目录
    OUTPUT_DIR = "./processed_data"
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # 保存滤波后的数据
    for exp_name, data in kelvin_wave_data.items():
        output_file = os.path.join(OUTPUT_DIR, 
                                   f'pr_kelvin_wave_{exp_name.lower()}_1980_1993.nc')
        data.to_netcdf(output_file)
        print(f"✅ 保存: {output_file}")
    
    # 保存标准差数据
    for exp_name, std_data in kelvin_std.items():
        std_file = os.path.join(OUTPUT_DIR, 
                               f'pr_kelvin_wave_std_{exp_name.lower()}_1980_1993.nc')
        std_data.to_netcdf(std_file)
        print(f"✅ 保存: {std_file}")
    
    print(f"\n✅ 所有数据已保存到 {OUTPUT_DIR}")


# ============================================================================
# 主函数
# ============================================================================

if __name__ == "__main__":
    print("""
    ╔══════════════════════════════════════════════════════════════════╗
    ║            CCKWFilter 使用示例集合                                 ║
    ║                                                                  ║
    ║  本脚本演示了 wave_tools.filters.CCKWFilter 的多种使用方式         ║
    ╚══════════════════════════════════════════════════════════════════╝
    """)
    
    # 运行示例（根据需要注释/取消注释）
    
    # 示例1: 基本使用
    # filtered_data = example1_basic_usage()
    
    # 示例2: 一键执行
    # filtered_data = example2_one_step()
    
    # 示例3: 批量处理
    kelvin_wave_data, kelvin_std = example3_batch_processing()
    
    # 示例4: ER波
    # er_data = example4_er_wave()
    
    # 示例5: 保存结果
    # example5_save_results()
    
    print("\n" + "="*70)
    print("✅ 所有示例运行完成！")
    print("="*70)
