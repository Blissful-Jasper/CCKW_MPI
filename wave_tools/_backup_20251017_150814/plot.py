import numpy as np
import matplotlib.ticker as ticker
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.stats import linregress
from typing import List, Optional, Tuple, Union
import os
import matplotlib.pyplot as plt
import xarray as xr

def make_space_fig(ax, title: str, box: Optional[List[float]] = None):
    """
    设置地图子图 (ax) 的刻度、格式和标题，并返回修改后的 ax。

    输入:
        ax: cartopy/matplotlib 的子图对象 (Axes)。
        title (str): 子图标题。
        box (Optional[List[float]]): 可选，地图的经纬度范围 [lon_min, lon_max, lat_min, lat_max]。
                                      如果未提供，默认使用 [0, 360, -20, 20]。

    输出:
        ax: 返回设置好的子图对象。

    示例:
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
        makefig(ax, title='Example Map', box=[0, 360, -30, 30])
        plt.show()
    """
    if box is None:
        box = [0, 360, -20, 20]
    
    # 设置经纬度刻度
    ax.set_xticks(np.linspace(box[0], box[1], 7), crs=ccrs.PlateCarree())
    ax.set_yticks(np.linspace(box[2], box[3], 5), crs=ccrs.PlateCarree())

    # 设置格式
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())

    # 设置主刻度和副刻度
    ax.xaxis.set_major_locator(ticker.MultipleLocator((box[1] - box[0]) / 6))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(6))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))

    # 标题和刻度样式
    ax.set_title(title, loc='left')
    ax.tick_params(which='major', length=5)

    return ax


def plot_space_data(
    data,
    ax,
    cmap,
    fmt: str,
    title: str,
    box: Optional[List[float]] = None,
    levels: Optional[np.ndarray] = None
):
    """
    在给定的子图 (ax) 上绘制等值填色图 (contourf)，并进行地图美化设置。

    输入:
        data: 待绘制的 xarray.DataArray 或类似对象，通常是二维地理数据。
        ax: matplotlib/cartopy 的子图对象 (Axes)。
        cmap: matplotlib 颜色映射 (Colormap)。
        fmt (str): 数据单位或说明，将用于 colorbar 标注或图示说明。
        title (str): 子图标题。
        box (Optional[List[float]]): 地图显示范围 [lon_min, lon_max, lat_min, lat_max]，默认 [0, 360, -20, 20]。
        levels (Optional[np.ndarray]): 等值线的划分等级，默认是 np.linspace(0, 4, 41)。

    输出:
        f: 返回绘制的 contourf 图层对象 (QuadContourSet)。

    示例:
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
        plot_data(data, ax, cmap='viridis', fmt='mm/day', title='Rainfall', box=[0, 360, -30, 30])
        plt.show()
    """
    if box is None:
        box = [0, 360, -20, 20]
    
    if levels is None:
        levels = np.linspace(0., 4., 41)
    
    # 绘制填色图
    f = data.plot.contourf(
        ax=ax,
        cmap=cmap,
        add_labels=False,
        levels=levels,
        add_colorbar=False,
        transform=ccrs.PlateCarree()
    )

    # 地图基本设置
    ax.coastlines()
    ax.set_extent(box, crs=ccrs.PlateCarree())
    ax.set_aspect(3)
    
    # 调用辅助函数设置标题与刻度（需存在 makefig 函数）
    make_space_fig(ax, title, box)
    
    # 美化子图边框
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)

    return f

def plot_multiple_wave_trends(
    wave_filters: List[xr.DataArray],
    wave_names: Optional[List[str]] = None,
    regions: Optional[List[str]] = None,
    ylims: Optional[List[Optional[Tuple[float, float]]]] = None,
    save_path: str = 'Fig_wave_trend_1979-2020.png'
) -> None:
    """
    绘制多个热带波动在不同区域的年际变化趋势，并保存图像。

    Parameters
    ----------
    wave_filters : list of xarray.DataArray
        每个波动的滤波后 OLR 数据。
    wave_names : list of str, optional
        每个波动的名称，用作图表标题。
    regions : list of str, optional
        要绘制的区域名称，默认为全球和几个关键区域。
    ylims : list of tuple or None, optional
        每个波动图的 y 轴范围。例如 [(7, 10), (5, 9), None]，None 表示自动设置。
    save_path : str
        图像保存的文件路径，默认保存为当前目录下的 'Fig_wave_trend_1979-2020.png'。
    """
    if wave_names is None:
        wave_names = ['Kelvin', 'ER', 'MRG', 'TD']
    if regions is None:
        regions = ['Global', 'Indian-Pacific', 'Indian', 'EastPacific', 'Africa']
    if ylims is None:
        ylims = [None] * len(wave_filters)

    region_bounds = {
        'Global': (0, 360),
        'Indian-Pacific': (130, 250),
        'Indian': (55, 130),
        'EastPacific': (250, 345),
        'Africa': ((345, 360), (0, 55))  # 两段拼接
    }

    plt.rcParams.update({
        'axes.linewidth': 1,
        'font.family': 'Arial',
        'font.size': 12
    })

    n = len(wave_filters)
    ncols = 2
    nrows = (n + 1) // 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 3 * nrows + 1), dpi=300)
    axes = axes.flatten()

    for i, (ds_filter, wave_name) in enumerate(zip(wave_filters, wave_names)):
        ax = axes[i]
        color_cycle = plt.cm.Pastel1.colors
        region_colors = dict(zip(regions, color_cycle))
        all_values = []

        for region in regions:
            if region == 'Africa':
                ds_region = xr.concat([
                    ds_filter.sel(lon=slice(345, 360)),
                    ds_filter.sel(lon=slice(0, 55))
                ], dim='lon')
            else:
                lon_bounds = region_bounds[region]
                ds_region = ds_filter.sel(lon=slice(*lon_bounds))

            ds_region = ds_region.sel(time=slice('1979', '2020'), lat=slice(-15, 15))
            ds_inter_annual = ds_region.groupby('time.year').std('time')
            ds_trend = ds_inter_annual.mean(['lat', 'lon'])

            years = ds_trend['year'].values
            values = ds_trend.values
            all_values.append(values)

            slope, intercept, r, p, stderr = linregress(years, values)
            fit_line = slope * years + intercept

            sig = '**' if p < 0.01 else '*' if p < 0.05 else '^' if p < 0.1 else ''
            color = region_colors[region]
            ax.plot(years, values, color=color, label=f'{region}: {slope:.4f}/yr (p={p:.4f}) {sig}')
            ax.plot(years, fit_line, '--', color=color, linewidth=1)

        ax.set_title(f'{wave_name} Interannual Trend')
        ax.set_xticks(np.arange(1980, 2025, 5))
        ax.grid(True, linestyle='--', linewidth=0.5)

        if i % ncols == 0:
            ax.set_ylabel('Interannual Std Dev')
        if i >= nrows * ncols - ncols:
            ax.set_xlabel('Year')

        if ylims[i] is not None:
            ax.set_ylim(ylims[i])
        else:
            all_flat = np.concatenate(all_values)
            buffer = 0.05 * (all_flat.max() - all_flat.min())
            ax.set_ylim(all_flat.min() - buffer, all_flat.max() + buffer)

        ax.legend(fontsize=10, loc='best', frameon=False)

    # 删除多余子图
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    fig.tight_layout()
    
    plt.show()

