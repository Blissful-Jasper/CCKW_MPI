# 参考文献

## 核心理论文献

### 热带大气波动
1. **Matsuno, T. (1966).** Quasi-geostrophic motions in the equatorial area. *Journal of the Meteorological Society of Japan. Ser. II*, 44(1), 25-43.
   - 首次提出赤道浅水方程波动解，建立了热带波动理论基础

2. **Wheeler, M., & Kiladis, G. N. (1999).** Convectively coupled equatorial waves: Analysis of clouds and temperature in the wavenumber–frequency domain. *Journal of the Atmospheric Sciences*, 56(3), 374-399.
   - 经典的Wheeler-Kiladis波数-频率谱分析方法

3. **Kiladis, G. N., Wheeler, M. C., Haertel, P. T., Straub, K. H., & Roundy, P. E. (2009).** Convectively coupled equatorial waves. *Reviews of Geophysics*, 47(2).
   - 对流耦合赤道波动的综述文章

### Madden-Julian振荡 (MJO)
4. **Madden, R. A., & Julian, P. R. (1971).** Detection of a 40–50 day oscillation in the zonal wind in the tropical Pacific. *Journal of the Atmospheric Sciences*, 28(5), 702-708.
   - MJO的首次发现

5. **Zhang, C. (2005).** Madden‐Julian oscillation. *Reviews of Geophysics*, 43(2).
   - MJO研究综述

### 湿稳定度理论
6. **Raymond, D. J., Sessions, S. L., & Fuchs, Ž. (2009).** A theory for the spinup of tropical depressions. *Quarterly Journal of the Royal Meteorological Society*, 135(643), 1397-1407.
   - Gross Moist Stability (GMS) 理论

7. **Sobel, A. H., & Maloney, E. D. (2013).** Moisture modes and the eastward propagation of the MJO. *Journal of the Atmospheric Sciences*, 70(1), 187-192.
   - 湿模态理论与MJO传播

## 方法学文献

### 时空滤波
8. **Hendon, H. H., & Wheeler, M. C. (2008).** Some space–time spectral analyses of tropical convection and planetary-scale waves. *Journal of the Atmospheric Sciences*, 65(9), 2936-2948.

9. **Roundy, P. E., & Frank, W. M. (2004).** A climatology of waves in the equatorial region. *Journal of the Atmospheric Sciences*, 61(17), 2105-2132.

### 相位分析
10. **Wheeler, M. C., & Hendon, H. H. (2004).** An all-season real-time multivariate MJO index: Development of an index for monitoring and prediction. *Monthly Weather Review*, 132(8), 1917-1932.
    - RMM指数，MJO相位分析标准方法

### 数值方法
11. **Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007).** *Numerical recipes 3rd edition: The art of scientific computing*. Cambridge University Press.

## 数据来源

### 再分析数据
- **ERA5**: Hersbach et al. (2020). The ERA5 global reanalysis. *Quarterly Journal of the Royal Meteorological Society*, 146(730), 1999-2049.
- **MERRA-2**: Gelaro et al. (2017). The Modern-Era Retrospective Analysis for Research and Applications, Version 2 (MERRA-2). *Journal of Climate*, 30(14), 5419-5454.

### 卫星观测
- **NOAA OLR**: Liebmann, B., & Smith, C. A. (1996). Description of a complete (interpolated) outgoing longwave radiation dataset. *Bulletin of the American Meteorological Society*, 77(6), 1275-1277.
- **GPCP降水**: Huffman et al. (2001). Global precipitation at one-degree daily resolution from multisatellite observations. *Journal of Hydrometeorology*, 2(1), 36-50.

### 模式数据
- **CMIP6**: Eyring et al. (2016). Overview of the Coupled Model Intercomparison Project Phase 6 (CMIP6) experimental design and organization. *Geoscientific Model Development*, 9(5), 1937-1958.

## 软件工具参考

### Python生态
- **xarray**: Hoyer, S., & Hamman, J. (2017). xarray: ND labeled arrays and datasets in Python. *Journal of Open Research Software*, 5(1).
- **SciPy**: Virtanen et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17(3), 261-272.
- **Cartopy**: Met Office (2010-2015). Cartopy: a cartographic python library with a matplotlib interface.

### 气象工具
- **NCL (NCAR Command Language)**: Brown et al. (2012). The NCAR Command Language (NCL) (Version 6.0.0).
- **MetPy**: May et al. (2021). MetPy: A Python Package for Meteorological Data.

## 相关代码库

1. **CCEWs diagnostic scripts** (NOAA PSL)
   - https://psl.noaa.gov/data/gridded/data.interp_OLR.html

2. **MJO diagnostics** (NCAR)
   - https://www.ncl.ucar.edu/Applications/mjoclivar.shtml

3. **PyKE (Python Kelvin wave Extraction)**
   - 类似功能的Python工具包

## 推荐阅读

### 教材
- **Holton, J. R., & Hakim, G. J. (2013).** *An introduction to dynamic meteorology* (5th ed.). Academic Press.
- **Vallis, G. K. (2017).** *Atmospheric and oceanic fluid dynamics*. Cambridge University Press.

### 综述文章
- **Zhang, C., Adames, Á. F., Khouider, B., Wang, B., & Yang, D. (2020).** Four theories of the Madden‐Julian oscillation. *Reviews of Geophysics*, 58(3), e2019RG000685.

---

## 引用本工具包

如果您在研究中使用了本工具包，请引用：

```
Jianpu. (2025). Wave Tools: A Python package for tropical atmospheric wave analysis. 
Hohai University. https://github.com/yourusername/wave_tools
```

**BibTeX格式**:
```bibtex
@software{jianpu2025wavetools,
  author = {Jianpu},
  title = {Wave Tools: A Python package for tropical atmospheric wave analysis},
  year = {2025},
  institution = {Hohai University},
  email = {xianpuji@hhu.edu.cn},
  url = {https://github.com/yourusername/wave_tools}
}
```

---

**持续更新中...**
