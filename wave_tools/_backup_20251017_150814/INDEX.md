# Wave Tools 文档导航

> **欢迎使用 Wave Tools!** 👋  
> 本页面提供所有文档的快速导航。

---

## 🚀 快速开始

**第一次使用？** 从这里开始：

1. 📖 [README.md](README.md) - 了解项目概况和核心功能
2. ⚡ [QUICKSTART.md](QUICKSTART.md) - 5分钟快速入门
3. 🎯 运行第一个示例 - 见下方"新手教程"

---

## 📚 文档目录

### 主要文档

| 文档 | 内容 | 适合人群 | 预计阅读时间 |
|------|------|----------|------------|
| [📘 README.md](README.md) | 项目完整说明书<br/>• 功能介绍<br/>• 安装指南<br/>• API文档<br/>• 使用示例 | 所有用户 | 30分钟 |
| [⚡ QUICKSTART.md](QUICKSTART.md) | 快速开始指南<br/>• 5分钟入门<br/>• 示例代码<br/>• FAQ | 新手用户 | 10分钟 |
| [🏗️ ARCHITECTURE.md](ARCHITECTURE.md) | 架构设计文档<br/>• 模块结构<br/>• 数据流程<br/>• 设计原则 | 开发者 | 15分钟 |
| [📋 CHANGELOG.md](CHANGELOG.md) | 版本更新历史<br/>• 新功能<br/>• 改进项<br/>• 已知问题 | 所有用户 | 5分钟 |
| [📚 REFERENCES.md](REFERENCES.md) | 参考文献列表<br/>• 理论基础<br/>• 方法学<br/>• 数据源 | 研究人员 | 20分钟 |
| [✅ REFACTORING_SUMMARY.md](REFACTORING_SUMMARY.md) | 重构总结<br/>• 改进内容<br/>• 新增文档<br/>• 下一步计划 | 开发者 | 10分钟 |

### 技术文件

| 文件 | 说明 |
|------|------|
| [requirements.txt](requirements.txt) | Python依赖包列表 |
| [__init__.py](__init__.py) | 包初始化文件（含导出API） |

---

## 🎓 学习路径

### 路径1️⃣: 新手入门 (1-2小时)

```
步骤1: 阅读简介
  └─→ README.md (简介和模块结构章节)

步骤2: 运行第一个例子
  └─→ QUICKSTART.md (例子1: Wheeler-Kiladis频谱)

步骤3: 理解架构
  └─→ ARCHITECTURE.md (模块架构图)

步骤4: 尝试其他例子
  └─→ QUICKSTART.md (例子2-3)
```

### 路径2️⃣: 进阶使用 (1周)

```
步骤1: 深入功能
  └─→ README.md (核心功能章节，详细阅读)

步骤2: 学习理论
  └─→ REFERENCES.md (核心理论文献)

步骤3: 实践应用
  └─→ 使用自己的数据运行完整分析流程

步骤4: 参数调优
  └─→ 修改配置参数，理解不同设置的影响
```

### 路径3️⃣: 高级开发 (持续)

```
步骤1: 阅读源码
  └─→ 各.py文件（带详细注释）

步骤2: 扩展功能
  └─→ 添加自定义波动类型或分析方法

步骤3: 贡献代码
  └─→ 提交改进或新功能

步骤4: 维护文档
  └─→ 更新CHANGELOG和README
```

---

## 🔍 按需查找

### 我想了解...

**...如何安装**  
→ [README.md § 安装要求](README.md#💻-安装要求)

**...有哪些功能**  
→ [README.md § 核心功能](README.md#🔧-核心功能)

**...如何计算频谱**  
→ [QUICKSTART.md § 例子1](QUICKSTART.md#例子-1-计算wheeler-kiladis频谱)  
→ [README.md § Wheeler-Kiladis频谱分析](README.md#2-wheeler-kiladis频谱分析-spectral_analysispy)

**...如何提取波动**  
→ [QUICKSTART.md § 例子2](QUICKSTART.md#例子-2-提取kelvin波)  
→ [README.md § 波动滤波与提取](README.md#3-波动滤波与提取-wave_filterpy)

**...如何绘制频谱图**  
→ [QUICKSTART.md § 例子3](QUICKSTART.md#例子-3-绘制带matsuno模态的频谱)  
→ [README.md § 频谱可视化](README.md#5-频谱可视化-spectrumplotterpy)

**...理论基础**  
→ [REFERENCES.md § 核心理论文献](REFERENCES.md#核心理论文献)

**...设计思路**  
→ [ARCHITECTURE.md](ARCHITECTURE.md)

**...版本更新**  
→ [CHANGELOG.md](CHANGELOG.md)

**...数据格式要求**  
→ [QUICKSTART.md § 数据格式要求](QUICKSTART.md#3-数据格式要求)

**...常见问题**  
→ [QUICKSTART.md § 常见问题](QUICKSTART.md#4-常见问题)

---

## 📦 模块速查

### 核心分析模块

| 模块 | 功能 | 文档链接 |
|------|------|----------|
| `ma.py` | Matsuno模态计算 | [README § Matsuno模态](README.md#1-matsuno模态与色散关系-mapy) |
| `Spectral_Analysis.py` | WK频谱分析 | [README § 频谱分析](README.md#2-wheeler-kiladis频谱分析-spectral_analysispy) |
| `Wave_Filter.py` | 波动滤波 | [README § 波动滤波](README.md#3-波动滤波与提取-wave_filterpy) |

### 诊断工具

| 模块 | 功能 | 文档链接 |
|------|------|----------|
| `Wave_Phase.py` | 相位检测 | [README § 波动相位](README.md#4-波动相位与峰值检测-wave_phasepy) |
| `GMS.py` | 湿稳定度 | [README § GMS](README.md#9-总体湿稳定度分析-gmspy) |
| `core.py` | CCKW包络 | [README § core.py](README.md#核心功能) |

### 可视化工具

| 模块 | 功能 | 文档链接 |
|------|------|----------|
| `SpectrumPlotter.py` | 频谱绘图 | [README § 频谱可视化](README.md#5-频谱可视化-spectrumplotterpy) |
| `plot.py` | 地图绘图 | [README § 地图绘图](README.md#7-地图绘图工具-plotpy) |
| `TaylorDiagram.py` | 泰勒图 | [README § 泰勒图](README.md#8-泰勒图-taylordiagrampy) |

### 工具函数

| 模块 | 功能 | 文档链接 |
|------|------|----------|
| `utils.py` | 通用工具 | [README § 通用工具](README.md#6-通用工具函数-utilspy) |
| `functions.py` | 频谱辅助 | [README § functions.py](README.md#核心功能) |

---

## 🛠️ 实用工具

### 代码示例索引

所有示例代码位于：
- [README.md § 使用示例](README.md#📚-使用示例) - 4个详细示例
- [QUICKSTART.md § 五分钟入门](QUICKSTART.md#2-五分钟入门) - 3个快速示例

### API文档索引

| 类/函数 | 说明 | 位置 |
|---------|------|------|
| `SpectralAnalysis` | 频谱分析类 | [README § API](README.md#spectralanalysis-类) |
| `WaveFilter` | 波动滤波类 | [README § API](README.md#wavefilter-类) |
| `SpectrumPlotter` | 频谱绘图类 | [README § API](README.md#spectrumplotter-类) |
| `TaylorDiagram` | 泰勒图类 | [README § API](README.md#8-泰勒图-taylordiagrampy) |

---

## 📞 获取帮助

### 遇到问题？

1. 🔍 先查看 [QUICKSTART.md § 常见问题](QUICKSTART.md#4-常见问题)
2. 📖 查阅相关模块的详细文档
3. 📧 联系作者：xianpuji@hhu.edu.cn
4. 🐛 提交Issue（如果使用GitHub）

### 贡献与反馈

- 发现Bug → 提交Issue或联系作者
- 功能建议 → 通过Email讨论
- 文档改进 → 欢迎提交修正
- 代码贡献 → 查看 [REFACTORING_SUMMARY.md § 下一步建议](REFACTORING_SUMMARY.md#🚀-下一步建议)

---

## 📄 文件清单

### 文档文件 (7个)
```
├── INDEX.md                    (本文件)
├── README.md                   (15 KB) 主文档
├── QUICKSTART.md               (3.5 KB) 快速开始
├── ARCHITECTURE.md             (5.9 KB) 架构设计
├── CHANGELOG.md                (3.4 KB) 更新日志
├── REFERENCES.md               (5.2 KB) 参考文献
└── REFACTORING_SUMMARY.md      (6.8 KB) 重构总结
```

### Python模块 (12个)
```
核心分析:
├── ma.py                       (19 KB)  Matsuno模态
├── Spectral_Analysis.py        (17 KB)  频谱分析
├── spectral_analysis_icon.py   (33 KB)  非结构网格频谱
└── Wave_Filter.py              (17 KB)  波动滤波

诊断工具:
├── Wave_Phase.py               (9.7 KB) 相位检测
├── GMS.py                      (3.5 KB) 湿稳定度
└── core.py                     (2.4 KB) CCKW工具

可视化:
├── SpectrumPlotter.py          (8.8 KB) 频谱绘图
├── plot.py                     (7.5 KB) 地图绘图
└── TaylorDiagram.py            (6.4 KB) 泰勒图

基础设施:
├── utils.py                    (23 KB)  工具函数
└── functions.py                (7.4 KB) 辅助函数
```

### 配置文件 (2个)
```
├── requirements.txt            (299 B)  依赖清单
└── __init__.py                 (1.2 KB) 包初始化
```

---

## 🎯 下一步做什么？

### 如果你是新用户
1. 📖 阅读 [README.md](README.md) 的简介部分
2. ⚡ 跟随 [QUICKSTART.md](QUICKSTART.md) 运行第一个示例
3. 🎓 选择上面的"学习路径"继续深入

### 如果你是开发者
1. 🏗️ 阅读 [ARCHITECTURE.md](ARCHITECTURE.md) 了解设计
2. 📋 查看 [CHANGELOG.md](CHANGELOG.md) 了解开发历史
3. ✅ 参考 [REFACTORING_SUMMARY.md](REFACTORING_SUMMARY.md) 了解改进方向

### 如果你是研究人员
1. 📚 阅读 [REFERENCES.md](REFERENCES.md) 了解理论基础
2. 🔬 运行完整的分析流程
3. 📊 使用工具包发表研究成果

---

## 📊 统计信息

- **总代码量**: ~4,400 行 Python代码
- **总文档量**: ~1,160+ 行 Markdown文档
- **模块数量**: 12 个 Python模块
- **文档数量**: 7 个说明文档
- **示例代码**: 7+ 个完整示例
- **参考文献**: 15+ 篇经典论文

---

**最后更新**: 2025年10月17日  
**版本**: 1.0.0  
**维护者**: Jianpu (xianpuji@hhu.edu.cn)

---

**祝使用愉快！有任何问题随时联系 📧**
