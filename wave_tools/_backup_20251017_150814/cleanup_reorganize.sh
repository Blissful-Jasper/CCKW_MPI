#!/bin/bash
# 清理和重组wave_tools目录

echo "开始清理和重组 wave_tools..."

cd /home/m/m301257/wave_tools

# ===== 1. 备份旧文件 =====
echo "1. 创建备份目录..."
mkdir -p _backup_old_files
cp *.py *.md _backup_old_files/ 2>/dev/null

# ===== 2. 删除多余的markdown文档 =====
echo "2. 删除多余的文档..."
rm -f INDEX.md QUICKSTART.md CHANGELOG.md REFERENCES.md ARCHITECTURE.md REFACTORING_SUMMARY.md

# ===== 3. 删除重复的Python文件 =====
echo "3. 删除重复和旧的Python文件..."
# 删除原始文件（已被重命名的新文件替代）
rm -f ma.py Wave_Filter.py Wave_Phase.py GMS.py

# 删除重复的频谱分析文件
rm -f functions.py Spectral_Analysis.py SpectrumPlotter.py spectral_analysis_icon.py

# 删除分散的绘图文件
rm -f core.py plot.py TaylorDiagram.py

# ===== 4. 重命名新文件 =====
echo "4. 启用新文件..."
mv README_NEW.md README.md
mv __init___NEW.py __init__.py

# ===== 5. 保留的文件列表 =====
echo "5. 最终文件结构:"
echo ""
echo "📁 wave_tools/"
echo "├── 📄 README.md           (新的统一文档)"
echo "├── 📄 requirements.txt"
echo "├── 📄 __init__.py         (新的导入结构)"
echo "├── 🐍 matsuno.py          (Matsuno理论)"
echo "├── 🐍 spectral.py         (频谱分析 - 新)"
echo "├── 🐍 filters.py          (波动滤波)"
echo "├── 🐍 phase.py            (相位分析)"
echo "├── 🐍 diagnostics.py      (诊断工具)"
echo "├── 🐍 plotting.py         (所有绘图 - 新)"
echo "└── 🐍 utils.py            (工具函数)"
echo ""

echo "✅ 清理完成！"
echo ""
echo "旧文件已备份至: _backup_old_files/"
echo ""
echo "现在可以使用:"
echo "  from wave_tools import *"
