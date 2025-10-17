#!/bin/bash
# 完整清理wave_tools目录，保留最精简的结构

set -e  # 遇到错误立即退出

echo "================================================"
echo "  Wave Tools 完整清理和重构"
echo "================================================"
echo ""

cd /home/m/m301257/wave_tools

# ===== 1. 创建备份 =====
echo "步骤 1/6: 创建备份目录..."
BACKUP_DIR="_backup_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$BACKUP_DIR"
cp *.py *.md *.sh *.txt "$BACKUP_DIR/" 2>/dev/null || true
echo "✓ 备份完成: $BACKUP_DIR/"
echo ""

# ===== 2. 删除所有多余的Markdown文档 =====
echo "步骤 2/6: 删除多余的文档..."
rm -f INDEX.md
rm -f QUICKSTART.md
rm -f CHANGELOG.md
rm -f REFERENCES.md
rm -f ARCHITECTURE.md
rm -f REFACTORING_SUMMARY.md
rm -f REFACTORING_PLAN.md
rm -f README.md  # 删除旧的README
echo "✓ 已删除 8 个多余文档"
echo ""

# ===== 3. 删除重复和旧的Python文件 =====
echo "步骤 3/6: 删除重复和旧的Python文件..."

# 删除已被重命名的原始文件
rm -f ma.py
rm -f Wave_Filter.py
rm -f Wave_Phase.py
rm -f GMS.py
echo "  ✓ 删除旧命名文件 (4个)"

# 删除重复的频谱分析文件
rm -f Spectral_Analysis.py
rm -f spectral_analysis_icon.py
rm -f functions.py
rm -f SpectrumPlotter.py
echo "  ✓ 删除重复的频谱文件 (4个)"

# 删除分散的绘图文件
rm -f core.py
rm -f plot.py
rm -f TaylorDiagram.py
echo "  ✓ 删除分散的绘图文件 (3个)"

# 删除旧的init文件
rm -f __init__.py
echo "  ✓ 删除旧的 __init__.py"

echo "✓ 共删除 12 个Python文件"
echo ""

# ===== 4. 启用新文件 =====
echo "步骤 4/6: 启用新文件..."
mv README_NEW.md README.md
mv __init___NEW.py __init__.py
echo "✓ 已启用新的 README.md 和 __init__.py"
echo ""

# ===== 5. 删除旧的清理脚本 =====
echo "步骤 5/6: 删除旧脚本..."
rm -f cleanup_reorganize.sh
echo "✓ 已删除旧的清理脚本"
echo ""

# ===== 6. 清理缓存 =====
echo "步骤 6/6: 清理Python缓存..."
rm -rf __pycache__
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
find . -type f -name "*.pyc" -delete 2>/dev/null || true
echo "✓ 缓存已清理"
echo ""

# ===== 最终结果 =====
echo "================================================"
echo "  清理完成！"
echo "================================================"
echo ""
echo "📁 最终文件结构:"
echo ""
ls -lh *.py *.md *.txt 2>/dev/null | awk '{printf "  %s  %s\n", $9, $5}'
echo ""
echo "📊 统计:"
echo "  - Python模块: $(ls -1 *.py 2>/dev/null | wc -l) 个"
echo "  - 文档文件: $(ls -1 *.md 2>/dev/null | wc -l) 个"
echo ""
echo "✅ Wave Tools 现在更简洁、清晰、有序！"
echo ""
echo "💡 使用方法:"
echo "   from wave_tools import *"
echo ""
echo "📖 查看完整文档:"
echo "   cat README.md"
echo ""
echo "📦 旧文件备份在: $BACKUP_DIR/"
echo ""
