#!/bin/bash
"""
run_all_databases.sh
运行所有三个分子库的批量对接

Author: zhangshd
Date: 2024-12-05

Usage:
    bash run_all_databases.sh [n_workers]
    
Example:
    bash run_all_databases.sh 16
"""

set -e

# ============================================
# 配置区域
# ============================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

# 蛋白文件
PROTEIN="$PROJECT_ROOT/data/protein/pig_CYP_merged.pdb"

# 配体库目录
LIGAND_DIR="$PROJECT_ROOT/data/mols"

# 输出根目录
OUTPUT_ROOT="$PROJECT_ROOT/docking_results"

# 血红素配置
HEME_RES_NUM="600"
CHAIN="A"

# 并行进程数（默认8，可通过命令行参数覆盖）
N_WORKERS=${1:-8}

# Conda配置
CONDA_ENV="dock"
CONDA_BASE="/opt/share/miniconda3"

# ============================================
# 初始化
# ============================================
echo "=============================================="
echo "GalaxyDock2_HEME 批量对接 - 全数据库"
echo "=============================================="
echo "蛋白: $PROTEIN"
echo "配体库目录: $LIGAND_DIR"
echo "输出目录: $OUTPUT_ROOT"
echo "并行进程数: $N_WORKERS"
echo "=============================================="

# 初始化conda
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate $CONDA_ENV
export PATH="$CONDA_BASE/envs/$CONDA_ENV/bin:$PATH"

# 创建输出目录
mkdir -p "$OUTPUT_ROOT"

# Python脚本路径
BATCH_SCRIPT="$SCRIPT_DIR/batch_docking.py"

# ============================================
# 运行对接
# ============================================

# 获取所有mol2文件
DATABASES=($(ls "$LIGAND_DIR"/*.mol2 2>/dev/null))

if [ ${#DATABASES[@]} -eq 0 ]; then
    echo "错误：未找到mol2文件在 $LIGAND_DIR"
    exit 1
fi

echo ""
echo "找到 ${#DATABASES[@]} 个数据库:"
for db in "${DATABASES[@]}"; do
    mol_count=$(grep -c "@<TRIPOS>MOLECULE" "$db")
    echo "  - $(basename "$db"): $mol_count 个分子"
done
echo ""

# 记录开始时间
START_TIME=$(date +%s)

# 逐个处理数据库
for db in "${DATABASES[@]}"; do
    db_name=$(basename "$db" .mol2)
    db_output="$OUTPUT_ROOT/$db_name"
    
    echo ""
    echo "############################################################"
    echo "处理数据库: $db_name"
    echo "############################################################"
    
    python "$BATCH_SCRIPT" \
        --protein "$PROTEIN" \
        --ligand_db "$db" \
        --output "$db_output" \
        --heme_res_num "$HEME_RES_NUM" \
        --chain "$CHAIN" \
        --n_workers "$N_WORKERS" \
        --conda_env "$CONDA_ENV" \
        --conda_base "$CONDA_BASE"
    
    echo ""
done

# 计算总耗时
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))
SECONDS=$((ELAPSED % 60))

echo ""
echo "=============================================="
echo "全部对接完成!"
echo "=============================================="
echo "总耗时: ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo ""
echo "结果目录: $OUTPUT_ROOT"
echo ""
echo "汇总文件:"
for db in "${DATABASES[@]}"; do
    db_name=$(basename "$db" .mol2)
    echo "  - $OUTPUT_ROOT/$db_name/${db_name}_summary.csv"
done
echo "=============================================="
