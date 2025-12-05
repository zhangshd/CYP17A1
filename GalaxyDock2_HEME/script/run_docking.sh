#!/bin/bash
"""
run_docking.sh
GalaxyDock2_HEME 分子对接自动化脚本
用于批量或单次分子对接

Author: zhangshd
Date: 2024-12-05

Usage:
    # 单个配体对接
    bash run_docking.sh -p protein.pdb -l ligand.mol2 -o output_dir
    
    # 批量对接（配体目录）
    bash run_docking.sh -p protein.pdb -L ligands_dir/ -o output_dir
    
    # 指定血红素残基号
    bash run_docking.sh -p protein.pdb -l ligand.mol2 -o output_dir --heme 600 --chain A
"""

set -e

# ============================================
# 配置区域 - 根据实际情况修改
# ============================================
CONDA_ENV="dock"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GALAXYDOCK_HOME="$(dirname "$SCRIPT_DIR")"

# Conda 路径（根据实际安装位置修改）
CONDA_BASE="/opt/share/miniconda3"
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    CONDA_BASE="$HOME/miniconda3"
fi

# ============================================
# 颜色和辅助函数
# ============================================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_header() {
    echo ""
    echo -e "${BLUE}============================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}============================================${NC}"
}

print_status() {
    echo -e "${GREEN}[✓]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[!]${NC} $1"
}

print_error() {
    echo -e "${RED}[✗]${NC} $1"
}

# ============================================
# 参数解析
# ============================================
usage() {
    echo "GalaxyDock2_HEME 分子对接脚本"
    echo ""
    echo "用法:"
    echo "  $0 -p <protein.pdb> -l <ligand.mol2> -o <output_dir> [options]"
    echo "  $0 -p <protein.pdb> -L <ligands_dir> -o <output_dir> [options]"
    echo ""
    echo "必需参数:"
    echo "  -p, --protein    蛋白PDB文件（含血红素）"
    echo "  -l, --ligand     单个配体MOL2文件"
    echo "  -L, --ligand-dir 配体目录（批量对接）"
    echo "  -o, --output     输出目录"
    echo ""
    echo "可选参数:"
    echo "  --heme           血红素残基号 (默认: 自动检测)"
    echo "  --chain          蛋白链ID (默认: 第一条链)"
    echo "  -x, -y, -z       对接盒中心坐标 (默认: 配体中心)"
    echo "  --seed           随机种子 (可重复性)"
    echo "  -h, --help       显示帮助"
    echo ""
    echo "示例:"
    echo "  $0 -p protein.pdb -l ligand.mol2 -o ./output --heme 600 --chain A"
    echo "  $0 -p protein.pdb -L ./ligands/ -o ./batch_output"
    exit 1
}

# 默认值
PROTEIN=""
LIGAND=""
LIGAND_DIR=""
OUTPUT_DIR=""
HEME_RES=""
CHAIN=""
CENTER_X=""
CENTER_Y=""
CENTER_Z=""
RANDOM_SEED=""

# 解析参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--protein)
            PROTEIN="$2"
            shift 2
            ;;
        -l|--ligand)
            LIGAND="$2"
            shift 2
            ;;
        -L|--ligand-dir)
            LIGAND_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --heme)
            HEME_RES="$2"
            shift 2
            ;;
        --chain)
            CHAIN="$2"
            shift 2
            ;;
        -x)
            CENTER_X="$2"
            shift 2
            ;;
        -y)
            CENTER_Y="$2"
            shift 2
            ;;
        -z)
            CENTER_Z="$2"
            shift 2
            ;;
        --seed)
            RANDOM_SEED="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            print_error "未知参数: $1"
            usage
            ;;
    esac
done

# 检查必需参数
if [ -z "$PROTEIN" ]; then
    print_error "缺少蛋白文件 (-p)"
    usage
fi

if [ -z "$LIGAND" ] && [ -z "$LIGAND_DIR" ]; then
    print_error "缺少配体文件 (-l) 或配体目录 (-L)"
    usage
fi

if [ -z "$OUTPUT_DIR" ]; then
    print_error "缺少输出目录 (-o)"
    usage
fi

# ============================================
# 环境初始化
# ============================================
print_header "GalaxyDock2_HEME 分子对接"

echo "初始化环境..."

# 初始化 conda
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate $CONDA_ENV

# 确保使用 conda 环境的 Python（而非 Chimera 的 Python）
export PATH="$CONDA_BASE/envs/$CONDA_ENV/bin:$PATH"

PYTHON_BIN="$CONDA_BASE/envs/$CONDA_ENV/bin/python"
RUN_SCRIPT="$SCRIPT_DIR/run_GalaxyDock2_heme.py"

print_status "Python: $($PYTHON_BIN --version)"
print_status "GalaxyDock HOME: $GALAXYDOCK_HOME"

# 检查文件
if [ ! -f "$PROTEIN" ]; then
    print_error "蛋白文件不存在: $PROTEIN"
    exit 1
fi

if [ ! -f "$RUN_SCRIPT" ]; then
    print_error "对接脚本不存在: $RUN_SCRIPT"
    exit 1
fi

# ============================================
# 对接函数
# ============================================
run_single_docking() {
    local lig_file=$1
    local out_dir=$2
    local lig_name=$(basename "$lig_file" .mol2)
    
    echo ""
    print_status "对接配体: $lig_name"
    
    # 构建命令
    CMD="$PYTHON_BIN $RUN_SCRIPT \
        -d $GALAXYDOCK_HOME \
        -p $PROTEIN \
        -l $lig_file \
        --out_dir $out_dir"
    
    # 添加可选参数
    [ -n "$HEME_RES" ] && CMD="$CMD --heme_res_num $HEME_RES"
    [ -n "$CHAIN" ] && CMD="$CMD --chain $CHAIN"
    [ -n "$CENTER_X" ] && CMD="$CMD -x $CENTER_X"
    [ -n "$CENTER_Y" ] && CMD="$CMD -y $CENTER_Y"
    [ -n "$CENTER_Z" ] && CMD="$CMD -z $CENTER_Z"
    [ -n "$RANDOM_SEED" ] && CMD="$CMD --random_seed $RANDOM_SEED"
    
    # 执行对接
    echo "执行命令: $CMD"
    eval $CMD
    
    # 检查结果
    if [ -f "$out_dir/GD2_HEME_fb.mol2" ]; then
        print_status "对接完成: $lig_name"
        # 显示 Top 3 结果
        echo "Top 3 构象:"
        head -6 "$out_dir/GD2_HEME_fb.E.info" | tail -3
    else
        print_error "对接失败: $lig_name"
    fi
}

# ============================================
# 执行对接
# ============================================
mkdir -p "$OUTPUT_DIR"

if [ -n "$LIGAND" ]; then
    # 单个配体对接
    if [ ! -f "$LIGAND" ]; then
        print_error "配体文件不存在: $LIGAND"
        exit 1
    fi
    
    run_single_docking "$LIGAND" "$OUTPUT_DIR"
    
elif [ -n "$LIGAND_DIR" ]; then
    # 批量对接
    if [ ! -d "$LIGAND_DIR" ]; then
        print_error "配体目录不存在: $LIGAND_DIR"
        exit 1
    fi
    
    # 获取所有 mol2 文件
    LIGAND_FILES=($(find "$LIGAND_DIR" -name "*.mol2" -type f))
    TOTAL=${#LIGAND_FILES[@]}
    
    if [ $TOTAL -eq 0 ]; then
        print_error "未找到 mol2 文件在: $LIGAND_DIR"
        exit 1
    fi
    
    print_status "找到 $TOTAL 个配体文件"
    
    # 逐个对接
    COUNT=0
    for lig_file in "${LIGAND_FILES[@]}"; do
        COUNT=$((COUNT + 1))
        lig_name=$(basename "$lig_file" .mol2)
        lig_output="$OUTPUT_DIR/$lig_name"
        
        print_header "[$COUNT/$TOTAL] 对接 $lig_name"
        
        mkdir -p "$lig_output"
        run_single_docking "$lig_file" "$lig_output"
    done
    
    # 汇总结果
    print_header "对接汇总"
    echo "结果目录: $OUTPUT_DIR"
    echo ""
    echo "配体名称              总能量      HM_E(血红素)"
    echo "----------------------------------------"
    
    for lig_file in "${LIGAND_FILES[@]}"; do
        lig_name=$(basename "$lig_file" .mol2)
        lig_output="$OUTPUT_DIR/$lig_name"
        
        if [ -f "$lig_output/GD2_HEME_fb.E.info" ]; then
            # 提取 Top 1 的能量
            TOP1=$(sed -n '4p' "$lig_output/GD2_HEME_fb.E.info")
            ENERGY=$(echo $TOP1 | awk '{print $2}')
            HM_E=$(echo $TOP1 | awk '{print $7}')
            printf "%-20s %10s %10s\n" "$lig_name" "$ENERGY" "$HM_E"
        else
            printf "%-20s %10s %10s\n" "$lig_name" "FAILED" "-"
        fi
    done
fi

print_header "完成"
echo "输出目录: $OUTPUT_DIR"
echo ""
