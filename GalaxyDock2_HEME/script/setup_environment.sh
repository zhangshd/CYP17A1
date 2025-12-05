#!/bin/bash
"""
setup_environment.sh
GalaxyDock2_HEME 环境配置脚本
自动配置 conda 环境和依赖

Author: zhangshd
Date: 2024-12-05

Usage:
    bash setup_environment.sh
"""

set -e  # 遇到错误立即退出

echo "=============================================="
echo "GalaxyDock2_HEME 环境配置脚本"
echo "=============================================="

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 配置变量
CONDA_ENV_NAME="dock"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# 函数：打印状态
print_status() {
    echo -e "${GREEN}[✓]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[!]${NC} $1"
}

print_error() {
    echo -e "${RED}[✗]${NC} $1"
}

# 1. 检查 Chimera
echo ""
echo "1. 检查 UCSF Chimera..."
if command -v chimera &> /dev/null; then
    CHIMERA_VERSION=$(chimera --version 2>&1 | head -1)
    print_status "Chimera 已安装: $CHIMERA_VERSION"
    CHIMERA_PATH=$(which chimera)
else
    print_error "Chimera 未安装!"
    echo "   请从 https://www.cgl.ucsf.edu/chimera/download.html 下载安装"
    echo "   安装后将 chimera/bin 添加到 PATH"
    exit 1
fi

# 2. 检查/创建 Conda 环境
echo ""
echo "2. 配置 Conda 环境 ($CONDA_ENV_NAME)..."

# 初始化 conda
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "/opt/share/miniconda3/etc/profile.d/conda.sh" ]; then
    source "/opt/share/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
fi

# 检查环境是否存在
if conda env list | grep -q "^${CONDA_ENV_NAME} "; then
    print_status "Conda 环境 '$CONDA_ENV_NAME' 已存在"
else
    print_warning "创建 Conda 环境 '$CONDA_ENV_NAME'..."
    conda create -n $CONDA_ENV_NAME python=3.12 -y
    print_status "Conda 环境创建完成"
fi

# 激活环境
conda activate $CONDA_ENV_NAME

# 3. 安装依赖
echo ""
echo "3. 安装 Python 依赖..."

# 检查并安装依赖
PACKAGES="rdkit openbabel numpy"
for pkg in $PACKAGES; do
    if conda list | grep -q "^${pkg} "; then
        print_status "$pkg 已安装"
    else
        print_warning "安装 $pkg..."
        conda install $pkg -c conda-forge -y
    fi
done

# 4. 获取路径信息
echo ""
echo "4. 获取环境路径..."
CONDA_PREFIX=$(conda info --base)/envs/$CONDA_ENV_NAME
OBABEL_PATH="$CONDA_PREFIX/bin/obabel"
PYTHON_PATH="$CONDA_PREFIX/bin/python"

print_status "Python: $PYTHON_PATH"
print_status "OpenBabel: $OBABEL_PATH"

# 5. 更新 run_GalaxyDock2_heme.py 中的路径
echo ""
echo "5. 更新脚本配置..."
RUN_SCRIPT="$SCRIPT_DIR/run_GalaxyDock2_heme.py"

if [ -f "$RUN_SCRIPT" ]; then
    # 备份原文件
    cp "$RUN_SCRIPT" "${RUN_SCRIPT}.bak"
    
    # 替换 OBABEL_PATH
    sed -i "s|OBABEL_PATH = Path('.*')|OBABEL_PATH = Path('$OBABEL_PATH')|g" "$RUN_SCRIPT"
    print_status "已更新 OBABEL_PATH: $OBABEL_PATH"
else
    print_error "未找到 run_GalaxyDock2_heme.py"
fi

# 6. 验证安装
echo ""
echo "6. 验证安装..."
echo "=============================================="

$PYTHON_PATH --version
$OBABEL_PATH -V 2>&1 | head -1
$PYTHON_PATH -c "import rdkit; print(f'RDKit: {rdkit.__version__}')"
$PYTHON_PATH -c "import numpy; print(f'NumPy: {numpy.__version__}')"

echo "=============================================="
print_status "环境配置完成!"
echo ""
echo "使用方法:"
echo "  conda activate $CONDA_ENV_NAME"
echo "  export PATH=\"$CONDA_PREFIX/bin:\$PATH\""
echo "  python $SCRIPT_DIR/run_GalaxyDock2_heme.py --help"
echo ""
