# GalaxyDock2_HEME 使用说明

## 概述

GalaxyDock2_HEME 是专门针对**含血红素(Heme)蛋白**的分子对接程序，能够处理血红素辅基与配体之间的特殊相互作用（如Fe配位），这是传统对接程序难以处理的。

---

## 1. 环境配置

### 1.1 系统要求

- **操作系统**: Linux (Ubuntu 22.04 推荐)
- **Python**: >= 3.6
- **UCSF Chimera**: 1.18 或更高版本
- **Conda**: 用于管理Python环境

### 1.2 安装 UCSF Chimera

Chimera 需要手动从官网下载（需接受非商业许可协议）：

1. 访问 https://www.cgl.ucsf.edu/chimera/download.html
2. 下载 Linux 64-bit 版本 (chimera-1.x-linux_x86_64.bin)
3. 安装：

```bash
chmod +x ~/Downloads/chimera-*-linux_x86_64.bin
~/Downloads/chimera-*-linux_x86_64.bin
# 安装路径建议: /home/<username>/chimera
```

4. 添加到 PATH（在 `~/.bashrc` 中添加）：

```bash
export PATH="/home/<username>/chimera/bin:$PATH"
```

5. 验证安装：

```bash
source ~/.bashrc
chimera --version
```

### 1.3 创建 Conda 环境

```bash
# 创建名为 dock 的环境
conda create -n dock python=3.12 -y

# 激活环境
conda activate dock

# 安装依赖
conda install rdkit openbabel numpy -c conda-forge -y
```

### 1.4 配置 GalaxyDock2_HEME

修改脚本中的 OpenBabel 路径（`GalaxyDock2_HEME/script/run_GalaxyDock2_heme.py`）：

```python
# 第14行，修改为你的实际路径
OBABEL_PATH = Path('/opt/share/miniconda3/envs/dock/bin/obabel')
```

---

## 2. 输入文件准备

### 2.1 蛋白结构 (PDB格式)

蛋白PDB文件需要包含：
- 蛋白链（ATOM记录）
- 血红素辅基（HETATM记录，残基名为 HEM 或 HEC）
- 共晶配体（可选，用于定义对接中心）

### 2.2 配体结构 (MOL2格式)

配体需要为 **mol2 格式**，包含：
- 3D坐标
- 正确的原子类型
- 部分电荷（可选，脚本会自动计算）

如果配体是其他格式，可用 OpenBabel 转换：

```bash
obabel ligand.sdf -O ligand.mol2 --gen3d
```

---

## 3. 运行对接

### 3.1 基本用法

```bash
conda activate dock

# 确保使用正确的Python（避免Chimera的Python覆盖）
export PATH="/opt/share/miniconda3/envs/dock/bin:$PATH"

python GalaxyDock2_HEME/script/run_GalaxyDock2_heme.py \
    -d ./GalaxyDock2_HEME/ \
    -p <蛋白PDB文件> \
    -l <配体MOL2文件> \
    --heme_res_num <血红素残基号> \
    --chain <链ID> \
    --out_dir <输出目录>
```

### 3.2 参数说明

| 参数 | 必需 | 说明 | 示例 |
|------|------|------|------|
| `-d` | 是 | GalaxyDock2_HEME 主目录 | `./GalaxyDock2_HEME/` |
| `-p` | 是 | 蛋白PDB文件（含HEM） | `protein.pdb` |
| `-l` | 是 | 配体MOL2文件 | `ligand.mol2` |
| `--heme_res_num` | 否 | 血红素残基编号 | `600` |
| `--chain` | 否 | 蛋白链ID | `A` |
| `--out_dir` | 否 | 输出目录 | `./output/` |
| `-x, -y, -z` | 否 | 对接盒中心坐标 | 默认使用配体中心 |
| `-size_x/y/z` | 否 | 对接盒尺寸(Å) | 默认 22.5 |
| `--random_seed` | 否 | 随机种子（可重复性） | `12345` |

### 3.3 示例命令

```bash
# 对接 pig_CYP17A1 与 MYT 配体
python GalaxyDock2_HEME/script/run_GalaxyDock2_heme.py \
    -d ./GalaxyDock2_HEME/ \
    -p data/protein/pig_CYP_merged.pdb \
    -l data/mols/MYT_sup.mol2 \
    --heme_res_num 600 \
    --chain A \
    --out_dir ./docking_output
```

---

## 4. 输出文件说明

| 文件 | 说明 |
|------|------|
| `GD2_HEME_fb.mol2` | **最终结果**：按能量排序的对接构象（默认100个） |
| `GD2_HEME_fb.E.info` | 各构象的能量打分详情 |
| `GD2_HEME_cl.mol2` | 聚类后的代表构象 |
| `contact.pdb` | 预处理后的蛋白+血红素结构 |
| `HEM.mol2` | 处理后的血红素mol2文件（含部分电荷） |
| `input_ligand.mol2` | 预处理后的配体 |
| `gd2_heme.in` | 对接参数配置文件 |
| `gd2_heme.log` | 运行日志 |

### 4.1 能量打分列说明 (E.info)

| 列名 | 说明 |
|------|------|
| Energy | 总打分能量 (越负越好) |
| ATDK_E | AutoDock 能量项 |
| INT_E | 内部能量 |
| DS_E | DrugScore 能量 |
| **HM_E** | **血红素相互作用能** (Fe配位) |
| PLP | PLP 打分 |

---

## 5. 批量对接（虚拟筛选）

### 5.1 概述

`batch_docking.py` 脚本支持对多分子mol2库进行并行对接，适用于大规模虚拟筛选。

**功能特点**：
- 自动拆分多分子mol2文件
- 多进程并行对接（基于subprocess）
- 自定义最大进程数
- 自动汇总结果到CSV和mol2文件
- 仅保留必要的输出文件（节省磁盘空间）

### 5.2 基本用法

```bash
conda activate dock
export PATH="/opt/share/miniconda3/envs/dock/bin:$PATH"

python GalaxyDock2_HEME/script/batch_docking.py \
    --protein data/protein/pig_CYP_merged.pdb \
    --ligand_db data/mols/library.mol2 \
    --output ./docking_results \
    --heme_res_num 600 \
    --chain A \
    --n_workers 8
```

### 5.3 参数说明

| 参数 | 必需 | 说明 | 默认值 |
|------|------|------|--------|
| `-p, --protein` | 是 | 蛋白PDB文件 | - |
| `-l, --ligand_db` | 是 | 配体库mol2文件（可含多分子） | - |
| `-o, --output` | 是 | 输出目录 | - |
| `-d, --galaxydock_home` | 否 | GalaxyDock2_HEME主目录 | 自动检测 |
| `--heme_res_num` | 否 | 血红素残基号 | None |
| `--chain` | 否 | 蛋白链ID | None |
| `-n, --n_workers` | 否 | 并行进程数 | 4 |
| `--conda_env` | 否 | Conda环境名 | dock |
| `--conda_base` | 否 | Conda安装路径 | 自动检测 |
| `--keep_temp` | 否 | 保留临时拆分文件 | False |

### 5.4 输出结构

```
docking_results/
├── library_name_summary.csv      # 汇总CSV（按能量排序）
├── library_name_top1_poses.mol2  # 所有分子的Top1构象
├── failed_molecules.txt          # 失败的分子列表（如有）
├── molecule_1/                   # 每个分子的子目录
│   ├── box.pdb
│   ├── contact.pdb
│   ├── GD2_HEME_fb.E.info
│   └── GD2_HEME_fb.mol2
├── molecule_2/
│   └── ...
└── ...
```

### 5.5 汇总CSV格式

| 列名 | 说明 |
|------|------|
| Rank | 排名（按能量） |
| Molecule | 分子名称 |
| Energy | 总能量 |
| ATDK_E | AutoDock能量 |
| INT_E | 内部能量 |
| DS_E | DrugScore能量 |
| HM_E | 血红素相互作用能 |
| PLP | PLP打分 |

### 5.6 运行所有数据库

如果有多个配体库文件，可使用wrapper脚本一次性运行：

```bash
# 运行所有数据库，使用16个并行进程
bash GalaxyDock2_HEME/script/run_all_databases.sh 16
```

---

## 6. 结果分析

### 6.1 查看能量排名

```bash
head -20 docking_output/GD2_HEME_fb.E.info
```

### 6.2 可视化对接结果

使用 PyMOL 查看：

```bash
pymol docking_output/contact.pdb docking_output/GD2_HEME_fb.mol2
```

在 PyMOL 中：
```
# 显示蛋白表面
show surface, contact
# 显示配体棍状模型
show sticks, GD2_HEME_fb
# 按构象编号选择
select pose1, GD2_HEME_fb and state 1
```

---

## 7. 常见问题

### Q1: 日志中出现 "ERROR: failed in getting ang parm"

这是关于 N-Fe-N 键角参数的警告，**不影响对接结果**。血红素作为刚性辅基处理，其内部几何不会改变。

### Q2: 日志中出现 "WARNING: unknown ligand atom name"

氢原子命名与内部数据库不完全匹配，**对结果影响很小**。

### Q3: Python 版本错误 (SyntaxError: f-string)

Chimera 自带 Python 2.7 可能覆盖了系统 Python。解决方法：

```bash
export PATH="/opt/share/miniconda3/envs/dock/bin:$PATH"
```

### Q4: PyMOL 中铁原子不显示

原始PDB文件中元素符号可能是 `Fe2+`（4字符），PyMOL 需要标准2字符格式。可在 PyMOL 中手动修复：

```
alter (name FE), elem='Fe'
```

---

## 8. 参考文献

- GalaxyDock2-HEME: [论文链接]
- UCSF Chimera: https://www.cgl.ucsf.edu/chimera/

---

## 作者

zhangshd  
日期: 2024-12-05
