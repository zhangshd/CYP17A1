# 猪CYP17A1虚拟筛选 - MD模拟与结合自由能分析

## 项目概述

本项目对猪CYP17A1与候选抑制剂的复合物进行分子动力学(MD)模拟，并通过MM/GBSA方法计算结合自由能，用于虚拟筛选饲料添加剂类CYP17A1抑制剂。

### 技术要点
- **软件**: AMBER22 (AmberTools + pmemd.cuda, GPU加速)
- **力场**: ff14SB (蛋白) + GAFF2 (配体) + OPC (水)
- **血红素参数**: MCPB.py生成的Fe-Heme-Cys共价键参数
- **配体约束**: 3原子质心距离约束防止解离（Fe ↔ 最近3个配体原子）
- **模拟时长**: 20 ns生产阶段（虚拟筛选优化）
- **分析方法**: MM/GBSA结合自由能计算

### 最新更新 (2024-12-24)

#### ✅ 距离约束策略升级
- **原方案**: 单原子约束（Fe ↔ C1）→ 配体在加热阶段解离
- **新方案**: 3原子质心约束（Fe ↔ 最近3个重原子的质心）
- **改进**: 自动包含配位原子（如Type II的N22），约束范围收紧至4.0/6.0 Å
- **效果**: 有效防止加热阶段的突发解离（Frame 23, 5.97 Å → 强惩罚能量 ~190 kcal/mol）

#### ✅ 模拟时长优化  
- **50 ns → 20 ns**（虚拟筛选场景）
- 节省60%计算时间（~5h vs 12h per system）
- 节省60%存储空间（2 GB vs 5 GB per system）
- MM/GBSA只用最后5 ns，20 ns已充分平衡

#### ✅ Pipeline简化
- 统一使用GPU配置（systems_gpu/目录）
- 新增批量提交脚本：`batch_submit_md.py`
- 移除临时测试脚本和GPU/CPU分离逻辑

---

## 快速开始

### 完整Pipeline

```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts

# 1. 准备配体参数（GAFF2 + AM1-BCC）
python3 step1_prepare_ligands.py

# 2. 构建蛋白-配体复合物拓扑
python3 step2_build_complex_topology.py

# 3. 生成MD输入文件（含新的3原子质心约束）
python3 step3_setup_md_simulation.py <ligand_id>

# 4. 批量提交GPU任务
python3 batch_submit_md.py

# 5. 监控任务
squeue -u $USER

# 6. MD完成后运行MM/GBSA
python3 step4_run_mmgbsa.py

# 7. 分析结果
python3 step5_analyze_results.py
```

---

## 目录结构

```
md_simulation/
├── scripts/                      # 所有处理脚本
│   ├── config.py                 # 全局配置（GPU版本，20ns）
│   ├── step1_prepare_ligands.py
│   ├── step2_build_complex_topology.py
│   ├── step3_setup_md_simulation.py  # 含3原子质心约束
│   ├── batch_submit_md.py        # 批量提交MD任务（新）
│   ├── step4_run_mmgbsa.py
│   └── step5_analyze_results.py
│
├── mcpb_heme/                    # Fe-Heme MCPB.py参数
│   ├── pig_CYP17A1_dry.pdb       # 输入蛋白结构
│   ├── pig_CYP17A1_mcpbpy.pdb    # 蛋白+Heme结构
│   ├── pig_CYP17A1_mcpbpy.frcmod # 力场参数
│   ├── CM1.mol2, HM1.mol2, FE1.mol2  # 残基模板
│   └── ...
│
├── complexes/
│   ├── ligand_params/            # 配体GAFF2参数
│   │   └── <ligand_id>/
│   │       ├── LIG.mol2, LIG.frcmod, LIG.lib
│   │
│   └── systems_gpu/              # MD模拟系统（GPU版本，主工作目录）
│       └── <ligand_id>/
│           ├── complex_solv.prmtop   # 溶剂化拓扑
│           ├── complex_solv.inpcrd   # 溶剂化坐标
│           ├── dist_restraint.RST    # 3原子质心约束
│           ├── run_md.sh             # SLURM提交脚本
│           ├── min*.in, heat.in, eq*.in, prod.in
│           ├── prod.nc               # 生产轨迹 (~2 GB, 20 ns)
│           └── FINAL_RESULTS_MMGBSA.dat
│
├── results/                      # 汇总结果
│   └── mmgbsa_summary.csv
│
└── README.md                     # 本文档
```

---

## 距离约束详解（防止配体解离）

### 问题背景
初始MD模拟发现配体在**加热阶段（Frame 23, ~20ps, 120K）突然解离**：
- Fe-N22距离从3 Å跳变到6 Å（单帧事件）
- 原因：单原子约束（Fe-C1）未包含配位原子N22
- 配体旋转满足C1约束但N22解离

### 新方案：3原子质心约束

**AMBER DISANG格式**:
```fortran
&rst
  iat=7668, -7691,-7690,-7692,
  r1=0.0, r2=0.0,
  r3=4.0, r4=6.0,
  rk2=50.0, rk3=50.0,
/
```

**参数解释**:
- `7668`: Fe原子
- `-7691,-7690,-7692`: 3个最近配体原子（负号=质心）
- `r3=4.0 Å`: Flat-bottom上限（无惩罚区）
- `r4=6.0 Å`: 线性惩罚上限
- `rk=50.0`: 力常数（kcal/mol/Å²）

**自动选择策略**:
1. 解析`complex_dry.pdb`获取所有配体重原子坐标
2. 计算Fe到各原子距离
3. 选择最近的3个原子（通常包含配位原子）
4. 生成质心约束

**典型示例（阿比特龙AER601）**:
```
Atom  7691 (N22): 2.19 Å  ← Type II配位氮
Atom  7690 (C21): 3.11 Å  ← 吡啶环
Atom  7692 (C23): 3.11 Å  ← 吡啶环
Mean distance: 2.80 Å
```

**能量惩罚**:
- 距离 < 4.0 Å: 0 kcal/mol（自由运动）
- 距离 = 6.0 Å: ~100 kcal/mol
- 距离 = 8.0 Å: ~290 kcal/mol（强力阻止）

---

## 已完成的MM/GBSA结果

| 配体ID | ΔG_bind (kcal/mol) | 标准差 | 标准误 | 来源 |
|--------|-------------------|--------|--------|------|
| ligand_2546 | -69.59 | 3.36 | 0.15 | gras_new |
| ligand_338 | -69.20 | 3.76 | 0.17 | gras_new |
| ligand_1015 | -69.09 | 3.39 | 0.15 | gras_new |
| ligand_2397 | -62.27 | 3.48 | 0.16 | gras_new |
| ligand_325 | -61.55 | 3.76 | 0.17 | gras_new |
| ligand_2010 | -60.39 | 3.29 | 0.33 | gras_new |
| ligand_573 | -54.13 | 4.04 | 0.40 | gras_new |
| ligand_987 | -39.23 | 3.52 | 0.16 | gras_new |
| ligand_519 | -11.84 | 2.74 | 0.12 | gras_new |

### 排名标准

| 评级 | ΔG (kcal/mol) | 说明 |
|------|---------------|------|
| 优秀 | < -50 | 强结合候选物 |
| 良好 | -40 ~ -50 | 值得进一步研究 |
| 一般 | -30 ~ -40 | 需要优化 |
| 较差 | > -30 | 不推荐 |

---

## 第一部分：Fe-Heme参数化 (MCPB.py)

### 1.1 背景说明

CYP17A1是一种含有血红素辅因子的细胞色素P450酶。血红素中的铁(Fe)与蛋白质中的半胱氨酸(Cys)硫原子形成共价配位键。标准的AMBER力场不包含这种金属-蛋白质共价键的参数，因此需要使用MCPB.py工具生成自定义参数。

### 1.2 输入文件准备

```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/mcpb_heme

# 输入文件:
# - pig_CYP17A1_dry.pdb: 去水的蛋白结构，包含:
#   * 蛋白质链 (残基1-473)
#   * HEM血红素 (残基474)
#   * Cys442与Fe的配位
```

### 1.3 MCPB.py参数化流程

**步骤1: 创建MCPB.py输入文件**

```bash
# pig_CYP17A1_mcpbpy.in 内容:
original_pdb pig_CYP17A1_dry.pdb
group_name pig_CYP17A1
cut_off 2.8
ion_ids 7596                    # Fe原子ID
ion_mol2files HEM.mol2          # 血红素mol2文件
naa_mol2files CYS.mol2          # 半胱氨酸mol2文件
frcmod_files hem.frcmod         # 血红素力场修正
large_opt 1                     # 使用大模型优化
software_version g16            # Gaussian 16
```

**步骤2: 运行MCPB.py (分步执行)**

```bash
source /opt/share/Amber/amber22/amber.sh

# Step 1: 生成模型结构
MCPB.py -i pig_CYP17A1_mcpbpy.in -s 1

# Step 2: 生成Gaussian输入文件
MCPB.py -i pig_CYP17A1_mcpbpy.in -s 2

# Step 3: 运行Gaussian计算 (需要高性能计算资源)
# 提交到计算集群运行 pig_CYP17A1_large_mk.com

# Step 4: 生成力场参数
MCPB.py -i pig_CYP17A1_mcpbpy.in -s 4
```

**步骤3: 生成AMBER拓扑**

```bash
# 运行tleap生成最终参数
tleap -f pig_CYP17A1_mcpbpy_tleap.in

# 输出文件:
# - pig_CYP17A1_mcpbpy.prmtop: 蛋白+Heme拓扑文件
# - pig_CYP17A1_mcpbpy.inpcrd: 蛋白+Heme坐标文件
```

### 1.4 关键参数说明

MCPB.py生成的参数包含:
- **M1 (Fe)**: 血红素铁原子
- **Y1-Y4 (卟啉氮)**: 卟啉环上配位的四个氮原子
- **Y5 (Cys-S)**: 半胱氨酸硫原子
- **Fe-N键**: ~2.0 Å
- **Fe-S键**: ~2.3 Å

---

## 第二部分：配体参数化

### 2.1 运行配体参数化脚本

使用antechamber生成GAFF2力场参数和AM1-BCC电荷：

```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts
source /opt/share/Amber/amber22/amber.sh

python3 step1_prepare_ligand_params.py <input.mol2> <ligand_id>
```

### 2.2 脚本功能说明

该脚本对每个配体执行以下操作:

1. **读取配体结构** (从对接结果mol2文件)
2. **添加氢原子并优化几何结构** (如需要)
3. **计算AM1-BCC电荷** (使用antechamber)
4. **生成GAFF2参数** (使用parmchk2)
5. **创建AMBER库文件** (使用tleap)

### 2.3 输入要求

- **输入mol2必须包含氢原子** - 对接结果的mol2可能缺少氢原子，需要从原始数据库提取或使用RDKit添加
- 分子结构完整，无开放价

### 2.4 输出文件

每个配体生成以下文件:
```
ligand_params/ligand_XXX/
├── LIG.mol2      # 带电荷的mol2文件
├── LIG.frcmod    # 力场修正参数
└── LIG.lib       # AMBER库文件
```

---

## 第三部分：复合物拓扑构建

### 3.1 运行拓扑构建脚本

```bash
python3 step2_build_complex_topology.py <ligand.mol2> <ligand_id>
```

### 3.2 关键步骤

1. 修正mol2残基名为"LIG"（tleap兼容）
2. 组合蛋白质PDB和配体坐标
3. 使用tleap加载MCPB.py参数
4. 创建Fe-Heme-Cys键
5. 溶剂化（OPC水模型，12Å缓冲）
6. 添加离子中和电荷

### 3.3 tleap脚本内容

脚本自动生成并执行以下tleap命令:

```bash
# 加载力场
source leaprc.protein.ff14SB    # 蛋白力场
source leaprc.gaff2             # 配体力场
source leaprc.water.opc         # OPC水模型

# 加载MCPB.py生成的参数
loadamberparams pig_CYP17A1_mcpbpy.frcmod
loadoff pig_CYP17A1_mcpbpy.lib

# 加载配体参数
loadamberparams ligand_XXX.frcmod
loadoff ligand_XXX.lib

# 加载蛋白和配体
receptor = loadpdb pig_CYP17A1_mcpbpy.pdb
ligand = loadmol2 ligand_XXX.mol2

# 合并
complex = combine {receptor ligand}

# 溶剂化 (12 Å缓冲区)
solvatebox complex OPCBOX 12.0

# 添加反离子中和体系
addions complex Na+ 0
addions complex Cl- 0

# 保存拓扑
saveamberparm complex complex_solv.prmtop complex_solv.inpcrd
```

### 3.4 输出文件

```
systems/ligand_XXX/
├── complex_solv.prmtop    # 溶剂化复合物拓扑 (~19 MB)
├── complex_solv.inpcrd    # 溶剂化复合物坐标 (~5 MB)
└── tleap.in               # tleap输入文件 (用于复现)
```

### 3.5 注意事项

- **配体mol2必须包含氢原子** - 对接结果的mol2可能缺少氢原子，需要从原始数据库提取或使用RDKit添加
- **PDB格式很重要** - 链ID必须在正确的列位置（第22列）
- **Fe-Heme键** - 必须在tleap中显式创建Fe与卟啉氮、Cys硫原子的共价键
- **GPU内存** - eq2阶段分成10个子阶段以避免GPU网格重建问题

---

## 第四部分：MD模拟设置与运行

### 4.1 运行MD设置脚本

生成输入文件和SLURM脚本：

```bash
python3 step3_setup_md_simulation.py <ligand_id>
```

### 4.2 模拟阶段说明

| 阶段 | 时长 | 系综 | 约束 | 说明 |
|------|------|------|------|------|
| min1 | 5000步 | - | 500 kcal/mol·Å² | 水和离子能量最小化 |
| min2 | 10000步 | - | 无 | 全系统能量最小化 |
| heat | 100 ps | NVT | 10 kcal/mol·Å² | 0→300 K加热 |
| eq1 | 500 ps | NPT | 10 kcal/mol·Å² | 初始密度平衡 |
| eq2 | 500 ps | NPT | 无 | **分10个子阶段运行** |
| prod | 50 ns | NPT | 无 | 生产阶段 |

### 4.3 关键技巧：eq2阶段分段运行

**问题**: 在NPT平衡阶段，由于密度变化较大，GPU版本pmemd.cuda可能报错:
```
Periodic box dimensions have changed too much from their initial values.
```

**解决方案**: 将eq2阶段(500 ps)分成10个子阶段(每个50 ps)运行:

```bash
# eq2_01.in 到 eq2_10.in
# 每个子阶段50 ps (25000步)
# 阶段之间重启，让GPU重建网格
```

### 4.4 SLURM提交脚本

```bash
#!/bin/bash
#SBATCH --job-name=ligand_XXX_md
#SBATCH --partition=G4090
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=72:00:00

source /opt/share/Amber/amber22/amber.sh
module load cuda/11.8
module load fftw/fftw-3.3.8

# 依次运行各阶段
pmemd.cuda -O -i min1.in -p complex_solv.prmtop -c complex_solv.inpcrd -r min1.rst7 -o min1.out -ref complex_solv.inpcrd
pmemd.cuda -O -i min2.in -p complex_solv.prmtop -c min1.rst7 -r min2.rst7 -o min2.out
pmemd.cuda -O -i heat.in -p complex_solv.prmtop -c min2.rst7 -r heat.rst7 -o heat.out -x heat.nc -ref min2.rst7
pmemd.cuda -O -i eq1.in -p complex_solv.prmtop -c heat.rst7 -r eq1.rst7 -o eq1.out -x eq1.nc -ref heat.rst7

# eq2分段运行 (关键!)
for i in 01 02 03 04 05 06 07 08 09 10; do
    if [ "$i" == "01" ]; then
        prev_rst="eq1.rst7"
    else
        prev_i=$(printf "%02d" $((10#$i - 1)))
        prev_rst="eq2_${prev_i}.rst7"
    fi
    pmemd.cuda -O -i eq2_${i}.in -p complex_solv.prmtop -c $prev_rst -r eq2_${i}.rst7 -o eq2_${i}.out -x eq2_${i}.nc
done

# 生产阶段
pmemd.cuda -O -i prod.in -p complex_solv.prmtop -c eq2_10.rst7 -r prod.rst7 -o prod.out -x prod.nc
```

### 4.5 提交MD模拟

```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/complexes/systems/<ligand_id>
sbatch run_md.sh
```

或交互式运行：

```bash
bash run_md.sh
```

### 4.6 批量提交作业

```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/complexes/systems

for sys in ligand_*; do
    cd $sys
    sbatch run_md.sh
    cd ..
done

# 检查作业状态
squeue -u zhangsd
```

### 4.7 输出文件

MD模拟完成后，每个系统目录将包含：

| 文件 | 说明 |
|------|------|
| min1.out, min2.out | 最小化输出 |
| heat.out, heat.nc | 加热输出和轨迹 |
| eq1.out, eq1.nc | 平衡1输出和轨迹 |
| eq2_*.out, eq2_*.nc | 平衡2各子阶段输出和轨迹 |
| prod.out | 生产模拟输出 |
| prod.nc | 生产模拟轨迹 |
| prod.rst7 | 最终重启文件 |

### 4.8 预计运行时间

使用GPU (pmemd.cuda)：
- 能量最小化: ~5分钟
- 加热: ~2分钟
- 平衡: ~10分钟
- 生产(50ns): ~24小时

总计: 约25小时/系统

---

## 第五部分：MM/GBSA结合自由能计算

### 5.1 运行分析设置脚本

```bash
python3 step4_run_mmgbsa.py <ligand_id>
cd /home/zhangsd/repos/CYP17A1/md_simulation/complexes/systems/<ligand_id>
sbatch run_mmgbsa.sh
```

### 5.2 MM/GBSA计算流程

**步骤1: 准备拓扑文件**

使用cpptraj生成三个拓扑:
- `complex_dry.prmtop`: 复合物(去水)
- `receptor.prmtop`: 受体(蛋白+Heme+Fe, 残基1-475)
- `ligand.prmtop`: 配体(残基476)

**步骤2: 去除周期性边界信息**

GB方法不兼容周期性边界条件，需使用parmed去除:
```python
import parmed as pmd
parm = pmd.load_file("complex_dry.prmtop")
parm.box = None
parm.save("complex_dry.prmtop", overwrite=True)
```

**步骤3: 去水轨迹**

```bash
cpptraj -p complex_solv.prmtop << EOF
trajin prod.nc
strip :WAT,Na+,Cl-
trajout prod_dry.nc netcdf
go
EOF
```

**步骤4: 运行MMPBSA.py**

```bash
MMPBSA.py -O \
    -i mmgbsa.in \
    -o FINAL_RESULTS_MMGBSA.dat \
    -cp complex_dry.prmtop \
    -rp receptor.prmtop \
    -lp ligand.prmtop \
    -y prod_dry.nc
```

### 5.3 MM/GBSA参数说明

```
# mmgbsa.in
&general
  startframe=1,
  endframe=500,      # 分析最后500帧 (最后5 ns)
  interval=1,
  verbose=2,
/
&gb
  igb=5,             # GB-OBC2模型 (推荐)
  saltcon=0.15,      # 150 mM盐浓度
/
```

### 5.4 批量提交MM/GBSA作业

```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/complexes/systems

for sys in ligand_*; do
    cd $sys
    sbatch run_mmgbsa.sh
    cd ..
done
```
MM/GBSA时间: 12小时

### 5.5 收集并汇总结果

```bash
# 查看单个结果
grep "DELTA TOTAL" systems/ligand_*/FINAL_RESULTS_MMGBSA.dat

# 汇总并排名
python3 step4_run_mmgbsa.py --collect
```

---

## 第六部分：轨迹可视化

### 6.1 生成可视化文件

```bash
cd systems/ligand_XXX
source /opt/share/Amber/amber22/amber.sh

cpptraj << EOF
parm complex_solv.prmtop
trajin prod.nc 1 last 100    # 每100帧取1帧

strip :WAT,Na+,Cl-           # 去水
autoimage                     # 成像处理
center :1-476 mass origin    # 居中
rms first :1-473@CA          # 对蛋白CA原子叠合

trajout traj_for_pymol.pdb pdb    # PyMOL格式
trajout traj_for_vmd.dcd dcd      # VMD格式
outtraj last_frame.pdb pdb onlyframes 50  # 最后一帧
go
EOF
```

### 6.2 可视化软件使用

**PyMOL:**
```bash
pymol last_frame.pdb
# 或查看轨迹
pymol traj_for_pymol.pdb
```

**VMD:**
```bash
vmd complex_dry.prmtop traj_for_vmd.dcd
```

### 6.3 关键观察点

1. **配体稳定性**: 配体是否稳定在活性位点
2. **Fe-配体距离**: 对于氮杂唑类抑制剂，Fe-N距离应为2.1-2.3 Å
3. **氢键网络**: 配体与关键残基(如Asn202)的氢键
4. **蛋白结构**: 整体结构是否稳定，无明显展开

---

## 第七部分：结果解读

### 7.1 MM/GBSA结果示例

```
DELTA TOTAL = -54.13 +/- 4.04 kcal/mol
```

- **DELTA TOTAL**: 结合自由能 (ΔG)
- **负值越大**: 结合越强
- **标准误差**: 波动大小

### 7.2 能量分解

| 能量项 | 含义 |
|--------|------|
| VDWAALS | 范德华相互作用 |
| EEL | 静电相互作用 |
| EGB | 极性溶剂化能 |
| ESURF | 非极性溶剂化能 |
| DELTA G gas | 气相结合能 |
| DELTA G solv | 溶剂化贡献 |

---

## 第八部分：常见问题与解决方案

### 问题1: Periodic box dimensions have changed too much

**原因**: GPU pmemd.cuda无法在运行中重建pair list网格

**解决**: 将NPT平衡阶段分成多个子阶段运行，每次重启时GPU会重建网格

### 问题2: GB计算报错 "gb>0 is incompatible with periodic boundary conditions"

**原因**: 拓扑文件中包含周期性边界信息(IFBOX标志)

**解决**: 使用parmed去除box信息:
```python
parm.box = None
parm.save(filename, overwrite=True)
```

### 问题3: antechamber计算AM1-BCC电荷失败

**原因**: 分子结构复杂或含有不支持的原子类型

**解决**: 
1. 检查分子结构是否正确
2. 尝试使用sqm替代mopac
3. 手动调整原子类型

### 问题4: MCPB.py Gaussian计算不收敛

**原因**: 初始几何结构不合理或基组选择不当

**解决**:
1. 优化初始结构
2. 使用更大的基组
3. 增加优化步数

### 问题5: tleap报错 "Unknown residue"

**原因**: mol2残基名不符合AMBER命名规范，或未正确加载模板

**解决**: 确保mol2残基名为3字符（如"LIG"），并在tleap中使用相同名称加载模板

### 问题6: OXT原子无类型

**原因**: C端残基的羧基氧(OXT)在某些情况下缺少力场参数

**解决**: 脚本中通过 `remove mol mol.473.OXT` 删除该原子解决

### 问题7: 原子名不匹配

**原因**: PDB文件中的原子名与mol2模板中的不一致

**解决**: 确保PDB中的原子名与mol2模板中的完全一致，必要时手动修正

---

## 附录

### A. 力场配置

- **蛋白质**: ff14SB
- **配体**: GAFF2 + AM1-BCC
- **水模型**: OPC
- **Fe-Heme**: MCPB.py生成的参数

### B. 环境配置

**软件版本:**
- AMBER 22
- CUDA 11.8
- Python 3.9+
- Gaussian 16 (用于MCPB.py)

**环境加载命令:**

```bash
source /opt/share/Amber/amber22/amber.sh
module load cuda/11.8
module load fftw/fftw-3.3.8
```

### C. SLURM集群配置

```bash
#SBATCH --partition=G4090
#SBATCH --gres=gpu:1       # MD模拟需要GPU
#SBATCH --time=72:00:00    # MD时间
#SBATCH --ntasks=8         # MM/GBSA使用MPI并行
```

---

## 参考资料

1. Case, D.A., et al. (2022). AMBER 2022. University of California, San Francisco.
2. Li, P. & Merz, K.M. (2016). MCPB.py: A Python Based Metal Center Parameter Builder. J. Chem. Inf. Model., 56(4), 599-604.
3. Miller, B.R., et al. (2012). MMPBSA.py: An Efficient Program for End-State Free Energy Calculations. J. Chem. Theory Comput., 8(9), 3314-3321.
4. [AMBER Manual](https://ambermd.org/Manuals.html)
5. [MCPB.py Tutorial](https://ambermd.org/tutorials/advanced/tutorial20/MCPBPY.php)
6. [MM/GBSA Tutorial](https://ambermd.org/tutorials/advanced/tutorial3/)

---

**文档版本**: 2.0  
**更新日期**: 2024-12-17  
**作者**: zhangshd
