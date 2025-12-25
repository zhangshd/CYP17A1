# 猪CYP17A1血红素参数的MCPB.py工作流程

本文档描述了使用MCPB.py为猪CYP17A1酶与共价结合的血红素生成AMBER力场参数的完整工作流程。

## 概述

**体系**: 猪CYP17A1，Fe(III)-血红素与Cys442配位（重编号后为CYM412）

**方法**: MCPB.py（金属中心参数构建器）+ Seminario方法计算力常数

**量子化学水平**: B3LYP/6-31G*，UltraFine积分格点

**电荷/自旋**: -2 / 6（Fe³⁺高自旋六重态）

## 前置条件

### 软件要求
- AMBER22（含AmberTools）
- Gaussian 16
- Python 3 + NumPy

### 环境配置
```bash
# AMBER环境
source /opt/share/Amber/amber22/amber.sh

# Gaussian环境
export g16root=/opt/share/gaussian
source $g16root/g16/bsd/g16.profile
```

## 目录结构

```
mcpb_heme/
├── README.md                          # 本文档
├── mcpb.in                            # MCPB.py输入文件
├── pig_CYP17A1_mcpb.pdb              # 原始准备好的蛋白PDB
├── FE.mol2                            # Fe离子mol2（输入）
├── HEM.mol2                           # 血红素mol2（输入）
├── HEM.frcmod                         # 血红素frcmod（输入）
├── pig_CYP17A1_small.pdb             # 小模型（Fe + Cys + Heme）
├── pig_CYP17A1_large.pdb             # 大模型（含封端基团）
├── pig_CYP17A1_standard.pdb          # 标准模型
├── pig_CYP17A1_small_opt.com         # Gaussian优化输入
├── pig_CYP17A1_small_fc.com          # Gaussian频率计算输入
├── pig_CYP17A1_large_mk.com          # Gaussian MK电荷输入
├── pig_CYP17A1_small_opt.fchk        # 格式化检查点文件（优化）
├── pig_CYP17A1_large_mk.fchk         # 格式化检查点文件（MK）
├── CM1.mol2                           # 带RESP电荷的CYM残基
├── HM1.mol2                           # 带RESP电荷的HEM残基
├── FE1.mol2                           # 带RESP电荷的Fe离子
├── pig_CYP17A1_mcpbpy.frcmod         # 生成的力场参数
├── pig_CYP17A1_mcpbpy.pdb            # 用于tleap的最终PDB
├── pig_CYP17A1_tleap.in              # tleap输入脚本
├── pig_CYP17A1_dry.prmtop            # 干体系拓扑
├── pig_CYP17A1_dry.inpcrd            # 干体系坐标
├── pig_CYP17A1_solv.prmtop           # 溶剂化体系拓扑
└── pig_CYP17A1_solv.inpcrd           # 溶剂化体系坐标
```

## 分步工作流程

### 步骤0: 准备输入文件

#### 0.1 准备蛋白PDB

蛋白必须满足：
- 配位半胱氨酸为CYM（去质子化）
- Fe作为单独残基
- HEM作为单独残基

```bash
# 示例：从原始结构提取和准备
# CYS442 → CYM（移除HG，更改残基名称）
# Fe从HEM中分离
```

#### 0.2 创建FE.mol2

```bash
cat > FE.mol2 << 'EOF'
@<TRIPOS>MOLECULE
FE
 1 0 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 FE          16.5220    7.9590  -62.7120 FE        1 FE         0.0000
@<TRIPOS>BOND
EOF
```

#### 0.3 生成HEM.mol2

```bash
# 从PDB提取HEM
grep "HEM A 474" pig_CYP17A1_mcpb.pdb > HEM_temp.pdb

# 使用antechamber生成GAFF2原子类型的mol2
antechamber -i HEM_temp.pdb -fi pdb -o HEM.mol2 -fo mol2 -c bcc -at gaff2 -pf y
```

#### 0.4 创建mcpb.in

```bash
cat > mcpb.in << 'EOF'
original_pdb pig_CYP17A1_mcpb.pdb
group_name pig_CYP17A1
cut_off 2.8
ion_ids 7675
ion_mol2files FE.mol2
naa_mol2files HEM.mol2
frcmod_files HEM.frcmod
large_opt 1
force_field ff14SB
water_model OPC
software_version g16
lgmodel_chg -2
lgmodel_spin 6
smmodel_chg -2
smmodel_spin 6
EOF
```

**参数说明**:
- `ion_ids 7675`: PDB中Fe的原子编号
- `cut_off 2.8`: 金属配位的距离截断值
- `lgmodel_chg -2`: 总电荷（HEM²⁻ + 中性蛋白 = -2）
- `lgmodel_spin 6`: 多重度6，对应Fe³⁺高自旋态（S=5/2）

### 步骤1: 生成Gaussian输入文件

```bash
source /opt/share/Amber/amber22/amber.sh
MCPB.py -i mcpb.in -s 1
```

**输出文件**:
- `pig_CYP17A1_small.pdb`（修复后81个原子）
- `pig_CYP17A1_large.pdb`（修复后95个原子）
- `pig_CYP17A1_small_opt.com`
- `pig_CYP17A1_small_fc.com`
- `pig_CYP17A1_large_mk.com`

### ⚠️ 关键修复：移除重复原子

MCPB.py生成的Gaussian输入文件包含重复的氢原子（HHA, HHB, HHC, HHD, HAB, HAC）。运行Gaussian前必须移除这些原子。

#### 修复脚本: `fix_duplicate_atoms.py`

```python
"""
修复MCPB.py生成的Gaussian输入文件中的重复原子。
移除导致段错误的重叠氢原子。
Author: zhangshd
Date: 2024-12-15
"""

import numpy as np

def parse_gaussian_com(filename):
    """解析Gaussian .com文件，返回原子及其坐标"""
    atoms = []
    header_lines = []
    footer_lines = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    in_coords = False
    coord_end = False
    
    for i, line in enumerate(lines):
        stripped = line.strip()
        
        if not in_coords:
            header_lines.append(line)
            parts = stripped.split()
            if len(parts) == 2:
                try:
                    int(parts[0])
                    int(parts[1])
                    in_coords = True
                except ValueError:
                    pass
        elif not coord_end:
            if stripped == '':
                coord_end = True
                footer_lines.append(line)
            else:
                parts = stripped.split()
                if len(parts) >= 4:
                    elem = parts[0]
                    if len(parts) == 5:
                        layer = parts[1]
                        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                        atoms.append({'element': elem, 'layer': layer, 
                                    'x': x, 'y': y, 'z': z, 'original': line})
                    else:
                        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                        atoms.append({'element': elem, 'layer': None,
                                    'x': x, 'y': y, 'z': z, 'original': line})
        else:
            footer_lines.append(line)
    
    return atoms, header_lines, footer_lines

def find_atoms_to_remove(atoms, threshold=0.5):
    """找出需要移除的原子（距离过近时移除后面的原子）"""
    to_remove = set()
    n = len(atoms)
    for i in range(n):
        if i in to_remove:
            continue
        for j in range(i+1, n):
            if j in to_remove:
                continue
            dx = atoms[i]['x'] - atoms[j]['x']
            dy = atoms[i]['y'] - atoms[j]['y']
            dz = atoms[i]['z'] - atoms[j]['z']
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)
            if dist < threshold:
                to_remove.add(j)
                print(f"移除原子 {j+1}: 与原子 {i+1} 的距离 = {dist:.3f} Å")
    return to_remove

def fix_gaussian_file(filename):
    """通过移除重复原子来修复Gaussian输入文件"""
    atoms, header, footer = parse_gaussian_com(filename)
    print(f"原始原子数: {len(atoms)}")
    
    to_remove = find_atoms_to_remove(atoms, threshold=0.5)
    print(f"需要移除的原子数: {len(to_remove)}")
    
    # 备份原始文件
    import shutil
    shutil.copy(filename, filename + '.bak')
    
    # 写入修复后的文件
    with open(filename, 'w') as f:
        for line in header:
            f.write(line)
        for i, atom in enumerate(atoms):
            if i not in to_remove:
                f.write(atom['original'])
        for line in footer:
            f.write(line)
    
    print(f"修复后的文件包含 {len(atoms) - len(to_remove)} 个原子")

# 修复两个文件
fix_gaussian_file('pig_CYP17A1_small_opt.com')
fix_gaussian_file('pig_CYP17A1_large_mk.com')
```

**运行修复**:
```bash
python fix_duplicate_atoms.py
```

预期结果:
- `pig_CYP17A1_small_opt.com`: 87 → 81 个原子
- `pig_CYP17A1_large_mk.com`: 101 → 95 个原子

### 步骤2: 运行Gaussian计算

#### 2.1 几何优化（小模型）

```bash
export g16root=/opt/share/gaussian
source $g16root/g16/bsd/g16.profile

# 如需要，在.com文件中更新资源设置
# %Mem=200GB
# %NProcShared=64

g16 pig_CYP17A1_small_opt.com
```

**运行时间**: 64核约1小时

#### 2.2 频率计算（小模型）

频率计算从优化后的检查点文件读取：

```bash
g16 pig_CYP17A1_small_fc.com
```

**运行时间**: 64核约25分钟

#### 2.3 MK电荷计算（大模型）

```bash
g16 pig_CYP17A1_large_mk.com
```

**运行时间**: 64核约1.5小时

#### 2.4 转换检查点文件

```bash
formchk pig_CYP17A1_small_opt.chk pig_CYP17A1_small_opt.fchk
formchk pig_CYP17A1_large_mk.chk pig_CYP17A1_large_mk.fchk
```

### 步骤3: 修复PDB和Fingerprint文件

从Gaussian输入中移除重复原子后，PDB和fingerprint文件也必须相应更新：

```bash
# 修复 small.pdb（81个原子）
grep -v "HHA \|HHB \|HHC \|HHD \|HAB \|HAC " pig_CYP17A1_small.pdb > temp && mv temp pig_CYP17A1_small.pdb

# 修复 large.pdb（95个原子）
grep -v "HHA \|HHB \|HHC \|HHD \|HAB \|HAC " pig_CYP17A1_large.pdb > temp && mv temp pig_CYP17A1_large.pdb

# 修复 standard.pdb（83个原子）
grep -v "HHA \|HHB \|HHC \|HHD \|HAB \|HAC " pig_CYP17A1_standard.pdb > temp && mv temp pig_CYP17A1_standard.pdb

# 修复 fingerprint 文件
grep -v -E "HHA|HHB|HHC|HHD|HAB|HAC" pig_CYP17A1_large.fingerprint > temp && mv temp pig_CYP17A1_large.fingerprint
grep -v -E "HHA|HHB|HHC|HHD|HAB|HAC" pig_CYP17A1_standard.fingerprint > temp && mv temp pig_CYP17A1_standard.fingerprint
```

### 步骤4: 生成力场参数（Seminario方法）

```bash
source /opt/share/Amber/amber22/amber.sh
MCPB.py -i mcpb.in -s 2
```

**输出**: `pig_CYP17A1_mcpbpy.frcmod`，包含键和角度参数：

```
BOND
Y1-M1   76.1    2.3156      Fe-S(Cys) 键
Y2-M1   51.7    2.0777      Fe-NA 键
Y3-M1   42.6    2.1155      Fe-NB 键
Y4-M1   42.3    2.1151      Fe-NC 键
Y5-M1   52.9    2.0737      Fe-ND 键
```

### 步骤5: RESP电荷拟合

```bash
MCPB.py -i mcpb.in -s 3
```

**输出文件**:
- `CM1.mol2` - 带RESP电荷的CYM残基
- `HM1.mol2` - 带RESP电荷的HEM残基
- `FE1.mol2` - 带RESP电荷的Fe离子

### 步骤6: 生成tleap输入

```bash
MCPB.py -i mcpb.in -s 4
```

**输出**: `pig_CYP17A1_tleap.in`

### 步骤7: 修复用于tleap的最终PDB

生成的 `pig_CYP17A1_mcpbpy.pdb` 可能包含：
1. HEM的重复H原子（HHA, HHB等）
2. HIE残基中错误的HD1原子

```bash
# 移除重复的HEM氢原子
grep -v " HHA \| HHB \| HHC \| HHD \| HAB \| HAC " pig_CYP17A1_mcpbpy.pdb > temp && mv temp pig_CYP17A1_mcpbpy.pdb

# 移除HIE中错误的HD1（HIE应该只有HE2，没有HD1）
grep -v " HD1 HIE\| HD1 CHIE" pig_CYP17A1_mcpbpy.pdb > temp && mv temp pig_CYP17A1_mcpbpy.pdb
```

### 步骤8: 使用tleap构建拓扑

```bash
source /opt/share/Amber/amber22/amber.sh
tleap -f pig_CYP17A1_tleap.in
```

**输出文件**:
- `pig_CYP17A1_dry.prmtop` / `pig_CYP17A1_dry.inpcrd` - 干体系
- `pig_CYP17A1_solv.prmtop` / `pig_CYP17A1_solv.inpcrd` - 溶剂化体系（OPC水模型，离子）

## 生成的力场参数

### 键参数（Seminario方法）

| 键类型 | 力常数 (kcal/mol/Å²) | 平衡键长 (Å) |
|--------|---------------------|--------------|
| Y1-M1 (S-Fe) | 76.1 | 2.316 |
| Y2-M1 (NA-Fe) | 51.7 | 2.078 |
| Y3-M1 (NB-Fe) | 42.6 | 2.116 |
| Y4-M1 (NC-Fe) | 42.3 | 2.115 |
| Y5-M1 (ND-Fe) | 52.9 | 2.074 |

### 角度参数（Seminario方法）

| 角度类型 | 力常数 (kcal/mol/rad²) | 平衡角度 (°) |
|----------|------------------------|--------------|
| CT-Y1-M1 | 61.55 | 106.41 |
| Y1-M1-Y2 | 36.27 | 105.20 |
| Y1-M1-Y3 | 37.47 | 100.52 |
| Y1-M1-Y4 | 39.10 | 102.10 |
| Y1-M1-Y5 | 33.83 | 106.61 |
| Y2-M1-Y3 | 42.59 | 87.50 |
| Y2-M1-Y4 | 45.64 | 152.63 |
| Y2-M1-Y5 | 42.78 | 86.67 |
| Y3-M1-Y4 | 45.12 | 85.69 |
| Y3-M1-Y5 | 44.89 | 152.84 |
| Y4-M1-Y5 | 43.25 | 87.40 |

### 原子类型映射

| 原始类型 | 新类型 | 描述 |
|----------|--------|------|
| SH (CYS) | Y1 | 与Fe配位的硫醇根硫原子 |
| n2 (HEM NA) | Y2 | 吡咯氮A |
| n2 (HEM NB) | Y3 | 吡咯氮B |
| n2 (HEM NC) | Y4 | 吡咯氮C |
| n2 (HEM ND) | Y5 | 吡咯氮D |
| FE | M1 | Fe(III)离子 |

## 最终体系总结

| 属性 | 值 |
|------|-----|
| 总原子数（干体系） | 7,662 |
| 总原子数（溶剂化） | 93,336 |
| 水模型 | OPC |
| 离子参数 | Li/Merz 12-6 OPC |
| 力场 | ff14SB + GAFF2 |
| 净电荷 | +3（用Cl⁻中和） |

## 故障排除

### 1. Gaussian在l202.exe中段错误

**原因**: 输入文件中存在重复/重叠原子

**解决方案**: 使用 `fix_duplicate_atoms.py` 脚本移除重叠的H原子

### 2. MCPB.py "Coordinates not consistent"（坐标不一致）错误

**原因**: PDB和fchk文件之间的原子数不匹配

**解决方案**: 确保在修复Gaussian输入后也更新PDB文件

### 3. tleap "Atom does not have a type"（原子没有类型）错误

**原因**: PDB中的原子在mol2模板中未定义

**解决方案**: 从 `pig_CYP17A1_mcpbpy.pdb` 中移除多余的原子：
- HM1中的HHA, HHB, HHC, HHD, HAB, HAC
- HIE残基中的HD1

### 4. HIE残基HD1问题

**原因**: HIE（ε-质子化组氨酸）应该只有HE2，不应有HD1

**解决方案**: 从最终PDB中移除HIE残基的HD1原子

## 参考文献

1. Li, P.; Merz, K.M. Jr. "MCPB.py: A Python Based Metal Center Parameter Builder." J. Chem. Inf. Model. 2016, 56, 599-604.

2. Seminario, J.M. "Calculation of intramolecular force fields from second-derivative tensors." Int. J. Quantum Chem. 1996, 60, 1271-1277.

3. Case, D.A. et al. "AMBER 22." University of California, San Francisco. 2022.

## 作者信息

- **作者**: zhangshd
- **日期**: 2024-12-15
- **项目**: 猪CYP17A1抑制剂虚拟筛选
