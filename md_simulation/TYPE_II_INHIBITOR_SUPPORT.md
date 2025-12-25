# Type II Inhibitor Support for MD Simulation Pipeline

**Author:** zhangshd  
**Date:** 2024-12-26  
**Purpose:** 解决Type II抑制剂（如AER601）在MD模拟中配体漂移的问题

---

## 问题背景

### 症状
- **AER601**（Abiraterone类似物，Type II抑制剂）在MD模拟中配体严重漂移
- 初始结构正确：Fe-N22距离 = 2.19 Å ✓
- eq1后配体漂移：Fe-N22距离 = 11.56 Å ✗
- 导致eq2阶段出现"Periodic box dimensions have changed too much"错误

### 根本原因
1. **AER601是Type II抑制剂**，N22原子应与Fe形成配位键（~2.1 Å）
2. **拓扑中无Fe-N22成键**，因为MCPB.py只为Fe-Heme生成参数
3. **尝试添加共价键失败**：缺少Fe-配体扭转参数，tleap报错

### 对比分析
| 配体 | 类型 | Fe-配体距离变化 | 结果 |
|------|------|----------------|------|
| ligand_2546 | Type I/底物（无氮配位） | 8.99→6.01 Å（靠近） | ✓ 成功 |
| AER601 | Type II（N22配位） | 2.19→11.56 Å（远离） | ✗ 失败 |

---

## 解决方案：距离约束策略

使用**平底势距离约束**（flat-bottom harmonic restraint）维持Fe-N22配位：

### 约束参数
```python
# config.py
LIGAND_FE_COORDINATION = {
    "AER601": {
        "atom": "N22",          # 配位氮原子
        "target_dist": 2.1,     # 目标距离（Å）
        "force_const": 50.0,    # 力常数（kcal/mol/Å²）
        "tolerance": 0.3,       # 平底宽度（±Å）
    },
}

# 分阶段调整力常数
FE_COORDINATION_RESTRAINT = {
    "min1": 50.0,      # 强约束
    "min2": 50.0,
    "heat": 50.0,      # 维持强约束
    "eq1": 30.0,       # 逐渐减弱
    "eq2": 10.0,       # 进一步减弱
    "prod": 0.0,       # 无约束（自然平衡）
}
```

### 约束文件格式（AMBER DISANG）
```
&rst
  iat=475@FE,476@N22,
  r1=0.0, r2=1.80, r3=2.40, r4=999.0,
  rk2=50.0, rk3=50.0,
&end
```

- **平底区间**：1.8-2.4 Å（无惩罚）
- **目标值**：2.1 Å（Type II配位键的典型距离）
- **力常数**：分阶段从50→30→10→0 kcal/mol/Å²

---

## 代码修改

### 1. `config.py`
添加Type II抑制剂配置：
```python
# Line 123-150
LIGAND_FE_COORDINATION = {...}
FE_COORDINATION_RESTRAINT = {...}
```

### 2. `step2_build_complex_topology.py`
- 添加PDB格式输入支持（Line 51-68）
- 添加`read_pdb_coordinates()`函数
- 修改`build_complex_topology()`支持PDB和mol2输入（Line 300-320）
- **不添加Fe-配体共价键**（注释说明使用距离约束）

### 3. `step3_setup_md_simulation.py`
- 添加`create_distance_restraint_file()`函数（Line 30-87）
- 修改所有输入文件生成函数，添加`ligand_id`参数：
  - `create_minimization_inputs()` - 添加nmropt和DISANG
  - `create_heating_input()` - 添加DISANG（与温度ramping共存）
  - `create_equilibration_inputs()` - 添加DISANG（eq1和eq2）
- Production阶段**不使用约束**，让系统自然平衡

---

## 使用方法

### 添加新的Type II抑制剂
在`config.py`中添加配置：
```python
LIGAND_FE_COORDINATION = {
    "AER601": {"atom": "N22", "target_dist": 2.1, "force_const": 50.0, "tolerance": 0.3},
    "新配体ID": {"atom": "配位原子名", "target_dist": 2.1, "force_const": 50.0, "tolerance": 0.3},
}
```

### 运行Pipeline
```bash
# Step 1: 准备配体参数（如果还未准备）
python3 step1_prepare_ligand_params.py 配体.pdb 配体ID

# Step 2: 构建拓扑（支持PDB输入）
python3 step2_build_complex_topology.py 配体.pdb 配体ID --force

# Step 3: 设置MD（自动检测Type II抑制剂并添加约束）
python3 step3_setup_md_simulation.py 配体ID --force

# Step 4: 提交任务
cd ../complexes/systems_gpu/配体ID
sbatch run_md.sh
```

### 监控Fe-配体距离
```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts
bash check_md_progress.sh /path/to/system_dir
```

---

## 预期结果

### Type II抑制剂（带约束）
- ✓ Fe-N22距离维持在1.8-2.4 Å范围
- ✓ 配体不会漂移
- ✓ 顺利通过所有equilibration阶段
- ✓ Production中约束移除，系统自然平衡

### Type I抑制剂/底物（无约束）
- ✓ 保持原有行为（无配位键）
- ✓ 配体可自由移动和调整
- ✓ 如ligand_2546正常运行

---

## 技术细节

### 为什么不用共价键？
尝试在tleap中添加`bond mol.475.FE mol.476.N22`失败：
```
Error! ** No torsion terms for atom types: c2-Y4-M1-nb
```
原因：MCPB.py只为Fe-Heme生成参数，没有Fe-配体的键、角、扭转参数。

### 为什么用平底势？
- **刚性约束**（r2=r3）：过于僵硬，可能导致结构不合理
- **谐波约束**：小偏移也会有惩罚
- **平底势**：在合理范围内（1.8-2.4 Å）无惩罚，超出后才施加力
  - 更接近真实配位键的柔性
  - 允许合理的热涨落

### 为什么逐渐减弱？
```
min/heat: 50 → eq1: 30 → eq2: 10 → prod: 0
```
- **强约束**：初期防止配体漂移
- **逐渐减弱**：让系统逐步适应
- **完全移除**：production中观察真实动力学行为
- 如果配位稳定，移除约束后仍应维持在合理距离

---

## 测试案例

### AER601 (Job 213467)
- **初始**: Fe-N22 = 2.19 Å
- **约束**: 50/30/10 kcal/mol/Å² (min/eq1/eq2)
- **预期**: 维持在2.1±0.3 Å，模拟正常完成
- **状态**: 已提交，监控中

---

## 未来扩展

1. **自动识别Type II抑制剂**：
   - 扫描配体结构中的N/O原子
   - 计算与Fe的距离
   - 自动添加配位约束

2. **更精确的约束参数**：
   - 根据量化计算优化目标距离
   - 调整力常数和平底范围

3. **支持多配位**：
   - 某些配体可能有多个配位原子
   - 扩展为列表格式

---

## 参考文献

1. Type II CYP17A1 inhibitors coordination chemistry
2. AMBER distance restraint (DISANG) documentation
3. MCPB.py for metalloprotein parametrization

