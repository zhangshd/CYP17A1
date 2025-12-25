# Changelog - 2024-12-25

## 🎯 主要更新

### 1. Pipeline脚本修改（复现成功案例）

根据 `systems/ligand_*` 目录下的成功MD模拟案例，修改了pipeline脚本：

#### `step3_setup_md_simulation.py` 主要变更

| 项目 | 修改前 (systems_gpu) | 修改后 (systems/ligand_*) |
|------|---------------------|---------------------------|
| **距离约束** | ✅ dist_restraint.RST (3原子质心) | ❌ 无距离约束 |
| **加热阶段** | 3阶段 (heat1→heat2→heat3) | **单阶段** (heat: 0→300K, 50ps) |
| **NVT预平衡** | ✅ nvt_long (5ns) | ❌ 移除 |
| **eq2约束** | 逐渐减弱 (2.0→1.5→...→0) | **全部无约束** (ntr=0) |
| **最小化步数** | 20000 | **10000** |
| **生产时长** | 20 ns | **50 ns** (用户改回20ns) |

**原因**：systems/ligand_* 下的模拟已成功完成且MMGBSA结果良好，采用其参数可提高成功率。

#### `config.py` 修改

```python
MD_PROTOCOL = {
    "prod": {
        "description": "Production NPT (20 ns)",  # 用户设置
        "nstlim": 10000000,  # 20 ns
    }
}
```

#### `step4_run_mmgbsa.py` 修改

- 添加 `--force` 参数支持
- 使用 parmed 移除周期性边界条件（匹配成功案例）
- 生成 prod_dry.nc 作为 MMGBSA 输入
- SLURM配置：ntasks=32（用户调整）

### 2. 监控脚本整合 ⭐

整合了三个功能重叠的监控脚本：

#### 旧脚本（已弃用）

- ❌ `verify_fe_n_bond.sh` - Fe-N键验证
- ❌ `quick_check_ligand_stability.sh` - 配体稳定性快速检查
- ❌ `monitor_fe_n_bond.sh` - 实时Fe-N键监控

#### 新工具

1. **`check_md_progress.sh`** - 统一监控脚本
   ```bash
   # 单次检查
   ./check_md_progress.sh <ligand_id>
   
   # 实时监控（30秒更新）
   ./check_md_progress.sh <ligand_id> --watch
   ```

2. **`batch_check_status.sh`** - 批量检查脚本
   ```bash
   # 检查所有系统
   ./batch_check_status.sh
   
   # 只显示有问题的
   ./batch_check_status.sh --problems-only
   ```

#### 新功能特性

| 功能 | 说明 |
|------|------|
| 🔍 自动阶段检测 | 识别min/heat/eq/prod当前阶段 |
| 📏 Fe-配体距离 | 分析Fe到配体质心距离（通用） |
| 🔗 Fe-N22键长 | Type II抑制剂配位键分析 |
| ✅ 稳定性评估 | STABLE / DRIFTING / DISSOCIATED |
| ❌ 错误检测 | 自动识别AMBER错误 |
| 📊 统计信息 | 平均值、范围、稳定帧百分比 |
| 🔄 实时监控 | --watch模式支持 |
| 📂 轨迹自适应 | 自动选择最新轨迹文件 |

### 3. 文档更新

新增文档：

- `MONITORING_TOOLS.md` - 监控工具详细说明
- `README_MONITORING.md` - 监控工具快速参考
- `CHANGELOG_20241225.md` - 本变更日志

## 📝 使用示例

### Pipeline使用

```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts

# 1. 为新配体生成MD输入文件
python3 step3_setup_md_simulation.py <ligand_id>

# 2. 提交MD任务
cd ../complexes/systems_gpu/<ligand_id>
sbatch run_md.sh

# 3. MD完成后生成MMGBSA文件
cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts
python3 step4_run_mmgbsa.py <ligand_id>

# 4. 提交MMGBSA任务
cd ../complexes/systems_gpu/<ligand_id>
sbatch run_mmgbsa.sh
```

### 监控使用

```bash
# 检查单个系统
./check_md_progress.sh AER601

# 实时监控
./check_md_progress.sh AER601 --watch

# 批量检查
./batch_check_status.sh

# 只看有问题的
./batch_check_status.sh --problems-only
```

## 🔧 技术细节

### MD输入文件变更

**min1.in/min2.in**：
- 移除 `nmropt=1` 和 `DISANG` 行
- maxcyc: 20000 → 10000

**heat.in**：
- 单阶段加热（0→300K，50ps）
- 移除 nvt_long.in

**eq1.in**：
- 100ps，带backbone约束

**eq2_01.in ~ eq2_10.in**：
- 每段50ps
- 全部 ntr=0（无约束）

**prod.in**：
- 50 ns (用户改为20ns)
- 无距离约束

### MMGBSA流程变更

```bash
# Step 1: 准备拓扑（cpptraj）
cpptraj -i prep_topologies.cpptraj

# Step 1b: 移除box信息（parmed）
python3 << 'EOF'
import parmed as pmd
for prmtop in ["complex_dry.prmtop", "receptor.prmtop", "ligand.prmtop"]:
    parm = pmd.load_file(prmtop)
    parm.box = None
    parm.save(prmtop, overwrite=True)
EOF

# Step 2: Strip轨迹
cpptraj -i strip_traj.cpptraj

# Step 3: 运行MMPBSA.py
MMPBSA.py -i mmgbsa.in -cp complex_dry.prmtop ... -y prod_dry.nc
```

## 🎉 成功案例对比

### systems/ligand_* (参考标准)

- ✅ 无距离约束
- ✅ 简单的heating protocol
- ✅ MM/GBSA成功完成
- ✅ 结合能合理（-69 ~ -12 kcal/mol）

### 当前AER601问题

- ❌ Fe-N键断裂（28 Å）
- ❌ Periodic box错误
- 📊 正在重新测试简化protocol

## 📋 待办事项

- [ ] 用新protocol重新提交失败的系统
- [ ] 验证新pipeline在3-5个系统上的效果
- [ ] 如果成功，批量处理所有系统
- [ ] 可选：清理旧的监控脚本到 `_deprecated/`

## 🔗 相关文件

**Pipeline脚本**：
- `scripts/step3_setup_md_simulation.py`
- `scripts/step4_run_mmgbsa.py`
- `scripts/config.py`

**监控脚本**：
- `scripts/check_md_progress.sh` ⭐ 新
- `scripts/batch_check_status.sh` ⭐ 新

**文档**：
- `scripts/MONITORING_TOOLS.md`
- `scripts/README_MONITORING.md`
- `README.md` (主文档)

## 👤 Author

zhangshd  
Date: 2024-12-25

## 📞 Support

如有问题，请参考：
1. `README.md` - 主文档
2. `MONITORING_TOOLS.md` - 监控工具说明
3. 检查 `systems/ligand_*` 目录下的成功案例

