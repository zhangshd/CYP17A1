# MD Simulation Monitoring Tools

## Quick Start

### Single Check
检查MD模拟当前状态：
```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts
./check_md_progress.sh <ligand_id>
```

### Real-time Monitoring
实时监控（每30秒更新）：
```bash
./check_md_progress.sh <ligand_id> --watch
# or
./check_md_progress.sh <ligand_id> -w
```

按 `Ctrl+C` 退出监控模式。

## Features

### 1. Job Status（任务状态）
- 显示SLURM作业状态
- 作业ID、运行时间、时间限制

### 2. MD Stage（模拟阶段）
- 自动检测当前运行的阶段
- 显示最后完成的阶段
- 检测错误信息

### 3. Fe-Ligand Distance（Fe-配体距离）
- 分析Fe原子到配体质心的距离
- 统计信息：平均值、范围、稳定帧百分比
- 稳定性判断：
  - ✅ STABLE: 平均距离 ≤ 6.0 Å
  - ⚠️  DRIFTING: 6.0 < 平均距离 ≤ 10.0 Å
  - ❌ DISSOCIATED: 平均距离 > 10.0 Å

### 4. Fe-N22 Coordination Bond（Fe-N配位键）
仅对Type II抑制剂（含N22原子）：
- Fe-N键长分析
- 配位稳定性评估：
  - ✅ EXCELLENT: 平均键长 ≤ 2.5 Å
  - ✅ STABLE: 2.5 < 平均键长 ≤ 3.5 Å
  - ⚠️  STRETCHED: 3.5 < 平均键长 ≤ 5.0 Å
  - ❌ BROKEN: 平均键长 > 5.0 Å

### 5. File Sizes（文件大小）
- 显示prod.nc大小
- 系统总大小

## Examples

### 检查AER601系统
```bash
./check_md_progress.sh AER601
```

输出示例：
```
════════════════════════════════════════════════════════════════════
  MD Simulation Progress Monitor
  System: AER601
  Time: 2024-12-25 19:15:53
════════════════════════════════════════════════════════════════════

【1. Job Status】
  Job ID: 213465
  State: RUNNING
  Runtime: 15:31 / UNLIMITED

【2. MD Stage】
  Log: md_213465.log
  Last completed: eq2_10
  Currently running: production MD
  ✅ No errors

【3. Fe-Ligand Distance Analysis】
  Trajectory: prod.nc (Production)
  Total frames: 500
  Average distance: 4.52 Å
  Range: 3.21 - 6.45 Å
  Stable frames (≤6Å): 485 (97.0%)
  Dissociated frames (>10Å): 0 (0.0%)
  ✅ Ligand is STABLE
  ...

【4. Fe-N22 Coordination Bond】
  (Type II inhibitor detected)
  Total frames: 500
  Average bond length: 2.245 Å (ideal: ~2.1 Å)
  Range: 2.089 - 2.512 Å
  Coordinated frames (≤3Å): 500 (100.0%)
  Broken frames (>5Å): 0 (0.0%)
  ✅ Coordination bond is EXCELLENT
  ...

【5. File Sizes】
  prod.nc: 2.4G
  Total system size: 6.1G

════════════════════════════════════════════════════════════════════
  SUMMARY
  ──────────────────────────────────────────────────────────────────
  ✅ System is stable and running normally
  📊 Ligand remains in binding pocket
════════════════════════════════════════════════════════════════════
```

### 实时监控运行中的任务
```bash
./check_md_progress.sh ChMed_1332 --watch
```

屏幕会每30秒自动刷新，实时显示最新状态。

## Troubleshooting

### "No trajectory files found yet"
模拟刚开始，还没有生成轨迹文件。等待几分钟后再检查。

### "Failed to analyze distance"
可能原因：
1. 轨迹文件正在写入（等待当前阶段完成）
2. 拓扑文件损坏
3. cpptraj执行失败（检查AMBER环境）

### "Errors detected"
脚本会自动显示最近的错误信息，常见错误：
- `illegal memory access`: GPU内存错误
- `Periodic box dimensions have changed too much`: NPT平衡问题
- `vlimit exceeded`: 原子速度过大

## Deprecated Scripts

以下旧脚本已被 `check_md_progress.sh` 整合替代：

- ❌ `verify_fe_n_bond.sh` - 功能已整合
- ❌ `quick_check_ligand_stability.sh` - 功能已整合  
- ❌ `monitor_fe_n_bond.sh` - 功能已整合

建议删除这些旧脚本以避免混淆。

## Tips

1. **批量检查多个系统**：
```bash
for sys in GRAS_* ChMed_*; do
    echo "=== Checking $sys ==="
    ./check_md_progress.sh $sys | grep -E "SUMMARY|STABLE|DISSOCIATED"
    echo ""
done
```

2. **只显示有问题的系统**：
```bash
for sys in GRAS_*; do
    result=$(./check_md_progress.sh $sys 2>/dev/null | grep -E "DISSOCIATED|BROKEN|ERROR")
    if [ -n "$result" ]; then
        echo "⚠️  $sys has issues:"
        echo "$result"
    fi
done
```

3. **导出检查结果到文件**：
```bash
./check_md_progress.sh AER601 > AER601_check_$(date +%Y%m%d_%H%M%S).log
```

## Author

zhangshd  
Date: 2024-12-25

