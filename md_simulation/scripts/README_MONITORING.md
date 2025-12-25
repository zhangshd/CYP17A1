# MD Simulation Monitoring - Quick Reference

## 🚀 整合完成

三个功能重叠的脚本已整合为一个统一的监控工具：

| 旧脚本 | 状态 | 说明 |
|--------|------|------|
| `verify_fe_n_bond.sh` | ✅ 已整合 | Fe-N键验证功能已集成 |
| `quick_check_ligand_stability.sh` | ✅ 已整合 | 配体稳定性检查已集成 |
| `monitor_fe_n_bond.sh` | ✅ 已整合 | 实时监控功能已集成 |

**新工具**：
- ✨ `check_md_progress.sh` - 统一的MD监控脚本
- 📊 `batch_check_status.sh` - 批量状态检查

## 📖 使用方法

### 1. 检查单个系统

```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts

# 单次检查
./check_md_progress.sh AER601

# 实时监控（30秒更新）
./check_md_progress.sh AER601 --watch
```

### 2. 批量检查多个系统

```bash
# 检查所有系统
./batch_check_status.sh

# 检查特定模式
./batch_check_status.sh "GRAS_*"

# 只显示有问题的系统
./batch_check_status.sh --problems-only
```

## 📊 功能对比

### check_md_progress.sh 特性

| 功能 | 说明 |
|------|------|
| **自动阶段检测** | 识别当前运行的MD阶段（min/heat/eq/prod） |
| **Fe-配体距离** | 分析Fe到配体质心距离（所有配体类型） |
| **Fe-N22键长** | 针对Type II抑制剂的配位键分析 |
| **稳定性评估** | STABLE / DRIFTING / DISSOCIATED |
| **错误检测** | 自动识别AMBER错误信息 |
| **实时监控** | --watch模式支持连续监控 |
| **轨迹自适应** | 自动选择最新可用轨迹 |

### 输出示例

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
  Stable frames (≤6Å): 485 (97.0%)
  ✅ Ligand is STABLE

【4. Fe-N22 Coordination Bond】
  Total frames: 500
  Average bond length: 2.245 Å (ideal: ~2.1 Å)
  ✅ Coordination bond is EXCELLENT

【5. File Sizes】
  prod.nc: 2.4G
  Total system size: 6.1G

════════════════════════════════════════════════════════════════════
  SUMMARY
  ──────────────────────────────────────────────────────────────════
  ✅ System is stable and running normally
  📊 Ligand remains in binding pocket
════════════════════════════════════════════════════════════════════
```

## 🔧 高级用法

### 并行检查多个系统

```bash
# 快速扫描所有GRAS系统
for sys in GRAS_*; do
    echo "=== $sys ==="
    ./check_md_progress.sh $sys | grep -E "STABLE|DISSOCIATED|ERROR"
done
```

### 导出检查报告

```bash
# 生成时间戳报告
./check_md_progress.sh AER601 > reports/AER601_$(date +%Y%m%d_%H%M%S).log
```

### 自动化监控

```bash
# 每小时检查一次，记录到日志
while true; do
    ./batch_check_status.sh --problems-only >> daily_check.log
    echo "---" >> daily_check.log
    sleep 3600
done
```

## 🗑️ 清理旧脚本

旧的监控脚本现在可以安全删除：

```bash
cd /home/zhangsd/repos/CYP17A1/md_simulation/scripts

# 可选：备份旧脚本
mkdir -p _deprecated
mv verify_fe_n_bond.sh _deprecated/
mv quick_check_ligand_stability.sh _deprecated/
mv monitor_fe_n_bond.sh _deprecated/
```

## 📝 注意事项

1. **脚本依赖**：
   - AMBER22 环境
   - cpptraj工具
   - squeue命令（SLURM作业管理）

2. **距离判断标准**：
   - Fe-配体质心 ≤ 6Å：稳定
   - Fe-N22键 ≤ 3Å：配位稳定
   - 标准来自文献和实验数据

3. **实时监控注意**：
   - 使用`--watch`模式时按Ctrl+C退出
   - 建议在tmux/screen会话中运行长时间监控

## 📚 更多信息

详细说明请参阅：`MONITORING_TOOLS.md`

## 作者

zhangshd  
2024-12-25

