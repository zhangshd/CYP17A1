# CYP17A1 Virtual Screening Pipeline - 目录结构

```
CYP17A1/
├── data/                           # 原始数据（只读）
│   ├── HEM/                        # 血红素参数文件
│   ├── mols/                       # 配体分子库
│   │   ├── Ch_medicinal_DB_ligprep-out.mol2
│   │   ├── feed_additive_ligprep-out.mol2
│   │   └── gras_new_ligprep-out.mol2
│   └── protein/                    # 蛋白结构
│
├── docking_results/                # 对接结果（只读）
│   ├── Ch_medicinal_DB_ligprep-out/
│   ├── feed_additive_ligprep-out/
│   └── gras_new_ligprep-out/
│
├── md_simulation/
│   ├── scripts/                    # 所有脚本（代码与结果分离）
│   │   ├── config.py               # 全局配置文件
│   │   ├── step1_prepare_ligand_params.py   # 配体参数化
│   │   ├── step2_build_complex_topology.py  # 拓扑构建
│   │   ├── step3_setup_md_simulation.py     # MD输入文件
│   │   ├── step4_run_mmgbsa.py              # MM/GBSA分析
│   │   ├── batch_process.py                  # 批处理主脚本
│   │   ├── slurm_md.sh.template             # MD作业模板
│   │   └── slurm_mmgbsa.sh.template         # MM/GBSA作业模板
│   │
│   ├── mcpb_heme/                  # Fe-Heme MCPB.py参数（已完成）
│   │   ├── pig_CYP17A1_mcpbpy.prmtop
│   │   ├── pig_CYP17A1_mcpbpy.inpcrd
│   │   ├── pig_CYP17A1_mcpbpy.pdb
│   │   └── templates/              # MCPB.py生成的模板
│   │       ├── CM1.mol2
│   │       ├── HM1.mol2
│   │       └── FE1.mol2
│   │
│   ├── ligand_params/              # 配体GAFF2参数
│   │   └── {ligand_id}/
│   │       ├── LIG.mol2            # 带电荷的mol2
│   │       └── LIG.frcmod          # 力场参数
│   │
│   ├── systems/                    # MD模拟系统
│   │   └── {ligand_id}/
│   │       ├── complex_solv.prmtop
│   │       ├── complex_solv.inpcrd
│   │       ├── *.in                # 输入文件
│   │       ├── *.nc                # 轨迹文件
│   │       ├── *.out               # 输出文件
│   │       └── FINAL_RESULTS_MMGBSA.dat
│   │
│   └── results/                    # 汇总结果
│       ├── mmgbsa_summary.csv
│       └── analysis/
│
└── docs/
    ├── 猪CYP17A1虚拟筛选辅基处理.md
    └── MD_Simulation_Protocol_CN.md
```

## 工作流程

### Phase 1: 数据准备
1. 从对接结果中提取配体（需要与原始数据库匹配以获取氢原子）
2. 验证mol2文件完整性

### Phase 2: 配体参数化
1. antechamber生成GAFF2原子类型和AM1-BCC电荷
2. parmchk2生成缺失的力场参数

### Phase 3: 复合物拓扑构建
1. 合并蛋白质PDB和配体坐标
2. 使用tleap加载MCPB.py参数
3. 创建Fe-Heme-Cys共价键
4. 溶剂化和离子中和

### Phase 4: MD模拟
1. 能量最小化 (min1, min2)
2. 加热 (heat)
3. 平衡 (eq1, eq2_01-eq2_10)
4. 生产模拟 (prod, 50ns)

### Phase 5: MM/GBSA分析
1. 去除水分子和离子
2. 生成receptor/ligand拓扑
3. 运行MMPBSA.py

## 命名规范

- 配体ID格式：`{source}_{id}` (例如: ChMed_1332, FeedAdd_44)
- 残基名称：统一使用"LIG"
- 系统目录：与配体ID相同
