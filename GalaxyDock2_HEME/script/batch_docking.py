"""
batch_docking.py
GalaxyDock2_HEME 批量对接脚本

支持多分子mol2库文件的并行对接，使用subprocess进行多进程并行计算。

Features:
- 自动拆分多分子mol2文件
- 多进程并行对接
- 自动汇总结果到CSV和mol2文件
- 仅保留必要的输出文件

Author: zhangshd
Date: 2024-12-05

Usage:
    python batch_docking.py \
        --protein data/protein/pig_CYP_merged.pdb \
        --ligand_db data/mols/library.mol2 \
        --output output_dir \
        --heme_res_num 600 \
        --chain A \
        --n_workers 8
"""

import argparse
import csv
import os
import shutil
import subprocess as sp
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple


@dataclass
class DockingResult:
    """对接结果数据类"""
    mol_name: str
    energy: float
    atdk_e: float
    int_e: float
    ds_e: float
    hm_e: float
    plp: float
    success: bool
    mol2_content: str = ""


def split_mol2_file(mol2_file: Path, output_dir: Path) -> List[Tuple[str, Path]]:
    """
    拆分多分子mol2文件为单个分子文件
    
    Args:
        mol2_file: 输入的多分子mol2文件
        output_dir: 输出目录
        
    Returns:
        List of (mol_name, mol2_path) tuples
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    molecules = []
    current_mol_lines = []
    current_mol_name = None
    mol_count = 0
    
    with open(mol2_file, 'r') as f:
        for line in f:
            if line.startswith('@<TRIPOS>MOLECULE'):
                # 保存上一个分子
                if current_mol_lines and current_mol_name:
                    mol_path = output_dir / f"{current_mol_name}.mol2"
                    with open(mol_path, 'w') as out_f:
                        out_f.writelines(current_mol_lines)
                    molecules.append((current_mol_name, mol_path))
                
                # 开始新分子
                current_mol_lines = [line]
                mol_count += 1
                current_mol_name = None
            elif current_mol_lines:
                current_mol_lines.append(line)
                # 获取分子名（第二行）
                if len(current_mol_lines) == 2:
                    raw_name = line.strip()
                    # 清理分子名，移除非法字符
                    current_mol_name = raw_name.replace(':', '_').replace('/', '_').replace(' ', '_')
                    if not current_mol_name:
                        current_mol_name = f"mol_{mol_count}"
    
    # 保存最后一个分子
    if current_mol_lines and current_mol_name:
        mol_path = output_dir / f"{current_mol_name}.mol2"
        with open(mol_path, 'w') as out_f:
            out_f.writelines(current_mol_lines)
        molecules.append((current_mol_name, mol_path))
    
    return molecules


def run_single_docking(
    mol_name: str,
    mol2_path: Path,
    protein_path: Path,
    galaxydock_home: Path,
    output_dir: Path,
    heme_res_num: Optional[str],
    chain: Optional[str],
    conda_env: str,
    conda_base: Path
) -> DockingResult:
    """
    运行单个分子的对接
    
    Args:
        mol_name: 分子名称
        mol2_path: 配体mol2文件路径
        protein_path: 蛋白PDB文件路径
        galaxydock_home: GalaxyDock2_HEME主目录
        output_dir: 输出目录
        heme_res_num: 血红素残基号
        chain: 链ID
        conda_env: conda环境名
        conda_base: conda安装路径
        
    Returns:
        DockingResult对象
    """
    mol_output_dir = output_dir / mol_name
    mol_output_dir.mkdir(parents=True, exist_ok=True)
    
    # 构建命令
    python_bin = conda_base / "envs" / conda_env / "bin" / "python"
    run_script = galaxydock_home / "script" / "run_GalaxyDock2_heme.py"
    
    cmd = [
        str(python_bin),
        str(run_script),
        "-d", str(galaxydock_home),
        "-p", str(protein_path),
        "-l", str(mol2_path),
        "--out_dir", str(mol_output_dir)
    ]
    
    if heme_res_num:
        cmd.extend(["--heme_res_num", heme_res_num])
    if chain:
        cmd.extend(["--chain", chain])
    
    # 设置环境变量，确保使用正确的Python
    env = os.environ.copy()
    env["PATH"] = f"{conda_base}/envs/{conda_env}/bin:" + env.get("PATH", "")
    
    # 运行对接
    try:
        result = sp.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,  # 10分钟超时
            env=env
        )
        
        # 检查结果文件
        e_info_file = mol_output_dir / "GD2_HEME_fb.E.info"
        mol2_result_file = mol_output_dir / "GD2_HEME_fb.mol2"
        
        if e_info_file.exists() and mol2_result_file.exists():
            # 解析能量信息（取Top 1）
            with open(e_info_file, 'r') as f:
                lines = f.readlines()
                # 跳过前3行header，取第4行（Top 1）
                if len(lines) >= 4:
                    parts = lines[3].split()
                    if len(parts) >= 7:
                        # 读取Top 1的mol2构象
                        mol2_content = read_first_mol2_entry(mol2_result_file)
                        
                        # 清理不需要的文件，只保留指定的4个
                        cleanup_output_dir(mol_output_dir)
                        
                        return DockingResult(
                            mol_name=mol_name,
                            energy=float(parts[1]),
                            atdk_e=float(parts[2]),
                            int_e=float(parts[3]),
                            ds_e=float(parts[4]),
                            hm_e=float(parts[5]),
                            plp=float(parts[6]),
                            success=True,
                            mol2_content=mol2_content
                        )
        
        # 对接失败
        return DockingResult(
            mol_name=mol_name,
            energy=0.0,
            atdk_e=0.0,
            int_e=0.0,
            ds_e=0.0,
            hm_e=0.0,
            plp=0.0,
            success=False
        )
        
    except sp.TimeoutExpired:
        return DockingResult(
            mol_name=mol_name,
            energy=0.0,
            atdk_e=0.0,
            int_e=0.0,
            ds_e=0.0,
            hm_e=0.0,
            plp=0.0,
            success=False
        )
    except Exception as e:
        print(f"Error docking {mol_name}: {e}", file=sys.stderr)
        return DockingResult(
            mol_name=mol_name,
            energy=0.0,
            atdk_e=0.0,
            int_e=0.0,
            ds_e=0.0,
            hm_e=0.0,
            plp=0.0,
            success=False
        )


def read_first_mol2_entry(mol2_file: Path) -> str:
    """读取mol2文件中的第一个分子条目"""
    lines = []
    with open(mol2_file, 'r') as f:
        started = False
        for line in f:
            if line.startswith('@<TRIPOS>MOLECULE'):
                if started:
                    # 遇到第二个分子，停止
                    break
                started = True
            if started:
                lines.append(line)
    return ''.join(lines)


def cleanup_output_dir(output_dir: Path) -> None:
    """清理输出目录，只保留指定的文件"""
    keep_files = {'box.pdb', 'contact.pdb', 'GD2_HEME_fb.E.info', 'GD2_HEME_fb.mol2'}
    
    for item in output_dir.iterdir():
        if item.is_file() and item.name not in keep_files:
            item.unlink()


def write_summary_csv(results: List[DockingResult], output_file: Path) -> None:
    """写入汇总CSV文件"""
    # 按能量排序（越低越好）
    sorted_results = sorted(
        [r for r in results if r.success],
        key=lambda x: x.energy
    )
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Rank', 'Molecule', 'Energy', 'ATDK_E', 'INT_E', 'DS_E', 'HM_E', 'PLP'
        ])
        
        for rank, result in enumerate(sorted_results, 1):
            writer.writerow([
                rank,
                result.mol_name,
                f"{result.energy:.3f}",
                f"{result.atdk_e:.3f}",
                f"{result.int_e:.3f}",
                f"{result.ds_e:.3f}",
                f"{result.hm_e:.3f}",
                f"{result.plp:.3f}"
            ])
    
    # 同时写入失败的分子列表
    failed = [r for r in results if not r.success]
    if failed:
        failed_file = output_file.parent / "failed_molecules.txt"
        with open(failed_file, 'w') as f:
            for result in failed:
                f.write(f"{result.mol_name}\n")


def write_summary_mol2(results: List[DockingResult], output_file: Path) -> None:
    """将所有Top 1构象合并写入一个mol2文件"""
    # 按能量排序
    sorted_results = sorted(
        [r for r in results if r.success and r.mol2_content],
        key=lambda x: x.energy
    )
    
    with open(output_file, 'w') as f:
        for result in sorted_results:
            f.write(result.mol2_content)
            if not result.mol2_content.endswith('\n'):
                f.write('\n')


def _worker_wrapper(args):
    """进程池worker包装函数"""
    return run_single_docking(*args)


def main():
    parser = argparse.ArgumentParser(
        description="GalaxyDock2_HEME 批量对接脚本",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
    python batch_docking.py \\
        --protein data/protein/pig_CYP_merged.pdb \\
        --ligand_db data/mols/library.mol2 \\
        --output ./docking_results \\
        --heme_res_num 600 \\
        --chain A \\
        --n_workers 8
        """
    )
    
    parser.add_argument(
        "--protein", "-p",
        type=str,
        required=True,
        help="蛋白PDB文件路径（含血红素）"
    )
    parser.add_argument(
        "--ligand_db", "-l",
        type=str,
        required=True,
        help="配体库mol2文件路径（可包含多个分子）"
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        required=True,
        help="输出目录"
    )
    parser.add_argument(
        "--galaxydock_home", "-d",
        type=str,
        default=None,
        help="GalaxyDock2_HEME主目录（默认：自动检测）"
    )
    parser.add_argument(
        "--heme_res_num",
        type=str,
        default=None,
        help="血红素残基号"
    )
    parser.add_argument(
        "--chain",
        type=str,
        default=None,
        help="蛋白链ID"
    )
    parser.add_argument(
        "--n_workers", "-n",
        type=int,
        default=4,
        help="并行进程数（默认：4）"
    )
    parser.add_argument(
        "--conda_env",
        type=str,
        default="dock",
        help="Conda环境名（默认：dock）"
    )
    parser.add_argument(
        "--conda_base",
        type=str,
        default=None,
        help="Conda安装路径（默认：自动检测）"
    )
    parser.add_argument(
        "--keep_temp",
        action="store_true",
        help="保留临时拆分的mol2文件"
    )
    
    args = parser.parse_args()
    
    # 路径处理
    protein_path = Path(args.protein).resolve()
    ligand_db_path = Path(args.ligand_db).resolve()
    output_dir = Path(args.output).resolve()
    
    # 自动检测GalaxyDock2_HEME主目录
    if args.galaxydock_home:
        galaxydock_home = Path(args.galaxydock_home).resolve()
    else:
        # 假设脚本在 GalaxyDock2_HEME/script/ 目录下
        script_dir = Path(__file__).parent
        galaxydock_home = script_dir.parent
    
    # 自动检测conda路径
    if args.conda_base:
        conda_base = Path(args.conda_base)
    else:
        # 尝试常见路径
        for base in [
            Path("/opt/share/miniconda3"),
            Path.home() / "miniconda3",
            Path.home() / "anaconda3"
        ]:
            if base.exists():
                conda_base = base
                break
        else:
            print("错误：无法自动检测conda路径，请使用 --conda_base 指定", file=sys.stderr)
            sys.exit(1)
    
    # 验证文件存在
    if not protein_path.exists():
        print(f"错误：蛋白文件不存在: {protein_path}", file=sys.stderr)
        sys.exit(1)
    
    if not ligand_db_path.exists():
        print(f"错误：配体库文件不存在: {ligand_db_path}", file=sys.stderr)
        sys.exit(1)
    
    # 创建输出目录
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 创建临时目录用于拆分mol2
    db_name = ligand_db_path.stem
    temp_mol2_dir = output_dir / f".temp_mol2_{db_name}"
    
    print("=" * 60)
    print("GalaxyDock2_HEME 批量对接")
    print("=" * 60)
    print(f"蛋白文件: {protein_path}")
    print(f"配体库: {ligand_db_path}")
    print(f"输出目录: {output_dir}")
    print(f"并行进程数: {args.n_workers}")
    print(f"Conda环境: {args.conda_env}")
    print("=" * 60)
    
    # 1. 拆分mol2文件
    print("\n[1/3] 拆分配体库文件...")
    molecules = split_mol2_file(ligand_db_path, temp_mol2_dir)
    print(f"      共拆分出 {len(molecules)} 个分子")
    
    # 2. 并行对接
    print(f"\n[2/3] 开始并行对接 (使用 {args.n_workers} 个进程)...")
    
    results = []
    completed = 0
    failed = 0
    
    # 准备任务参数
    tasks = [
        (
            mol_name,
            mol2_path,
            protein_path,
            galaxydock_home,
            output_dir,
            args.heme_res_num,
            args.chain,
            args.conda_env,
            conda_base
        )
        for mol_name, mol2_path in molecules
    ]
    
    with ProcessPoolExecutor(max_workers=args.n_workers) as executor:
        futures = {
            executor.submit(_worker_wrapper, task): task[0]
            for task in tasks
        }
        
        for future in as_completed(futures):
            mol_name = futures[future]
            try:
                result = future.result()
                results.append(result)
                completed += 1
                
                if result.success:
                    print(f"      [{completed}/{len(molecules)}] {mol_name}: "
                          f"Energy={result.energy:.3f}, HM_E={result.hm_e:.3f}")
                else:
                    failed += 1
                    print(f"      [{completed}/{len(molecules)}] {mol_name}: FAILED")
                    
            except Exception as e:
                failed += 1
                completed += 1
                print(f"      [{completed}/{len(molecules)}] {mol_name}: ERROR - {e}")
                results.append(DockingResult(
                    mol_name=mol_name,
                    energy=0.0,
                    atdk_e=0.0,
                    int_e=0.0,
                    ds_e=0.0,
                    hm_e=0.0,
                    plp=0.0,
                    success=False
                ))
    
    # 3. 汇总结果
    print(f"\n[3/3] 汇总结果...")
    
    summary_csv = output_dir / f"{db_name}_summary.csv"
    summary_mol2 = output_dir / f"{db_name}_top1_poses.mol2"
    
    write_summary_csv(results, summary_csv)
    write_summary_mol2(results, summary_mol2)
    
    # 清理临时文件
    if not args.keep_temp and temp_mol2_dir.exists():
        shutil.rmtree(temp_mol2_dir)
    
    # 打印统计
    successful = len([r for r in results if r.success])
    print("\n" + "=" * 60)
    print("对接完成!")
    print("=" * 60)
    print(f"总分子数: {len(molecules)}")
    print(f"成功: {successful}")
    print(f"失败: {failed}")
    print(f"成功率: {successful/len(molecules)*100:.1f}%")
    print(f"\n结果文件:")
    print(f"  - 汇总CSV: {summary_csv}")
    print(f"  - Top1构象: {summary_mol2}")
    if failed > 0:
        print(f"  - 失败列表: {output_dir / 'failed_molecules.txt'}")
    print("=" * 60)


if __name__ == "__main__":
    main()
