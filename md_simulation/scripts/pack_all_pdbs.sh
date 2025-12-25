#!/bin/bash
# 每个分子的所有PDB文件放入独立文件夹，再整体打包

cd ../complexes/systems_gpu || exit 1
mkdir -p all_pdbs

for id in $(cat ligand_list.txt); do
  if [ -d "$id" ]; then
    mkdir -p all_pdbs/$id
    cp $id/*.pdb all_pdbs/$id/ 2>/dev/null
  fi
done

tar -czf all_pdbs_keyframes.tar.gz all_pdbs

echo "已将所有分子的PDB文件按文件夹归类并整体打包为all_pdbs_keyframes.tar.gz。"
