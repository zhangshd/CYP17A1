#!/bin/bash
# export_key_frames.sh
# Batch export key frames (min, heat, nvt, eq, prod) as PDB without water for all ligands
# Author: zhangshd
# Date: 2024-06-15

LIGAND_LIST="../complexes/systems_gpu/ligand_list.txt"
STAGES=(min1 min2 heat eq1 eq2_01 eq2_02 eq2_03 eq2_04 eq2_05 eq2_06 eq2_07 eq2_08 eq2_09 eq2_10)

while read -r id; do
    cd ../complexes/systems_gpu/"$id" || continue
    echo "Processing $id..."
    for stage in "${STAGES[@]}"; do
        rst="${stage}.rst7"
        if [ -f "$rst" ]; then
            cpptraj -p complex_solv.prmtop > /dev/null <<EOF
trajin $rst
strip :WAT,HOH
trajout ${stage}.pdb pdb
go
EOF
        fi
    done
    if [ -f prod.nc ]; then
        cpptraj -p complex_solv.prmtop > /dev/null <<EOF
trajin prod.nc 1 5000 500
strip :WAT,HOH
trajout prod_sampled.pdb pdb
go
EOF
    fi
    cd - > /dev/null
done < "$LIGAND_LIST"