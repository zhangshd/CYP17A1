#!/usr/bin/env python

import sys
import os
import argparse
from pathlib import Path
import numpy as np

import lcs_modified_libpdb as libpdb
import subprocess as sp

script_dir = Path(__file__).parent

# Openbabel path should be set
OBABEL_PATH = Path('/opt/share/miniconda3/envs/dock/bin/obabel')

CHIMERA_SCRIPT = (script_dir/'chimera_proc_ligand.py').resolve()

assert CHIMERA_SCRIPT.exists()

def make_trivial_ligand_set():
    ion_list = ['K','CA','MG','FE','MN','BA','ZN','HG','NI','GD','IOD','BR','CL','SO4','SO3','PO4','NO3','CO3','SCN','CAC']
    cofactor_list = ['FMN','NAD','NDP']
    inorganic_list = ['FES','FAD']
    protein_residue_list = ['CSO','SEP','TRP','GLN','CAS']
    minor_ligand_list = ['ACT','ACE','MLI','PLM','GAI','EPE','OXY','FMT','UNX','CIT','BTB','TRS','DMS','MES','CAD','BME']    
    solvent_list = ['HOH','MOH','GOL','EDO','PEG','BU3','P6G','PGE','ACE','MPD','1PE','PG4','PGE','N8E','MRD','OCT','IPA','HEZ','POL']
    sugar_list = ['BOG','NAG','CM5','BNG','GLC','FRU','BMA','MAN','LMT']
    #special_ligand_list = ['CPS','H4B']
    subligand_list = ion_list + cofactor_list + inorganic_list + protein_residue_list + minor_ligand_list + solvent_list + sugar_list  
    galaxysite_subligand_list = 'HOH, PO3, BME, IMD, FMT, UNL, FES, EOH, CAS, OH, ACY, TRS, MPD, IOD, MLI, MLT, DTT, ACT, MSE, GOL, EDO, PEG, SO4, CL, ABA, TPO, SEP, UNX, MAN, NAG, NDG, NLE, YCM, DAL, IPA, HEZ, BU1, PG4, P6G, BR, OCS, CSO, MLY, PGE, PE4, 2PE, MES, KCX'.split(', ')
    subligand_set = set(subligand_list)
    galaxysite_subligand_set = set(galaxysite_subligand_list)
    subligand_set = subligand_set | galaxysite_subligand_set
    subligand_set.update({'NA','TOX','OMT','TLA'})
    
    return subligand_set

def preprocess_ligand(lig_fn, out_dir):
    out_ligand = out_dir/'input_ligand.mol2'
    
    if out_ligand.suffix != '.mol2':
        OBABEL_PATH = OBABEL_PATH.resolve()
        assert OBABEL_PATH.exists()
        sp.run([OBABEL_PATH, str(lig_fn), '-O', str(out_ligand)], check=True)
    else:
        sp.run(['cp', str(lig_fn), str(out_ligand)], check=True)
    
    sp.run(['chimera', '--nogui', CHIMERA_SCRIPT, str(out_ligand)], check=True)
    
    return

def preprocess_protein(pdb_fn, ligand_res_num = None, heme_res_num = None, chain = None):
    trivial_ligand_set = make_trivial_ligand_set()
    pdb_lines = filter(None,pdb_fn.read_text().splitlines())
    out_pdb_lines = []
    out_heme_lines = []
    out_ligand_lines = []
    select_chain = chain
    first_alternative_res = None
    heme_res_name = None
    out_heme_res_num = None
    for line in pdb_lines:
        if (not line.startswith('ATOM')) and (not line.startswith('HETATM')) and (not line.startswith('TER')):
            continue
        chain = line[21]
        if select_chain is None:
            select_chain = chain
        if chain != select_chain:
            continue
        
        if line.startswith('ATOM'):
            out_pdb_lines.append(line)
        elif line.startswith('TER'):
            out_pdb_lines.append(line)
        elif line.startswith('HETATM'):
            res_name = line[17:20].strip()
            res_num = line[22:26].strip()
            alternative_res = line[16]
            if res_name in ['HEM','HEC']:
                if (heme_res_num is not None) and int(res_num) != int(heme_res_num):
                    continue
                if heme_res_num is None:
                    heme_res_num = res_num
                heme_res_name = res_name
                out_heme_res_num = heme_res_num
                if alternative_res not in [' '] and first_alternative_res is None:
                    first_alternative_res = alternative_res
                    if first_alternative_res != 'A':
                        new_line = f'{line[:16]}A{line[17:]}'
                        out_pdb_lines.append(new_line)
                        out_heme_lines.append(new_line)
                        continue
                out_pdb_lines.append(line)
                out_heme_lines.append(line)
            elif ligand_res_num and int(res_num) == int(ligand_res_num):
                if alternative_res not in [''] and first_alternative_res is None:
                    first_alternative_res = alternative_res 
                    if first_alternative_res != 'A':
                        new_line = f'{line[:16]}A{line[17:]}'
                        out_ligand_lines.append(new_line)
                        continue
                out_ligand_lines.append(line)
            elif ligand_res_num is None and res_name not in trivial_ligand_set:
                ligand_res_num = res_num
                if alternative_res not in [''] and first_alternative_res is None:
                    first_alternative_res = alternative_res 
                    if first_alternative_res != 'A':
                        new_line = f'{line[:16]}A{line[17:]}'
                        out_ligand_lines.append(new_line)
                        continue
                out_ligand_lines.append(line)
        
    if out_heme_res_num is None or heme_res_name is None:
        raise Exception('Case 1. The heme cofactor is not in the chain selected or the first chain\n Case 2. The name of a heme cofactor should be HEM or HEC')
    
    return out_pdb_lines, out_heme_lines, out_ligand_lines, heme_res_name, out_heme_res_num, select_chain

def get_center_crd_from_holo_protein(out_ligand_lines):
    crds_list = []
    for line in out_ligand_lines:
        crd_x, crd_y, crd_z = map(float,(line[30:38],line[38:46],line[46:54]))
        crds_list.append([crd_x,crd_y,crd_z])
    
    crds_array = np.array(crds_list)
    
    return crds_array.mean(axis=0)

def replace_atom_name(out_heme_mol2_lines):
    atom_name_mol2_dict = {'NA': 'N.pl3', 'NB': 'N.pl3', 'NC': 'N.pl3', 'ND': 'N.pl3'}
    processed_out_heme_mol2_lines = []
    read_count = False
    for line in out_heme_mol2_lines:
        if line.startswith('@<TRIPOS>BOND'):
            read_count = False        
        elif read_count and line.split()[1] in ['HA','HB','HC','HD']:                
            line = line.replace(line.split()[1] + ' ', line.split()[1] + 'E')
        elif read_count and line.split()[1] in atom_name_mol2_dict:
            atom_name = line.split()[1]
            line = line.replace(line.split()[5],atom_name_mol2_dict[atom_name])                
        elif line.startswith('@<TRIPOS>ATOM'):
            read_count = True
        processed_out_heme_mol2_lines.append(line)
        
    return processed_out_heme_mol2_lines

def replace_partial_charge(out_heme_mol2_lines):
    read_count = False
    processed_out_heme_mol2_lines = []
    
    charge_dict = {'O1A':-0.76, 'O2A':-0.76, 'CGA':0.62, 'HBA1':0.09, 'CBA':-0.28, 'HBA2':0.09, 'HAA1':0.09,
            'CAA':-0.18, 'HAA2':0.09, 'HMA1':0.09, 'HMA2':0.09, 'HMA3':0.09, 'CMA':-0.27, 'HMB1':0.09, 'HBM2':0.09,
            'HBM3':0.09, 'CMB':-0.27, 'HBB1':0.21, 'HBB2':0.21, 'CBB':-0.42, 'CAB':-0.15, 'HAB':0.15, 'CMC':-0.27,
            'HMC1':0.09, 'HMC2':0.09, 'HMC3':0.09, 'CAC':-0.15, 'HAC':0.15, 'CBC':-0.42, 'HBC1':0.21, 'HBC2':0.21,
            'CMD':-0.27, 'HMD1':0.09, 'HMD2':0.09, 'HMD3':0.09, 'CAD':-0.18, 'HAD1':0.09, 'HAD2':0.09, 'CBD':-0.28,
            'HBD1':0.09, 'HBD2':0.09, 'CGD':0.62, 'O1D':-0.76, 'O2D':-0.76, 'C1A':0.12, 'C2A':-0.06, 'C3A':-0.06,
            'C4A':0.12, 'NA':-0.18, 'C1B':0.12, 'C2B':-0.06, 'C3B':-0.06, 'C4B':0.12, 'ND':-0.18, 'C1C':0.12,
            'C2C':-0.06, 'C3C':-0.06, 'C4C':0.12, 'NC':-0.18, 'C1D':0.12, 'C2D':-0.06, 'C3D':-0.06, 'C4D':0.12,
            'NB':-0.18, 'CHA':-0.10, 'CHB':-0.10, 'CHC':-0.10, 'CHD':-0.10, 'HAE':0.10, 'HBE':0.10, 'HCE':0.10,
            'HDE':0.10, 'FE':0.24}
    
    for line in out_heme_mol2_lines:
        if line.startswith('@<TRIPOS>BOND'):
            read_count = False                   
        elif read_count:
            for atm_name in charge_dict:
                if atm_name == line.split()[1].strip():
                    replaced_charge = format(charge_dict[atm_name],"5.4f")
                    previous_charge = line.split()[-1]
                    line = replaced_charge.join(line.rsplit(previous_charge,1))
        elif line.startswith('@<TRIPOS>ATOM'):
            read_count = True
        processed_out_heme_mol2_lines.append(line)
        
    return processed_out_heme_mol2_lines

def make_hydrogen_lines(processed_out_heme_mol2_lines, out_heme_res_num, residue_name, select_chain):
    read_count = False
    hydrogen_lines= []
    for line in processed_out_heme_mol2_lines:                    
        if line.startswith('@<TRIPOS>BOND'):
            break
        elif read_count and line[47] == 'H':
            hydrogen_lines.append(line)             
        elif line.startswith('@<TRIPOS>ATOM'):
            read_count = True
    
    added_hydrogen_lines=[]
    for line in hydrogen_lines:
        atomnumber, atom, x_crd, y_crd, z_crd = line.split()[:5]
        x_crd = format(float(x_crd),".3f")
        y_crd = format(float(y_crd),".3f")
        z_crd = format(float(z_crd),".3f")
        # hydrogen_line = 'HETATM   44  HHA HEM A                                  1.00  0.00           H  '
        hydrogen_line = ''.join(['HETATM ',
                                 f'{atomnumber:>4}',
                                 ' ',
                                 f'{atom:>4}',
                                 ' ',
                                 f'{residue_name:>3}',
                                 ' ',
                                 select_chain,
                                 f'{out_heme_res_num:>4}',
                                 ' '*4,
                                 f'{x_crd:>8}',
                                 f'{y_crd:>8}',
                                 f'{z_crd:>8}'
                                 ])

        added_hydrogen_lines.append(hydrogen_line)

    return added_hydrogen_lines

def create_heme_mol2_file(out_heme_lines, heme_res_name, out_heme_res_num, out_dir, select_chain):
    out_heme_pdb = out_dir/f'{heme_res_name}.pdb'
    out_heme_mol2 = out_heme_pdb.with_suffix('.mol2')
    
    out_heme_pdb.write_text('\n'.join(out_heme_lines) + '\n')
    
    sp.run(['chimera', '--nogui', CHIMERA_SCRIPT, str(out_heme_pdb)], check=True)

    out_heme_mol2_lines = out_heme_mol2.read_text().splitlines()
    processed_out_heme_mol2_lines = replace_atom_name(out_heme_mol2_lines)
    processed_out_heme_mol2_lines = replace_partial_charge(processed_out_heme_mol2_lines)
    added_hydrogen_lines = make_hydrogen_lines(processed_out_heme_mol2_lines, out_heme_res_num, heme_res_name, select_chain)
    
    out_heme_mol2.write_text('\n'.join(processed_out_heme_mol2_lines))
    
    return added_hydrogen_lines

def clean_pdb(out_contact_pdb):
    pdb = libpdb.PDB(str(out_contact_pdb))
    
    missing_bb = False
    for residue in pdb[0].get_residues():
        if residue.isHetatm():
            continue
        else:
            if not residue.check_bb():
                missing_bb = True
    if missing_bb:
        wrt = pdb.write()
        with out_contact_pdb.open('w') as f:
            f.writelines(wrt)
            
    return

def make_contact_pdb(out_pdb_lines, added_hydrogen_lines, out_dir):
    out_contact_pdb = out_dir/'contact.pdb'
    out_contact_pdb.write_text('\n'.join(out_pdb_lines+added_hydrogen_lines+['END']))
    clean_pdb(out_contact_pdb)
    
    return

def create_gd2_heme_script(args, grid_box_cntr, heme_res_name, random_seed, home_dir, out_dir):
    out_script = out_dir/'gd2_heme.in'
    heme_mol2 = f'{heme_res_name}.mol2'
    infile_pdb = 'contact.pdb'
    infile_ligand = 'input_ligand.mol2'
    box_cntr = ' '.join(grid_box_cntr)
    script_lines = ['!==============================================',
                    '! I/O Parameters',
                    '!==============================================',
                    f'data_directory    {str(home_dir)}/data/',
                    f'infile_pdb        {str(infile_pdb)}',
                    f'infile_ligand     {str(infile_ligand)}',
                    f'top_type          polarh',
                    f'fix_type          all',
                    f'ligdock_prefix    GD2_HEME',
                    '!==============================================',
                    '! Grid Options',
                    '!==============================================',
                    f'grid_box_cntr     {box_cntr}',
                    f'grid_n_elem       {args.n_elem_x} {args.n_elem_y} {args.n_elem_z}',
                    f'grid_width        0.375',
                    '!==============================================',
                    '! Energy Parameters',
                    '!==============================================',
                    f'weight_type       GalaxyDock2',
                    '!==============================================',
                    '! Initial Bank Parameters',
                    '!==============================================',
                    f'first_bank               rand',
                    f'max_trial                50000',
                    f'e0max                    1000.0',
                    f'e1max                    1000000.0',
                    f'infile_mol2_topo {heme_mol2} {heme_res_name}',
                    f'random {random_seed}']
    
    if random_seed is None:
        script_lines.pop()
    out_script.write_text('\n'.join(script_lines))
    
    return

def prepare_input_files(args):
    home_dir = Path(args.home_dir).resolve()
    pdb_fn = Path(args.pdb_fn).resolve()
    lig_fn = Path(args.lig_fn).resolve()
    out_dir = Path(args.out_dir).resolve()
    preprocess_ligand(lig_fn, out_dir)
    out_pdb_lines, out_heme_lines, out_ligand_lines, heme_res_name, out_heme_res_num, select_chain = preprocess_protein(pdb_fn, args.lig_res_num, args.heme_res_num, args.chain)

    added_hydrogen_lines = create_heme_mol2_file(out_heme_lines, heme_res_name, out_heme_res_num, out_dir, select_chain)
    make_contact_pdb(out_pdb_lines, added_hydrogen_lines, out_dir)
    
    grid_box_cntr = (args.cntr_x, args.cntr_y, args.cntr_z)
    
    if all(grid_box_cntr):
        grid_box_cntr = map(float, grid_box_cntr)
    else:
        grid_box_cntr = get_center_crd_from_holo_protein(out_ligand_lines)
    
    grid_box_cntr = map(lambda x: round(x,3), grid_box_cntr)
    grid_box_cntr = map(str, grid_box_cntr)
    
    create_gd2_heme_script(args, grid_box_cntr, heme_res_name, args.random_seed, home_dir, out_dir)
    
    return

def run_GalaxyDock2_HEME(args):
    sys.stdout.write("Running GalaxyDock....\n") 
    home_dir = Path(args.home_dir).resolve()
    pdb_fn = Path(args.pdb_fn).resolve()
    lig_fn = Path(args.lig_fn).resolve()
    out_dir = Path(args.out_dir).resolve()
    gd2_heme_in_fn = Path(out_dir/'gd2_heme.in')
    gd2_heme_log = Path(out_dir/'gd2_heme.log')
    
    if not home_dir.exists():
        print(f'ERROR: Please check path of GalaxyDock2 directory, {args.home_dir}')
        return False
    if not pdb_fn.exists():
        print(f'ERROR: Please check path of your receptor structure file, {args.pdb_fn}')
        return False
    if not lig_fn.exists():
        print(f'ERROR: Please check path of your ligand structure file, {args.lig_fn}')
        return False
    
    if args.prep:
        prepare_input_files(args)
    
    exec_gd2_heme = f'{home_dir}/bin/ligdock'
    #
    with gd2_heme_log.open('w') as f:
        sp.run([exec_gd2_heme,str(gd2_heme_in_fn)], stdout=f, check=True, cwd=str(out_dir))
    #
    sys.stdout.write("Done.\n") 

def main():
    current_path = os.getcwd()
    opt = argparse.ArgumentParser(description=\
            '''GalaxyDock2-HEME: Protein-Ligand Docking for heme proteins''',
            formatter_class=argparse.RawTextHelpFormatter)
    #
    opt.add_argument('-d', '--dir', dest='home_dir', metavar='GalaxyDock_HOME', required=True,\
                     help='Path to GalaxyDock_BP2 directory')
    opt.add_argument('-p', '--pdb', dest='pdb_fn', metavar='PDB', required=True,\
                     help='Protein receptor structure file in PDB format')
    opt.add_argument('-l', '--lig', dest='lig_fn', metavar='MOL2', required=True,\
                     help='Ligand structure file in MOL2 format')
    #
    opt.add_argument('-size_x', dest='size_x', metavar='SIZE_X', default='22.5',\
                     help='size in the X dimension (Angstroms) [Default: 22.5]')
    opt.add_argument('-size_y', dest='size_y', metavar='SIZE_Y', default='22.5',\
                     help='size in the Y dimension (Angstroms) [Default: 22.5]')
    opt.add_argument('-size_z', dest='size_z', metavar='SIZE_Z', default='22.5',\
                     help='size in the Z dimension (Angstroms) [Default: 22.5]')
    #
    opt.add_argument('-x', dest='cntr_x', metavar='CENTER_X', default=None,\
                     help='X coordinate of the binding site center')
    opt.add_argument('-y', dest='cntr_y', metavar='CENTER_Y', default=None,\
                     help='Y coordinate of the binding site center')
    opt.add_argument('-z', dest='cntr_z', metavar='CENTER_Z', default=None,\
                     help='Z coordinate of the binding site center')
    #
    opt.add_argument('--random_seed', dest='random_seed', metavar='RANDOM_SEED', default=None,\
                     help='Docking random seed for reproducibility. [Default: None]')
    opt.add_argument('--out_dir', dest='out_dir', metavar='OUT_DIR', default=current_path,\
                     help='Output directory. [Default: current_path]')
    opt.add_argument('--prep', dest='prep', metavar='PREP', default=True,\
                     help='Prepare input data, True or False [Default: True]')
    #
    opt.add_argument('--lig_res_num', dest='lig_res_num', default=None,\
                     help='Choose ligand_residue_number for center of search space without manual center of box, default for the first ligand [Default: None]')
    opt.add_argument('--heme_res_num', dest='heme_res_num', default=None,\
                     help='Choose heme_residue_number in the case of multiple heme exist [Default: None]')
    opt.add_argument('--chain', dest='chain', default=None,\
                     help='Choose protein chain identifier, default for the first chain [Default: None]')
    #
    if len(sys.argv) == 1:
        opt.print_help()
        return
    #
    args = opt.parse_args()
    #
    args.n_elem_x = '%d'%(int(float(args.size_x)/0.375) + 1)
    args.n_elem_y = '%d'%(int(float(args.size_y)/0.375) + 1)
    args.n_elem_z = '%d'%(int(float(args.size_z)/0.375) + 1)
    #
    run_GalaxyDock2_HEME(args)

    os.chdir(os.path.join(current_path))
    #check_call('rm -r %s'%(temp_dir),shell = True) 
    
if __name__=='__main__':
    main()
