import glob
import os
from subprocess import check_output
import sys

import torch
import torch.optim as optim
import torch.nn.functional as F
import pickle

def read_info(info_file):
    info_tuple_list = []
    with open(info_file,'r') as f:
        info_file_lines = f.readlines()
    
    for line in info_file_lines:
        if line.startswith('!') or line.startswith('Bank'):
            continue
        else:
            tokens=line.split()
            
            # RMSD (Calculated by Hungarian algorithm implemented in GalaxyDock2),
            # total score (total energy),
            # non-aromatic nitrogen coordination bond potential,
            # aromatic nitrogen coordination bond potential,
            # sulfur coordination bond potential,
            # oxygen coordination bond potential
            info_list = [tokens[2], tokens[1], tokens[8], tokens[9], tokens[10], tokens[11]]

            for i in range(len(info_list)):
                if '*' in info_list[i]:
                    info_list[i] = 9999999

            info_tuple = tuple(float(s) for s in info_list)
            info_tuple_list.append(info_tuple)
    
    if len(info_tuple_list) != 101:
        raise Exception('bank number is incorrect')
    
    return info_tuple_list

def coord_check(info_tuple):
    coord_tuple = info_tuple[2:6]
    min_val = 0.0
    min_index = 0
    for i, val in enumerate(coord_tuple):
        if val < min_val:
            min_val = val
            min_index = i+1
    
    return min_index

def get_local_opt_energy(info_file):
    ref_dir = '/'.join(info_file.split('/')[:-1]) + '/'
    opt_energy_file = ref_dir + 'opt_out.txt' # output log of local optimization
    file_name = ref_dir + 'out.mol2' # locally optimized pose using perturbation and local minimization by GalaxyDock2
    merge = ref_dir + 'merged_ligand.mol2' # input ligand
    
    cmd_rel_path = os.path.relpath('/applic/OpenBabel/2.4.1/bin/obrms') # OpenBabel's obrms path   
    cmd1 = '%s -firstonly'%(cmd_rel_path)
    output = check_output('%s %s %s'%(cmd1,file_name,merge), shell=True)
    with open(opt_energy_file) as f:
        lines = f.readlines()
    read_check = False
    for line in lines:        
        if line.startswith('After optimization'):
            read_check = True
            total_energy = line.split()[3]
        elif read_check and len(line.split()) == 4:
            tokens = line.split()
            break

    info_list = [0.0, total_energy, tokens[0], tokens[1], tokens[2], tokens[3]]
    true_info_list = [float(s) for s in info_list]
    
    # Check whether RMSD of locally optimized pose is in 2.0
    if float(output.split()[-1]) < 2.0:
        rmsd_check = True
    else:
        rmsd_check = False

    return rmsd_check, true_info_list

def train_check(info_tuple_list):
    rmsd_cutoff = 1.8
    if coord_check(info_tuple_list[100]):
        if info_tuple_list[0][0] < rmsd_cutoff:
            return 0
        else:
            return 1
    else:
        if info_tuple_list[0][0] < rmsd_cutoff:
            return 2
        else:
            return 3

def prep_for_obj_func(train_check_value, info_tuple_list):    
    rmsd_cutoff = 1.8
    ref_index_list = []   
    decoy_index_list = []
    # 100: GalaxyDock2-HEME score of locally optimized crystal pose only for training
    coord_index = coord_check(info_tuple_list[100])
    max_decoy_num = 5

    if train_check_value == 0:
        ref_index_list.append(0)        
        for i, info_tuple in enumerate(info_tuple_list[:100]):
            if coord_check(info_tuple) != coord_index:         
                decoy_index_list.append(i)
            if len(decoy_index_list) == max_decoy_num:
                break                

        return ref_index_list, decoy_index_list
    
    elif train_check_value == 1:
        for i, info_tuple in enumerate(info_tuple_list[:100]):
            if coord_check(info_tuple) == coord_index and info_tuple[0] < rmsd_cutoff:         
                ref_index_list.append(i)
                break      
        
        for i, info_tuple in enumerate(info_tuple_list[:100]):
            if coord_check(info_tuple) != coord_index:         
                decoy_index_list.append(i)
            if len(decoy_index_list) == max_decoy_num:
                break        

        if not ref_index_list:
            ref_index_list.append(100)

        return ref_index_list, decoy_index_list

    elif train_check_value == 2:
        ref_index_list.append(0)
        for i, info_tuple in enumerate(info_tuple_list[:100]):
            if coord_check(info_tuple) != coord_index:         
                decoy_index_list.append(i)
            if len(decoy_index_list) == max_decoy_num:
                break

        return ref_index_list, decoy_index_list

    elif train_check_value == 3:
        for i, info_tuple in enumerate(info_tuple_list[:100]):
            if coord_check(info_tuple) == coord_index and info_tuple[0] < rmsd_cutoff:         
                ref_index_list.append(i)
                break
        
        for i, info_tuple in enumerate(info_tuple_list[:100]):
            if coord_check(info_tuple) != coord_index:         
                decoy_index_list.append(i)
            if len(decoy_index_list) == max_decoy_num:
                break

        if not ref_index_list:
            ref_index_list.append(100)

        return ref_index_list, decoy_index_list

def info_tuple_to_formula(info_tuple, initial_weight):    
    const = info_tuple[1] - sum(info_tuple[2:]) # sum of other score (energy) terms constant during weight parameter optmization 
    
    raw_energy_nc = info_tuple[2]/abs(initial_weight[0])
    raw_energy_nr = info_tuple[3]/abs(initial_weight[1])
    raw_energy_s = info_tuple[4]/abs(initial_weight[2])
    raw_energy_o = info_tuple[5]/abs(initial_weight[3])       

    return raw_energy_nc, raw_energy_nr, raw_energy_s, raw_energy_o, const

def main():        
    ref_dir = '../%s_docking/'%(sys.argv[1]) # main directory of output files

    info_file_list = glob.glob(os.path.join(ref_dir,'*/GD2_fb.E.info')) # prefix_fb.E.info
    
    with open('pickles/temp_weight_list.pickle', 'rb') as f:
        initial_weight = pickle.load(f)

    W = torch.tensor([[float(s)] for s in initial_weight], requires_grad=True)
    
    raw_energy_diff_list = []
    const_diff = []

    for info_file in info_file_list:       
        info_tuple_list = read_info(info_file)
        train_check_value = train_check(info_tuple_list)

        ref_index_list, decoy_index_list = prep_for_obj_func(train_check_value, info_tuple_list)
        if ref_index_list[0] == 100:
            rmsd_check, true_info_list = get_local_opt_energy(info_file)

            if rmsd_check:
                info_tuple_list[100] = true_info_list

        if decoy_index_list and ref_index_list:
            out_diff = []
            for ref_index in ref_index_list:
                first_out = info_tuple_to_formula(info_tuple_list[ref_index], initial_weight)
                for decoy_index in decoy_index_list:
                    second_out = info_tuple_to_formula(info_tuple_list[decoy_index], initial_weight)
                    out_diff = []
                    for f, s in zip(first_out, second_out):
                        out_diff.append(f-s)

                    const_diff.append([out_diff[4]])
                    raw_energy_diff_list.append(out_diff[:4])
                        
    energy_var = torch.FloatTensor(raw_energy_diff_list)
    energy_const = torch.FloatTensor(const_diff)   
    optimizer = optim.SGD([W], lr=0.00001)
    
    nb_epochs = 100    
    margin = 2.0
   
    for _ in range(1, nb_epochs+1):        
        precost = energy_var.matmul(torch.abs(W)) + energy_const + margin
        
        cost1 = torch.sum(F.relu(precost))           

        optimizer.zero_grad()
        cost1.backward()
        optimizer.step()

        min_cost1 = 99999
        if cost1.item() < min_cost1:
            min_cost1 = cost1.item()
            min_weight = torch.abs(W).reshape(4).tolist()

    print(min_weight)
    print('min_cost: %.5f'%(min_cost1))
    
    with open('pickles/temp_weight_list.pickle', 'wb') as f:
        pickle.dump(min_weight, f)

main()