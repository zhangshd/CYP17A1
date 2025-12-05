Input for python script (In the case of apo and holo protein)

(For binary file, you should use apo protein with heme group, prepare heme_group mol2 file for partial charge assignment, add polar hydrogen to apo protein, and clean broken backbone.
Excecution using binary file is difficult so we recommend use python script (GalaxyDock2_HEME/script/run_GalaxyDock2_heme.py)

Please use python whose version is later than 3.4

-d path of GalaxyDock2-HEME directory (d), which use files, such as d/data and d/bin/ligdock

-p holo_protein; pdb format
If there is no specification of ligand residue number for center of binding pocket (docking box, default 22.5*22.5*22.5 angstrom), the first ligand residue number (except HOH water) is selected for center of binding pocket

-l ligand; mol2 format

-x {x coordinate of docking box} -y {y coordinate of docking box} -z {z coordinate of docking box}

--heme_res_num {a residue number of a heme cofactor}

--chain {a chain identifier}

Please check whether the chain identifier and the residue muber of heme cofactor match

if ligand format is not mol2 format, convert ligand format from smi or sdf or pdb to mol2 using openbabel

############################################################################
executeion of openbabel_requirement -> you should write obabel path like OBABEL_PATH = /.../obabel

For ligand preparation in python script, UCSF chimera should be installed
##########################################################################

other options can be found using "python run_GalaxyDock2_heme.py"

You can refer to examples/run_using_python/random_seed_0/linux_command.txt

You can test running docking in common cases by executing "python ../../script/run_GalaxyDock2_heme.py -p 2vhd.pdb -l MYT_sup.mol2 -d ../../../GalaxyDock2_HEME/ -x 14.107 -y -10.346 -z 19.504 --heme_res_num 401 --chain A" in examples/common_test/ directory

########################################################################
Outputs of python version and binary version can be slightly different because of a little bit different assignment of atom types and added hydrogen location of HEM.mol2

Finally sampled 100 ligand conformations sorted by their scores (outputs) are saved in GD2_HEME_fb.mol2, information of output ligand conformations is saved in GD2_HEME_fb.E.info (No need to care about RMSD values in this file)

Result of GalaxyDock2-HEME paper was conducted using preprocessed input files and the binary file.

If no random seed is assigned, random value of random seed will used for docking

If absolute path is too long, error can occur because of parsing constant length of absolute path. In that case, we recommend the relative path.
