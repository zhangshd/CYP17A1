import os
import glob
import sys
import chimera
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages
import AddH
import AddCharge
import WriteMol2


def main(lig_file):
    #open ligand file for reading
    chimera.openModels.open(lig_file)

    #chimera.update.checkForChanges()
    lig_file_stem = lig_file.split('/')[-1].split('.')[0]
    lig_file_parent = '/'.join(lig_file.split('/')[:-1])
    

    rc("del H")
    rc("addh")
    rc("addcharge all method gas")
    rc("write format mol2 #0 %s/%s.mol2"%(lig_file_parent,lig_file_stem))
    #rc("write format mol2 #0 test.mol2")

    rc('close all')

    print ('completed')
    
lig_file = sys.argv[3]

main(lig_file)
