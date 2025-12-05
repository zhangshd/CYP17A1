#!/usr/bin/env python

import sys
import numpy as np

PDBfmt     = '%-6s%5d %-4s %3s %1s%5s   %8.3f%8.3f%8.3f\n'
class PDB:
    def __init__(self, pdb_fn, read=True, read_het=True):
        self.pdb_fn = pdb_fn
        self.model_s = []
        self.use_model = False
        if read: self.read(read_het=read_het)
    def __repr__(self):
        return self.pdb_fn
    def __len__(self):
        return len(self.model_s)
    def __getitem__(self, i):
        return self.model_s[i]
    def read(self, read_het=True):
        read_model = True
        self.model_s = []
        i_model = 0
        resNo_prev = None
        chain_prev = None
        model = Model()
        self.model_s.append(model)
        #
        fp = open("%s"%self.pdb_fn)
        lines = fp.readlines()
        fp.close()
        #
        for line in lines:
            if line.startswith("MODEL"):
                self.use_model = True
                model_no = int(line.strip().split()[1])
                if len(model) != 0 and model.use:
                    i_model += 1
                if len(model) != 0 and model.use:
                    model = Model(model_no=model_no)
                    self.model_s.append(model)
                else:
                    model.model_no = model_no
                model.append(PDBline(line))
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                if line[16] not in ['A', ' ']:
                    continue
                if line.startswith("HETATM") and (not read_het):
                    continue
                #
                model.use = True
                #
                resNo = line[22:27]
                chain_id = line[21]
                #
                if resNo != resNo_prev or chain_id != chain_prev:
                    resNo_prev = resNo
                    chain_prev = chain_id
                    #
                    residue = Residue(line)
                    model.append(residue)
                residue.append(line)
            elif read_model:
                model.append(PDBline(line))
    def write(self, exclude_remark=True, exclude_symm=False, exclude_missing_bb=True, model_index=[],\
              remark_s=[]):
        model_index = range(len(self))
        #
        wrt = []
        wrt.extend(remark_s)
        #
        model_no = 0
        for j in model_index:
            if self.use_model:
                model_no += 1
                wrt.append("MODEL %4d\n"%(model_no))
            wrt.extend(self.model_s[j].write(exclude_remark=exclude_remark, exclude_symm=exclude_symm,\
                                             exclude_missing_bb=exclude_missing_bb))
            if self.use_model:
                wrt.append("ENDMDL\n")
        wrt.append("END\n")
        return wrt

class Model:
    def __init__(self, model_no=0):
        self.use = False
        self.model_no = model_no
        self.lines = []
        self.res_index = {}
    def append(self, X):
        self.lines.append(X)
        if X.isResidue():
            self.res_index[(X.chainID(), X.resNo_char())] = len(self.lines)-1
    def index(self, key):
        return self.res_index[key]
    def __getitem__(self, i):
        return self.lines[i]
    def __len__(self):
        return len(self.get_residues())
    def write(self, exclude_remark=False, exclude_symm=False, exclude_missing_bb=False,\
              exclude_nucl=False, exclude_SSbond=False, remark_s=[], chain_id=None):
        wrt = []
        wrt.extend(remark_s)
        for line in self.lines:
            if line.isResidue():
                if chain_id != None and chain_id != line.chainID():
                    continue
                if exclude_nucl and line.isAtom() and \
                        (line.resName().strip() in ['DA','DC','DG','DT','DU','A','C','G','T','U']):
                    continue
                if exclude_missing_bb and (not line.check_bb()):
                    continue
                wrt.append('%s'%line)
            else:
                if line.startswith("MODEL"):
                    continue
                elif line.startswith("END"):
                    continue
                elif line.startswith("TER"):
                    if len(wrt) != 0 and (not wrt[-1].startswith("TER")):
                        wrt.append("TER\n")
                    continue
                if line.startswith('REMARK 350') and (not exclude_symm):
                    wrt.append('%s'%line)
                    continue
                elif line.startswith('SSBOND') and (not exclude_SSbond):
                    wrt.append('%s'%line)
                    continue
                elif exclude_remark:
                    continue
                wrt.append('%s'%line)
        return wrt
    def __repr__(self):
        return ''.join(self.write())
    def get_residues(self, res_range=[]):
        lines = []
        for line in self.lines:
            if not line.isResidue():
                continue
            if len(res_range) != 0 and \
                    (line.resNo() not in res_range) and \
                    (line.resNo_char() not in res_range):
                continue
            lines.append(line)
        return lines
    def get_residue_lines(self, res_range=[]):
        lines = []
        for line in self.get_residues(res_range=res_range):
            lines.append('%s'%line)
        return lines

class Residue:
    def __init__(self, line):
        self._diso = False
        self._header = line[:6]
        self._resName = line[17:20]
        self._resNo = line[22:27]
        self._chainID = line[21]
        #
        self._R = []
        self._i_atm = []
        self._atmName = []
    def __len__(self):
        return len(self._R)
    def append(self, line):
        atmName = line[12:16].strip()
        if len(atmName) == 4:
            atmName = '%s%s'%(atmName[1:], atmName[0])
        self._atmName.append(atmName)
        self._i_atm.append(int(line[6:11]))
        self._R.append((float(line[30:38]),\
                        float(line[38:46]),\
                        float(line[46:54])))
    def isResidue(self):
        return True
    def isAtom(self):
        return self._header[:4] == 'ATOM'
    def isHetatm(self):
        return self._header == 'HETATM'
    def exists(self, atmName):
        return atmName in self._atmName
    def check_bb(self):
        stat = [False, False, False, False]
        if 'N' in self._atmName:
            stat[0] = True
        if 'CA' in self._atmName:
            stat[1] = True
        if 'C' in self._atmName:
            stat[2] = True
        if 'O' in self._atmName:
            stat[3] = True
        #
        if False in stat and not self._header == 'HETATM':
            return False
        else:
            return True
    def write(self):
        wrt = []
        if not self._diso:
            for i in range(len(self._R)):
                wrt.append(self[i])
        return ''.join(wrt)
    def __repr__(self):
        return self.write()
    def __getitem__(self, i):
        if isinstance(i, str):
            i = self.atmIndex(i)
        if len(self._atmName[i]) == 4:
            atmName = '%s%s'%(self._atmName[i][-1],self._atmName[i][:3])
        else:
            atmName = ' %s'%self._atmName[i]
        line = PDBfmt%(self._header, self._i_atm[i], atmName,\
                       self._resName, self._chainID, self._resNo,\
                       self._R[i][0], self._R[i][1], self._R[i][2])
        return line
    def resName(self):
        return self._resName
    def resNo(self):
        return int(self._resNo[:4])
    def resNo_char(self):
        return self._resNo
    def chainID(self):
        return self._chainID
    def atmName(self):
        return self._atmName
    def i_atm(self, atmName=None, atmIndex=None):
        if atmName != None:
            return self._i_atm[self.atmIndex(atmName)]
        elif atmIndex != None:
            return self._i_atm[atmIndex]
        else:
            return self._i_atm
    def R(self, atmName=None, atmIndex=None):
        if atmName != None:
            return self._R[self.atmIndex(atmName)]
        elif atmIndex != None:
            return self._R[atmIndex]
        else:
            return self._R
    def atmIndex(self, atmName):
        return self._atmName.index(atmName)
    def get_backbone(self):
        return [self._atmName.index("N"), self._atmName.index("CA"),\
                self._atmName.index("C"), self._atmName.index("O")]
    def get_heavy(self):
        heavy = []
        for i,atm in enumerate(self._atmName):
            if atm[0] != 'H':
                heavy.append(i)
        return heavy
    def get_sc(self):
        sc = []
        for i,atm in enumerate(self._atmName):
            if atm in ["N", "CA", "C", "O"]:
                continue
            if atm[0] == 'H':
                continue
            sc.append(i)
        return sc
    def get_CB(self):
        if self._resName == 'GLY':
            return [self._atmName.index("CA")]
        else:
            return [self._atmName.index("CB")]

class Atom(Residue):
    def __init__(self, line):
        self._header = line[:6]
        self._resName = line[17:20]
        self._resNo = line[22:27]
        self._chainID = line[21]
        #
        atmName = line[12:16].strip()
        if len(atmName) == 4:
            atmName = '%s%s'%(atmName[1:], atmName[0])
        self._atmName = atmName
        self._i_atm = int(line[6:11])
        self._R = np.array((float(line[30:38]),\
                            float(line[38:46]),\
                            float(line[46:54])))

    def __repr__(self):
        if len(self._atmName) == 4:
            atmName = '%s%s'%(self._atmName[-1],self._atmName[:3])
        else:
            atmName = ' %s'%self._atmName
        line = PDBfmt%(self._header, self._i_atm, atmName,\
                   self._resName, self._chainID, self._resNo,\
                   self._R[0], self._R[1], self._R[2])
        return line
    def R(self):
        return self._R
    def i_atm(self):
        return self._i_atm
    def atmName(self):
        return self._atmName

class PDBline:
    def __init__(self, line):
        self.line = line
    def __repr__(self):
        return self.line
    def isResidue(self):
        return False
    def isAtom(self):
        return False
    def isHetatm(self):
        return False
    def startswith(self, key):
        return self.line.startswith(key)
