import numpy as np
from ECFP0 import SurfaceStr


class BmxStr(SurfaceStr):
    def __init__(self):
        self.atom = []
        self.surfaceatom = []
        return

    def Inherit(self, str):
        self.natom = str.natom
        self.abc = str.abc
        self.lat = str.lat
        self.surfaceatom = str.surfaceatom[:]

        self.organicatom = []
        self.norganicatom =0
        for i in range(self.natom):
            self.atom.append(str.atom[i])
            if str.atom[i].ele < 18:
                self.organicatom.append(i)
                self.norganicatom = self.norganicatom +1

        self.bmx2D =np.zeros((self.natom,self.natom),dtype =np.int)
        for i in range(self.natom):
            for j in range(self.natom):
                self.bmx2D[i][j] = str.bmx2D[i][j]

        return

   

    def AddMetalBmx(self,str,nmetal,metalbmx,surfacemetal, metalsupercellinfo):
        if not hasattr(str, 'norganicatom'):
            str.norganicatom =0
            str.organicatom = []
            for i in range(str.natom):
                if str.atom[i].ele < 18:
                    str.organicatom.append(i)
                    str.norganicatom = str.norganicatom +1


        self.natom = str.norganicatom + nmetal
        self.lat = str.lat
        self.surfaceatom = [0 for i in range(self.natom)]

        self.norganicatom = str.norganicatom

        for i in str.organicatom:
            self.atom.append(str.atom[i])
        for i in range(nmetal):
            self.atom.append(surfacemetal[i])

        self.bmx2D =np.zeros((self.natom,self.natom),dtype =np.int)
        for i in range(0,str.norganicatom):
            for j in range(0,str.norganicatom):
                self.bmx2D[i][j] = str.bmx2D[i][j]
        for i in range(str.norganicatom,self.natom):
            for j in range(str.norganicatom,self.natom):
                self.bmx2D[i][j] = metalbmx[i-str.norganicatom][j-str.norganicatom]
       
        self.supercellinfo = {}
        for i in range(0,self.natom):
            for j in range(i+1,self.natom):
                if i >= str.norganicatom:
                    self.supercellinfo[i,j]= metalsupercellinfo[i-str.norganicatom,j-str.norganicatom]
                    self.supercellinfo[j,i]= metalsupercellinfo[j-str.norganicatom,i-str.norganicatom]
                else:
                    self.supercellinfo[i,j]= [999,[0,0,0]]
                    self.supercellinfo[j,i]= [999,[0,0,0]]



    def InitBondneed(self):
        bondneed = [0 for i in range(self.natom)]
        for i in range(self.natom):
            atom = self.atom[i]
            if atom.ele > 18:
                bondneed[i] = 0
            else:
                bondneed[i] = min(abs(atom.ele-2),abs(atom.ele-10),abs(atom.ele-18))
        self.bondneed = bondneed

    def CalBondneed(self):
        self.InitBondneed()
        for i in range(self.natom):
            if self.atom[i].ele > 18:
                continue
            else:
                self.bondneed[i] = self.bondneed[i] -sum(self.bmx2D[i][:])

    def RemoveMetalBond(self):
        self.surfaceatom = [0 for i in range(self.natom)]
        #self.InitBondneed()
        self.bmxsave =np.zeros((self.natom,self.natom),dtype =np.int)

        for i in range(self.natom):
            for j in range(self.natom):
                if self.bmx2D[i][j] > 0:
                    self.bmxsave[i][j] = self.bmx2D[i][j]
                    if self.atom[i].ele > 18:#or self.atom[j].ele > 18:
                        self.bmx2D[i][j] =0
                        self.bmx2D[j][i] =0
                        if self.atom[j].ele < 18:
                            self.surfaceatom[j]= 1
                    if self.atom[j].ele > 18:
                        self.bmx2D[i][j] =0
                        self.bmx2D[j][i] =0
                        if self.atom[i].ele < 18:
                            self.surfaceatom[i]= 1

        self.CalBondneed()
        self.bmx1D = []
        for line in self.bmx2D:
            self.bmx1D.extend(line)



    def GetSaveBmx(self):
        self.bmx2D =self.bmxsave
        self.bmx1D = []
        for line in self.bmx2D:
            self.bmx1D.extend(line)



    def GenSmilesName(self):

        substr = [[] for i in np.unique(self.group)]
        for id,atom in enumerate(str.atom):
            atom.id = id
            substr[str.group[id]-1].append(atom)
        substr = sorted(substr, key=lambda x:calmass(x), reverse=True)
        if flag == 1: allgroup.append((substr, str.bmx2D, str.lat,[],1,[]))
        if flag == 2: allgroup.append((substr, str.bmx2D, str.lat,str.bondneed,2,str.surfaceatom))

        return
