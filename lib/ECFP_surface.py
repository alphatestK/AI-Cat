
from allstr_new import allstr
from structure_new import Str
import numpy as np
import hashlib
import zlib
import ctypes
from ctypes import pointer
from babelfunc import *


def wrapGenECFPname(datain):
    Str,savefp = datain
    #print Str.sminame
    ECFPname,FP = Str.GenECFPname()
    if savefp:
        return ECFPname,FP
    else:
        return ECFPname

def GetInfofromID(atomid):
    _atominfo = atomid.split(',')
    _id = int(_atominfo[0])
    _supercellinfo =[int(_atominfo[1]), int(_atominfo[2]),int(_atominfo[3])]
    return _id,_supercellinfo


def one_of_k_encoding(x, allowable_set):
    if x not in allowable_set:
        raise Exception("input {0} not in allowable set{1}:".format(x, allowable_set))
    return map(lambda s: x == s, allowable_set)

def one_of_k_encoding_unk(x, allowable_set):
    """Maps inputs not in the allowable set to the last element."""
    if x not in allowable_set:
        x = allowable_set[-1]
    return map(lambda s: x == s, allowable_set)


class SurfaceStr(Str):
    def GenFeature(self):
        self.AddAtomID()
        self.GetAtomInfo()
        for iatom,atom in enumerate(self.atom):
            atom.surface = self.surfaceatom[iatom]
            atom.surfacemetal = self.atom[-1].elesymbol
            atom.feature = self.Atom_Features(atom)

        #self.NeighbourList()
        self.NeighbourList_periodic()

    def Atom_Features(self,atom):
        return np.array(one_of_k_encoding_unk(atom.elesymbol,['C', 'O', 'H','N', 'Cu','Pt','Rh','Zn', 'Unknown']) +
        #return np.array(one_of_k_encoding_unk(atom.elesymbol,['C', 'O', 'H', 'Cu', 'Unknown']) +
                    one_of_k_encoding_unk(atom.expbond, [0, 1, 2, 3, 4, 5,6,7,8,9,10,11,12,'Unknown']) +
                    one_of_k_encoding_unk(atom.imph, [0, 1, 2, 3, 4, 5,'Unknown'])+
                    one_of_k_encoding_unk(atom.surfacemetal, ['Cu','Pt','Rh','Zn','Unknown'])+
                    #one_of_k_encoding_unk(atom.value, [-2, -1, 0, 1, 2, 3])+
                    #one_of_k_encoding_unk(atom.surface, ['single', 'bridge','hcp', 'fcc','4','unkonwn'])+
                    [0])


    def Bond_Features(self,bondorder):
        return np.array(one_of_k_encoding(bondorder,[1,2,3,4,5]))

#    def GetNeighbour(self,iatom):
#        for i in range(self.natom):
#            if self.bmx2D[iatom][i]>0:
#                self.atom[iatom].neighbour.append()
#


    def NeighbourList(self):
        for i in range(self.natom):
            self.atom[i].nbatom=[]
            self.atom[i].bond=[]
            self.atom[i].nb=[]

        for i in range(self.natom):
            for j in xrange(i+1,self.natom):
                if self.bmx2D[i][j]>0:
                    # do i need skip H?
                    self.atom[j].nbatom.append(self.atom[i])
                    self.atom[i].nbatom.append(self.atom[j])
                    #self.atom[j].bond.append(self.Bond_Features(self.bmx2D[i][j]))
                    #self.atom[i].bond.append(self.Bond_Features(self.bmx2D[i][j]))
                    self.atom[j].nb.append([self.bmx2D[i][j],i])
                    self.atom[i].nb.append([self.bmx2D[i][j],j])

    def NeighbourList_periodic(self):
        if not hasattr(self, 'supercellinfo'):
            self.PeriodicInfo()

        for i in range(self.natom):
            for j in xrange(i+1,self.natom):
                if not (self.atom[i].ele > 18 and self.atom[j].ele > 18):
                    self.supercellinfo[i,j][1] = [0,0,0]
                    self.supercellinfo[j,i][1] = [0,0,0]

        self.PeriodicSitedefination()

        for i in range(self.natom):
            self.atom[i].nbatom=[]
            self.atom[i].bond=[]
            self.atom[i].nb=[]
            self.atom[i].id_periodic  = str(self.atom[i].id)+ str(',0') + str(',0')+ str(',0')

        for i in range(self.natom):
            for j in xrange(i+1,self.natom):
                if self.bmx2D[i][j]>0:
                    # do i need skip H?
                    self.atom[j].nbatom.append(self.atom[i])
                    self.atom[i].nbatom.append(self.atom[j])
                    #self.atom[j].bond.append(self.Bond_Features(self.bmx2D[i][j]))
                    #self.atom[i].bond.append(self.Bond_Features(self.bmx2D[i][j]))
#                    if (self.atom[i].ele > 18 and self.atom[j].ele > 18):
                    self.atom[i].nb.append([self.bmx2D[i][j],j,self.supercellinfo[i,j][1]])
                    self.atom[j].nb.append([self.bmx2D[j][i],i,self.supercellinfo[j,i][1]])


    def GetAllAdsorpsite(self):
        self.adsorpatom = {}
        for i in range(self.natom):
            for j in range(self.natom):
                if self.bmx2D[i][j] > 0:
                    if (self.surfaceatom[i] != 0) and (self.atom[j].ele >18):
                        if i not in self.adsorpatom.keys():
                            self.adsorpatom[i]= [j]
                        elif j not in self.adsorpatom[i]:
                            self.adsorpatom[i].append(j)


    def PeriodicSitedefination(self):
        if not hasattr(self,'adsorpatom'):
            self.GetAllAdsorpsite()
            #print ('no adsorpatom info')
            return
        #print self.adsorpatom
        for adid in self.adsorpatom.keys():
            scoresave = 999
            for satom in self.adsorpatom[adid]:
                _id = str(self.atom[satom].id)+ str(',0') + str(',0')+ str(',0')
                _vec = [0,0,0]
                _combine = [_id]
                for j in self.adsorpatom[adid]:
                    if j == satom : continue
                    _vec = np.add(_vec,[abs(x) for x in self.supercellinfo[satom,j][1]])
                    supercellinfo = self.supercellinfo[satom,j][1]
                    _newid = str(j)+',' +str(supercellinfo[0])+',' +str(supercellinfo[1])+',' +str(supercellinfo[2])
                    _combine.append(_newid)
                score = sum(_vec)
                if score < scoresave:
                    combinesave = _combine
                    scoresave = score

            for atomid in combinesave:
                _id,_supercellinfo = GetInfofromID(atomid)
                self.supercellinfo[adid,_id][1] = [ x  for x in _supercellinfo]
                self.supercellinfo[_id,adid][1] = [ -x for x in _supercellinfo]
            #print adid,combinesave


    def SortbyDegree(self,):
        self.atom.sort(key = lambda x:len(x.neighbour))

    def AllAtomFeature(self):
        atomfeature = []
        for atom in self.atom:
            atomfeature.append(atom.feature)
            atomfeature.append(atom.bond)
      
    def AllBondFeature(self):
        atomfeature = []
        for atom in self.atom:
            for bond in atom.bond:
                bondfeature.append(bond)
                bondfeature.append(bond)

    def GenECFPname(self,depth=3,dim=10000000,savefp = False):
        FP = ECFP(self,depth,dim)
        FP.GenECFP()
        string = ''
        for index in FP.allindex:
            #print index
            string = string + str(len(FP.IndextoFp[index]))
            string = string+ str(index)
        self.ECFPname= string
        if savefp:
            self.FP = FP
        
        #print self.ECFPname
        return string,FP

    def CheckMin(self,flag = 1):
        sqnatm = self.natom**2

        self.calCtypes()
        #program='/home10/kpl/pymodule/Lib/Lib_fillbond_new/checkminbond.so'
        program='/home10/kpl/pymodule/Lib/Lib_fillbond_surface/checkminbond.so'
        Lminstr = pointer(ctypes.c_bool(0))
        bmatrix = pointer((ctypes.c_int*sqnatm)(*[0 for i in range(sqnatm)]))
        bondneed = pointer((ctypes.c_int*(self.natom))(*[0 for i in range(self.natom)]))
        surface = pointer((ctypes.c_int*(self.natom))(*[0 for i in range(self.natom)]))

        checkmin = ctypes.cdll.LoadLibrary(program)
        if flag == 0:
            checkmin.judgebond_(self.c_natm,self.c_iza,self.c_xa,self.c_rv,Lminstr,bmatrix,bondneed)
        elif flag ==1 :
            #print 'into surface'
            checkmin.judgebondsurface_(self.c_natm,self.c_iza,self.c_xa,self.c_rv,Lminstr,bmatrix,bondneed,surface)
            #print 'finish fffffffffffffff'

        bmx = list(bmatrix.contents)
        #self.bmx2D = np.array(bmx).reshape(self.natom, self.natom)
        #self.bmx1D = bmx

        self.bondneed = list(bondneed.contents)
        self.Lminstr = bool(Lminstr.contents)
        return self.Lminstr,bmx,list(bondneed.contents),list(surface.contents)


class AllStr(allstr):
    def readfile(self,inputfile, forcefile=False, allformat = 0):
        f= open(inputfile,'r')
        currentStr = -1
        for line in f:
            if ('Energy' in line\
                or 'React' in line\
                or ('TS' in line and 'sect' not in line)\
                or 'SSW-Crystl' in line\
                or 'Str' in line):
                self.append(SurfaceStr())
                currentStr +=  1
                self[currentStr].Lfor = False
                try:
                    self[currentStr].energy = float(line.split()[-1])
                    try :
                        self[currentStr].maxFF = float(item.split()[-2])
                    except :
                        self[currentStr].maxFF = 0

                    if self[currentStr].energy.is_integer():
                        self[currentStr].energy = float(line.split()[-2])
                except:
                    self[currentStr].energy = float(line.split()[-2])
                    self[currentStr].maxFF = 0
            elif 'DATE' in line:
                try:
                    self[currentStr].label = "%s" %reduce(lambda a,b:a+b , ["%s "%s   for s in line.split()[1:]])
                except:
                    self[currentStr].label = ''
            elif 'CORE' in line:
                self[currentStr].addatom(line,1 )
            elif ('PBC' in line )and ('ON' not in line):
                self[currentStr].abc= [float(x) for x in line.split()[1:]]
        f.close()
        for Str in self:
            Str.sortatombyele()
            Str.calAtomnum()
            Str.abc2lat()

        if forcefile:
            f = open(forcefile,'r')
            currentStr= -1
            for line in f:
                if "For" in line:
                    self[currentStr].Lfor = True
                    currentStr += 1
                    iatom = 0
                    for atom in self[currentStr].atom: atom.force = [0.0, 0.0, 0.0]
                elif len(line.split()) == 6:
                    self[currentStr].addStress(line)
                elif len(line.split()) == 3:
                    if "****" not in line: self[currentStr].addForce(line, iatom)
                    else:                  self[currentStr].addForce('0.0 0.0 0.0', iatom)
                    iatom += 1

        if allformat:
            for Str in self:
                Str.TransferToXYZcoordStr()

 
    def GetAllECFPname(self, numproc=24,savefp = False):
        _tmp = []
        for Str in self:
            _tmp.append((Str,savefp))

        pool = Pool(processes=numproc)
        result = pool.map_async(wrapGenECFPname, _tmp)
        pool.close(); pool.join()

        for Str,r in zip(self,result.get()):
            if savefp:
                Str.ECFPname = r[0]
                Str.FP = r[1]
            else:
                Str.ECFPname = r

    def GetTmpFakebmx(self,numproc=24,flag =2,colorflag =  1):
        if flag == 1:
            self.calAllBondMatrix(numproc=numproc)
        if flag == 2:
            self.calAllFakeBondMaxtrix(numproc=numproc)
        self.calAllSegmolecular (numproc=numproc)

    def AddGibbsEnergy(self):
    
        TS_CO2 = 298.15*213.795*0.001*0.01036
        TS_CO  = 298.15*197.653*0.001*0.01036
        TS_H2  = 298.15*130.680*0.001*0.01036
    
    
        print TS_CO2
        print TS_CO
        print TS_H2
    
        for Str in self:
            print Str.sminame
            if 'Carbon_Dioxide' or 'Hydrogen_Gas' in Str.sminame:
                nCO2 = 0
                nH2 = 0
                nCO = 0
                fragments = Str.sminame.split('+')
                for frag in fragments:
                    if 'Carbon_Dioxide' in frag and 'ads' not in frag:
                        nCO2 = nCO2 + 1
                    if 'Hydrogen_Gas' in frag and 'ads' not in frag:
                        nH2 = nH2 + 1
                    if 'Carbon_Monoxide' in frag and 'ads' not in frag:
                        nCO = nCO + 1
                Str.energy = Str.energy - (TS_CO2*nCO2)- (TS_H2*nH2) - (TS_CO * nCO)


    def GetAllsminame_fromreactstr_serial(self):
        for Str in self:
            Str.RemoveMetalBond()
            Str.group = Str.Segmolecular()
            substr = [[] for i in np.unique(Str.group)]
            for id,atom in enumerate(Str.atom):
                atom.id = id
                substr[Str.group[id]-1].append(atom)
            substr = sorted(substr, key=lambda x:calmass(x), reverse=True)
            _allmol = [sminame(sub, Str.bmx2D, Str.lat,Str.bondneed,2,Str.surfaceatom) for sub in substr]
            Str.sminame, strflag = glueSegStr(_allmol)
            Str.bmx2D= Str.bmxsave


    def GetAllsminame_fromreactstr(self,numproc=24,flag =2,colorflag =  1):
        for Str in self:
            Str.RemoveMetalBond()
        self.calAllSegmolecular (numproc=numproc)

        allgroup = []
        for Str in self:
            substr = [[] for i in np.unique(Str.group)]
            for id,atom in enumerate(Str.atom):
                atom.id = id
                substr[Str.group[id]-1].append(atom)
            substr = sorted(substr, key=lambda x:calmass(x), reverse=True)
            if flag == 1: allgroup.append((substr, Str.bmx2D, Str.lat,[],1,[]))
            if flag == 2: allgroup.append((substr, Str.bmx2D, Str.lat,Str.bondneed,2,Str.surfaceatom))
            #print str.bondneed

        pool = Pool(processes=numproc)
        result = pool.map_async(calAllName, allgroup)
        pool.close(); pool.join()

        for istr,(str,re) in enumerate(zip(self,result.get())):
            str.allmol = re
            str.id     = istr

        for str in self:
            if colorflag:
                str.sminame, strflag = glueSegStr(str.allmol)
            else:
                str.sminame, strflag = glueSegStr_pure(str.allmol)

        for Str in self:
            Str.bmx2D = Str.bmxsave

class singlefp(object):
    def __init__(self):
        self.index = 0
        self.neighbour= []
        self.allatom = []
        self.trace = []
        self.save = 1
        self.coreele = 0
        self.coreid = 0

    def Inherit(self,parent):
        self.index = parent.index
        self.allatom = parent.allatom[:]
        self.trace = parent.trace[:]
        self.coreele = parent.coreele
        self.coreid = parent.coreid

    def genidentifier(self,dim=10000000):
        self.neighbour.sort(key= lambda x:x[0])
        
        #if self.coreele == 29:
        #    print len(self.neighbour)

        for nb in self.neighbour:
            self.trace.append(nb[1])

        self.allatom.sort()
        self.allatomset = np.unique(self.allatom)

        tmpall= [(1,self.index)]
        for nb in self.neighbour:
            tmpall.append(nb[0])
            #print nb[0]
            #tmpall.append(nb[0][0])
            #tmpall.append(nb[0][1])
        self.index = hashlist(tmpall,dim)

        #string = '1:%s'(%self.index)
        #for item in self.neighbour:
        #    string = string+ str(item)
        


class ECFP(object):
    def __init__(self,structure, nlayer, dim=10000000):
        self.Str = structure
        self.nlayer = nlayer
        self.allfp = []
        self.dim = dim

#    def PreparebyLayer(self):
#        for atom in self.str.atom:
#            atom.layeratom = {}
#            atom.layeratom[0] = [atom]
#
#        neighbourdict = {}
#        for ilayer in range(nlayer):
#            atom.layeratom[ilayer+1] = []
#            for iatom,atom in enumerate(self.str.atom):
#                #neighbourdict[(iatom,ilayer)] = []
#                for subatom in atom.layeratom[ilayer]:
#                    #neighbourdict[(iatom,ilayer)].append(subatom.id)
#                    for inb,newatom in enumerate(subatom.neighbour):
#                        atom.layeratom[ilayer+1].append(newatom)
#                        self.cyclenb[(iatom,ilayer+1)].append([subatom.bond[inb],newatom.id])

    def GenECFP(self):
        self.InitIdentifier()
        for ilayer in xrange(1,self.nlayer):
            self.UpdateLayer(ilayer)
        allindex = []
        for fp in self.allfp:
            allindex.append(fp.index)
        allindex.sort()
        self.allindex  = np.unique(allindex)
        #self.allindex,self.allindexconut =list(np.unique(allindex,return_counts=True))
        self.GenDict()
        return

    def GenDict(self):
        #AtomsettoIndex = {}
        #IndextoAtomset = {}
        IndextoFp = {}
        for fp in self.allfp:
            #AtomsettoIndex[set(fp.allatomset)] = fp.index
            #IndextoAtomset[fp.index] = fp.allatomset
            if fp.index not in IndextoFp.keys():
                IndextoFp[fp.index] = [fp]
            else:
                IndextoFp[fp.index].append(fp)

        #self.IndextoAtomset = IndextoAtomset
        self.IndextoFp = IndextoFp


    def InitIdentifier(self):
        self.Str.GenFeature()
        self.genbylayer=[[0 for x in range(self.Str.natom)]]
        for iatom,atom in enumerate(self.Str.atom):
            newfp= singlefp()
            newfp.allatom.append(atom.id)
            newfp.trace.append(atom.id)
            newfp.coreele = atom.ele
            newfp.coreid = iatom
            #print atom.feature
            newfp.index =hashlist(atom.feature,self.dim)
            #print newfp.index
            #print newfp.index
            newfp.allatom.sort()
            newfp.allatomset = np.unique(newfp.allatom)
            #newfplist.append(newfp)
            self.genbylayer[0][iatom] = newfp
            #if newfp.coreele < 18: continue
            #if newfp.coreele == 1: continue
            if  (newfp.coreele< 18):
                self.allfp.append(newfp)


    def UpdateLayer(self,layer):
        newfplist =[]
        for iatom,atom in enumerate(self.Str.atom):
            newfp = singlefp()
            #newfp.append((1,self.genbylayer[layer-1][atom.id]))
            newfp.Inherit(self.genbylayer[layer-1][atom.id])
            for subatomid in self.genbylayer[layer-1][atom.id].allatom:
                #if self.Str.atom[subatomid].ele > 18: continue
                for nb in self.Str.atom[subatomid].nb:
                    if nb[1] not in newfp.allatom :# and self.str.atom[nb[1]].ele != 1:
                        newfp.neighbour.append([(nb[0], self.genbylayer[0][nb[1]].index),nb[1]])
                        newfp.allatom.append(nb[1])
            
            #newfp.allatomset = set(newfp.allatom)
            newfp.genidentifier(self.dim)
            newfplist.append(newfp)
        #for fp in  newfplist:
            #print fp.save, fp.index
        newfplist =self.CheckDuplicate(newfplist)
        #print 'CheckDuplicate'
        #for fp in  newfplist:
        #    print fp.save, fp.index

        self.genbylayer.append([0 for x in range(self.Str.natom)])
        for ifp,fp in enumerate(newfplist):
            if fp.save ==1:
                self.genbylayer[layer][ifp]= fp
                #if fp.coreele < 18: continue
                #if fp.coreele == 1: continue
                if fp.coreele< 18:
                    self.allfp.append(fp)
            else:
                self.genbylayer[layer][ifp] = self.genbylayer[layer-1][ifp]
        return

    def CheckDuplicate(self,newfplist):
        for newfp in newfplist:
            for fp in self.allfp:
                if set(newfp.allatomset) == set(fp.allatomset):
                    # for unity
                    #newfp.index = fp.index
                    # old mode
                    newfp.save = 0
                    break

        for newfp1 in newfplist:
            for newfp2 in newfplist:
                if set(newfp1.allatomset) == set(newfp2.allatomset):
                    if newfp1.index > newfp2.index:
                        #newfp1.index = newfp2.index
                        newfp1.save = 0
                    elif newfp1.index < newfp2.index:
                        #newfp2.index = newfp1.index
                        newfp2.save = 0
        return newfplist


#class NNFingerPrint(object):
#    def __init__(self):
#        return
#
#    def neighbouractive(self):
#        for degree in degrees:
#            atom_neighbour= self.atom[j].neighbour
#            bond_neighbour = self.bond[j].neighbour
#
#
#    def layerupdate(self):
#        return

def hashlist(list,dim =10000000):
    #print list
    string = ''
    for i in list:
        string = string +str(i)

    #print string
    md5 =hashlib.md5()

    #print md5.digest_size
    #print md5.block_size
    md5.update(string)
    #print md5.digest
    #print string

    #return md5.hexdigest()
    _tmp= md5.hexdigest()
    #print int((int(_tmp,16))%dim)
    return int((int(_tmp,16))%dim)
    #return zlib.adler32(string)


if (__name__ == "__main__"):
    test =[1,2332421,1,14757242,2,124253533]
#    string = ''
#    for i in test:
#        string = string +str(i) 
#    sha1 =hashlib.sha1()
#    sha1.update(string)
#    print sha1.hexdigest()
#
#
#
#    md5 =hashlib.md5()
#    md5.update(string)
#    print md5.hexdigest()
#    
#
#
#    print zlib.adler32(string)
    print (hashlist(test))


    test =AllStr()
    test.readfile('target.arc')
    test.calAllFakeBondMaxtrix()
    fp= ECFP(test[0],3)
    fp.GenECFP()
    #print len(fp.allfp)

    #for sfp in fp.allfp:
    #    print sfp.index

    print (fp.allindex)
    print (fp.IndextoAtomset)
    for ilayer in range(3):
        for iatom in range(test[0].natom):
            print (ilayer,iatom,fp.genbylayer[ilayer][iatom].index, fp.genbylayer[ilayer][iatom].trace)
#    all.GetAllpuresminame()
#    all.GenCanoicalSmiles(chiral= False)
#
#    ndim = 1000000
#    info= {}
#    info2 = {}
#    for i in range(len(test)):
#        str1 = test[i]
#        m1 = Chem.MolFromSmiles(str1.canname)
#        Chem.SanitizeMol(m1)
#        fptest= AllChem.GetHashedMorganFingerprint(m1,2,nBits=ndim, bitInfo = info,useChirality=False)
#        if i ==0:
#            fpall =fptest
#        fpall =fpall + fptest
#    print info
