import numpy as np
import openbabel
import pybel
from bondOrder import Bond
from multiprocessing import Pool
from Font import outPink, outBlue, outSky, outGray, endMark
import ctypes
from ctypes import pointer
from PeriodicTable import Elemass,Eletable
from copy import deepcopy
def puresminame(str, bmatrix, lat, bondneed=[], modeflag=1, surface =[]):
    surfaceflag = '-'

    if modeflag == 1:

        iza = [atom.ele for atom in str]
        sna = len(str)**2
        latinv    = np.linalg.inv(lat)

        bmx = []
        for atom in str:
            bmx.append([bmatrix[atom.id][jatm.id] for jatm in str])
        bmx = np.array(bmx).reshape(sna)

        radical = [0 for atom in str]
        na = len(str)
        fa = reduce(lambda a,b:a+b, [list(np.matmul(atom.xyz,latinv)) for atom in str])

        c_na  = pointer(ctypes.c_int(na))
        c_iza = pointer((ctypes.c_int*na)(*iza))
        c_rad = pointer((ctypes.c_int*na)(*radical))
        c_bmx = pointer((ctypes.c_int*(sna))(*bmx))
        c_fa  = pointer((ctypes.c_double*(3*na))(*fa))

        c_lat = pointer((ctypes.c_double*9)(*reduce(lambda a,b:a+b,lat)))
        #print 'test1'
        program = ctypes.cdll.LoadLibrary('/home10/kpl/pymodule/script/Lib_bondmatrix/bondmatrix.so')
        #print 'test2'
        program.hradicalorder_(c_na, c_iza, c_bmx, c_rad, c_fa, c_lat)
        #print np.array(c_rad.contents)
        if np.array(c_rad.contents).max() > 0:
            flag = False
            radflag = '-'
            for iatom,irad in enumerate(np.array(c_rad.contents)) :
                if irad > 0:
                    radflag= radflag+'%s%d'%(Eletable[iza[iatom]-1],c_rad.contents[iatom])
        else:
            flag = True
    if modeflag ==2 :
        radflag = '-'
        flag = True
        for atom in str:
            if bondneed[atom.id]!= 0:
                flag = False
                radflag= radflag+'%s%d'%(Eletable[atom.ele-1],bondneed[atom.id])

        for atom in str:
            if surface[atom.id] == 1:
                surfaceflag =surfaceflag+'%sads'%(Eletable[atom.ele-1])


    mol = openbabel.OBMol()
    for i,atom in enumerate(str):
        newatom = mol.NewAtom()
        newatom.SetAtomicNum(atom.ele)
        newatom.SetVector(float(atom.xyz[0]), float(atom.xyz[1]), float(atom.xyz[2]))
        for j,jatm in enumerate(str):
            od = bmatrix[atom.id][jatm.id]
            if od > 0: mol.AddBond(i+1,j+1,od)
    newmol = pybel.Molecule(mol)

    molecular = newmol.write("can").strip()
    if molecular == 'C=O':
        if len(str) == 2:
            molecular = '[C]=O'; flag = True

    #if not flag:
    #    molecular = molecular+radflag
    #if surfaceflag != '-':
    #    molecular = molecular+surfaceflag

    for i,atom in enumerate(str):
        if atom.ele > 18:
            printmolecular =''
            return printmolecular, flag
    return molecular, flag


def sminame (str, bmatrix, lat, bondneed=[], modeflag=1, surface =[]):
    # judge whether a well molecular
    surfaceflag = '-'

    if modeflag == 1:

        iza = [atom.ele for atom in str]
        sna = len(str)**2
        latinv    = np.linalg.inv(lat)

        bmx = []
        for atom in str:
            bmx.append([bmatrix[atom.id][jatm.id] for jatm in str])
        bmx = np.array(bmx).reshape(sna)

        radical = [0 for atom in str]
        na = len(str)
        fa = reduce(lambda a,b:a+b, [list(np.matmul(atom.xyz,latinv)) for atom in str])

        c_na  = pointer(ctypes.c_int(na))
        c_iza = pointer((ctypes.c_int*na)(*iza))
        c_rad = pointer((ctypes.c_int*na)(*radical))
        c_bmx = pointer((ctypes.c_int*(sna))(*bmx))
        c_fa  = pointer((ctypes.c_double*(3*na))(*fa))

        c_lat = pointer((ctypes.c_double*9)(*reduce(lambda a,b:a+b,lat)))
        #print 'test1'
        program = ctypes.cdll.LoadLibrary('/home10/kpl/pymodule/script/Lib_bondmatrix/bondmatrix.so')
        #print 'test2'
        program.hradicalorder_(c_na, c_iza, c_bmx, c_rad, c_fa, c_lat)
        #print np.array(c_rad.contents)
        if np.array(c_rad.contents).max() > 0:
            flag = False
            radflag = '-'
            for iatom,irad in enumerate(np.array(c_rad.contents)) :
                if irad > 0:
                    radflag= radflag+'%s%d'%(Eletable[iza[iatom]-1],c_rad.contents[iatom])
        else:
            flag = True

    if modeflag ==2 :
        radflag = '-'
        flag = True
        for atom in str:
            if bondneed[atom.id]!= 0:
                flag = False
                radflag= radflag+'%s%d'%(Eletable[atom.ele-1],bondneed[atom.id])

        for atom in str:
            if surface[atom.id] == 1:
                surfaceflag =surfaceflag+'%sads'%(Eletable[atom.ele-1])


    mol = openbabel.OBMol()
    for i,atom in enumerate(str):
        newatom = mol.NewAtom()
        newatom.SetAtomicNum(atom.ele)
        newatom.SetVector(float(atom.xyz[0]), float(atom.xyz[1]), float(atom.xyz[2]))
        for j,jatm in enumerate(str):
            od = bmatrix[atom.id][jatm.id]
            if od > 0: mol.AddBond(i+1,j+1,od)
    newmol = pybel.Molecule(mol)

    molecular = newmol.write("smi").strip()
    if   molecular == '[H][H]':  molecular = 'Hydrogen_Gas'
    elif molecular == 'C(=O)=O':
        molecular = 'Carbon_Dioxide'
    elif molecular == 'C(=O)O': 
        if len(str) == 5:
            molecular = 'Formic_Acid'
        elif len(str) == 4:
            molecular = 'Formate/Carboxyl'
    elif molecular == 'C':
        if len(str) == 5:
            molecular = 'Methane'; flag = True
        elif len(str) == 4:
            molecular = 'Methyl';# flag = False
    elif molecular == 'O':
        if len(str) == 3:
            molecular = 'Water'; flag = True
        elif len(str) == 2:
            molecular = 'Hydroxy';# flag = False
    elif molecular == 'C#[O]':
        molecular = 'Carbon_Monoxide'; flag = True
    elif molecular == 'C(=O)(O)O' and len(str)==6:
        molecular = 'Carbonic-Acid'
    elif molecular == 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O':
        molecular = 'alpha_D_glucopyranose'
    elif molecular == '[C@@H]1([C@@H]([C@H]([C@@H]([C@@H](CO)O1)O)O)O)O':
        molecular = 'beta_D_glucopyranose'
    elif molecular == 'C(=O)[C@@H]([C@H]([C@@H]([C@@H](CO)O)O)O)O' :
        molecular = 'P1 : D-glucose'
    elif molecular == 'C(C(=O)[C@H]([C@@H]([C@@H](CO)O)O)O)O' :
        molecular = 'D-fructose'
    elif '[C@H]12[C@@H]([C@H]([C@@H]([C@@H](CO2)O1)O)O)O' in molecular:
    #elif molecular == '[C@H]12[C@@H]([C@H]([C@@H]([C@@H](CO2)O1)O)O)O + Water'        :
        molecular = 'P2 : Levoglucosan'
    elif '[C@@H]1([C@@H]([C@@H]2[C@@H]([C@@H](CO2)O1)O)O)O' in molecular:
        molecular = 'P6 : Levoglucosan-2'
    elif molecular == '[C@@H]1([C@@H]([C@H]([C@@H]([C@@H](CO)O)O1)O)O)O' :
        molecular = 'P3 :  beta-D-glucofuranose '
    elif molecular == '[C@@H]1([C@@H]([C@H]([C@@H]([C@@H](CO1)O)O)O)O)O' :
        molecular = 'P4 : beta-D-glucoseptanose'
    elif 'C(=O)[C@@H]1[C@H]([C@@H]([C@@H](CO)O1)O)O' in molecular:
        molecular = 'P5 : 2,5-anhydro-D-mnnose'

    elif 'C(=O)[C@H]1[C@H]([C@@H]([C@@H](CO)O1)O)O' in molecular:
        molecular = 'P5 : 2,5-anhydro-D-mnnose-2'
    #elif molecular == '' :
    #    molecular =
    elif molecular == 'CC':
        molecular = "Ethane"
        if len(str) == 7: molecular = 'Ethyl'
    elif molecular == 'O=O': molecular = 'Oxygen_Gas'
    elif molecular == 'N#N': molecular = 'Nitrogen_Gas'
    elif molecular == 'C=C':
        molecular = "Ethene"
        if len(str) == 5: molecular = 'Ethenyl'
    elif molecular == 'CO': 
        if len(str) ==6:
            molecular = "Methanol"
    elif molecular == 'CCO' or molecular == 'C(C)O' : molecular = "Ethanol"
    elif molecular == 'CCC' or molecular == 'C(C)C' :
        molecular = "Propane"
        if len(str) == 10: molecular = 'Propyl'
    elif molecular == 'C=O':
        if len(str) == 2:
            molecular = 'Carbon_Monoxide'; flag = True
        elif len(str) == 4:
            molecular = 'Formaldehyde'; flag = True
    elif molecular == '[H]':
        flag =True
    elif molecular == '[OH3]':
        flag = True


    if not flag:
        molecular = molecular+radflag
    if surfaceflag != '-':
        molecular = molecular+surfaceflag

    for i,atom in enumerate(str):
        if atom.ele > 18:
            printmolecular =''
            return printmolecular, flag
    return molecular, flag

def calAllName (data):
    group, bondmatrix, lat, bondneed,flag,surface = data
    return [sminame(sub, bondmatrix, lat,bondneed,flag,surface) for sub in group]

def calAllpureName (data):
    group, bondmatrix, lat, bondneed,flag,surface = data
    return [puresminame(sub, bondmatrix, lat,bondneed,flag,surface) for sub in group]
    #return [puresminame(group, bondmatrix, lat,bondneed,flag,surface)]

def calmass (str):
    allmass = np.array([Elemass[atm.ele-1] for atm in str])
    return allmass.sum()

#def judgeReaction (strin1, strin2):

def singleFindName (strin):
    str = deepcopy(strin)
    bmx   = str.Bondmatrix()
    bmx2D = np.array(bmx).reshape(str.natm, str.natm)
    str.bmx1D = bmx
    group = str.Segmolecular(recal=False)
    substr = [[] for i in np.unique(group)]
    for id,atom in enumerate(str.atom):
        atom.id = id; substr[group[id]-1].append(atom)
    substr = sorted(substr, key=lambda x:calmass(x), reverse=True)
    return [sminame(sub, bmx2D, str.lat) for sub in substr], bmx2D



def glueSegStr (allmol):
    outstr,cellflag = "", True
    p =0
    for i,(name,flag) in enumerate(allmol):
        if flag: font = outSky
        else:    font = outBlue
        if name == '' : continue
        if p>0:    outstr += "+%s%s%s"%(font, name, endMark)
        elif p==0:
            outstr += "%s%s%s"%(font, name, endMark)
            p=p+1
        if not flag: cellflag = False
    return outstr.strip(), cellflag

def glueSegStr_pure (allmol):
    outstr,cellflag = "", True
    p =0
    for i,(name,flag) in enumerate(allmol):
        if flag: font = outSky
        else:    font = outBlue
        if name == '' : continue
        if p>0:    outstr += ".%s"%(name)
        elif p==0:
            outstr += "%s"%(name)
            p=p+1
        if not flag: cellflag = False
    return outstr.strip(), cellflag

