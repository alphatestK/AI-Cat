from structure_new import Str
import re
import numpy as np
from reactpattern_pure import ReactPatternNew
from ECFP0 import ECFP,AllStr
from reactstr import BmxStr

class Pair(AllStr):
    def __init__(self):
        list.__init__(self)

    def AddStr(self,str1,str2,TS=0,TSstr=False,sortflag = True):
        self.append(str1)
        self.append(str2)
        #self.append(TSstr)
        self.TS= TS
        #if sortflag: self.sort(key= lambda X: X.name)
        #self.name = '%s'%(self[0].name+self[1].name)
        self.TSstr =TSstr
        return

    def GenReactPattern(self,rcdepth = 1):

        r = self.GetReactCenter(rcdepth)
        if isinstance(r,int):
            print 'error in extract react center'
            return -2
        else:
            _tmpstr1,_tmpstr2,groupdict1,groupdict2,groupdict01,groupdict02,atomlist0,atomlist = r


        #print groupdict1
        #print groupdict01
        FP =ECFP(self[0],3)
        FP.GenECFP()
        self[0].ECFP= FP
        #print groupdict1
        self.GetSubstrECFP_list(FP,groupdict1,groupdict01)
        self.rp.BmxSort(self[0].bmx2D,self[1].bmx2D)
        rp1 = self.rp


        FP =ECFP(self[1],3)
        FP.GenECFP()
        self[1].ECFP= FP
        #print groupdict2
        self.GetSubstrECFP_list(FP,groupdict2,groupdict02)
        self.rp.BmxSort(self[1].bmx2D,self[0].bmx2D)
        rp2 = self.rp

        rp1.GenName0(self[0].ECFP,self[1].ECFP,atomlist0)
        rp2.GenName0(self[1].ECFP,self[0].ECFP,atomlist0)

        return  rp1,rp2,atomlist0,atomlist

    def GetReactCenter(self, depth =1, Ltype = False , LaddfakeZPE = False):
        # must have caled bmx

        str1 = self[0]
        str2 = self[1]

        str1.GetAtomInfo()
        str2.GetAtomInfo()

        bmxc = str1.bmx2D - str2.bmx2D
        atomlist = []


        iformation = 0
        ibreak = 0
        for i in range(str1.natom):
            for j in xrange(i,str1.natom):
                if bmxc[i][j] != 0:
                    if ((str1.atom[i].ele > 18) or (str1.atom[j].ele > 18)): continue
                    if i not in atomlist:  atomlist.append(i)
                    if j not in atomlist:  atomlist.append(j)
                    if (str1.atom[i].ele == 1 and str2.atom[j].ele ==6) or (str1.atom[i].ele == 1 and str2.atom[j].ele ==8):
                        if bmxc[i][j] == -1:
                            iformation = iformation + 1
                        elif bmxc[i][j] == 1:
                            ibreak = ibreak + 1

        if (ibreak >= iformation):
            ZPE_TS = -0.15 *(ibreak-iformation)
            ZPE_heat = -0.15 *(ibreak-iformation)
        else:
            ZPE_TS = 0
            ZPE_heat = 0.15*(iformation-ibreak)

        if len(atomlist) ==0 :
            return -2

        atomlist0 = atomlist[:]
        self.nrc0 = len(atomlist0)

        if Ltype:
            if (self.nrc0 == 2):
                i = atomlist0[0]
                j = atomlist0[1]
                if str1.atom[i].ele == 6 and str1.atom[j].ele == 6 and bmxc[i][j] == -1 :
                    reacttype = "C-C formation"
                elif str1.atom[i].ele == 6 and str1.atom[j].ele == 6 and bmxc[i][j] == 1 :
                    reacttype = "C-C break"
                elif str1.atom[i].ele == 6 and str1.atom[j].ele == 8 and bmxc[i][j] == 1 :
                    reacttype = "C-O break"
                elif str1.atom[i].ele == 6 and str1.atom[j].ele == 8 and bmxc[i][j] == -1 :
                    reacttype = "C-O formation"
                elif str1.atom[i].ele == 8 and str1.atom[j].ele == 8 and bmxc[i][j] == 1 :
                    reacttype = "O-O break"
                elif str1.atom[i].ele == 8 and str1.atom[j].ele == 8 and bmxc[i][j] == -1 :
                    reacttype = "O-O formation"
                elif str1.atom[i].ele == 1 and str1.atom[j].ele == 6 and bmxc[i][j] == -1 :
                    reacttype = "C-H formation"
                elif str1.atom[i].ele == 1 and str1.atom[j].ele == 6 and bmxc[i][j] == 1 :
                    reacttype = "C-H break"
                elif str1.atom[i].ele == 1 and str1.atom[j].ele == 1 and bmxc[i][j] == -1 :
                    reacttype = "H-H formation"
                elif str1.atom[i].ele == 1 and str1.atom[j].ele == 1 and bmxc[i][j] == 1 :
                    reacttype = "H-H break"
                elif str1.atom[i].ele == 1 and str1.atom[j].ele == 8 and bmxc[i][j] == -1 :
                    reacttype = "O-H formation"
                elif str1.atom[i].ele == 1 and str1.atom[j].ele == 8 and bmxc[i][j] == 1 :
                    reacttype = "O-H break"
                else:
                    reacttype = "unknown"
            else:
                reacttype = "unknown"

        if depth ==1:
            for i in atomlist0:
                for j in range(str1.natom):
                    if str1.bmx2D[i][j] >0  and j not in atomlist :  atomlist.append(j)
                    if str2.bmx2D[i][j] >0  and j not in atomlist :  atomlist.append(j)

        #atomlist.sort()

        _tmpstr1 = Str()
        _tmpstr2 = Str()
        rcatomindice= {}

        for i in atomlist:
            if str1.atom[i].ele > 18: continue
            _tmpstr1.atom.append(str1.atom[i])
            _tmpstr2.atom.append(str2.atom[i])
            rcatomindice[len(_tmpstr1.atom)-1]=i

        _tmpstr1.abc = str1.abc[:]
        _tmpstr2.abc = str2.abc[:]

        _tmpstr1.calAtomnum()
        _tmpstr1.abc2lat()
        _tmpstr2.calAtomnum()
        _tmpstr2.abc2lat()

        _tmpstr1.bmx2D=np.array([[0 for i in range(_tmpstr1.natom)] for i in range(_tmpstr1.natom)])
        _tmpstr2.bmx2D=np.array([[0 for i in range(_tmpstr1.natom)] for i in range(_tmpstr1.natom)])

        for i in range(_tmpstr1.natom):
            for j in range(_tmpstr1.natom):
                _tmpstr1.bmx2D[i][j]= str1.bmx2D[rcatomindice[i]][rcatomindice[j]]
                _tmpstr2.bmx2D[i][j]= str2.bmx2D[rcatomindice[i]][rcatomindice[j]]



        reactgroup1=[]
        groupdict1 = {}
        groupdict01 = {}
        for i,atom in enumerate(_tmpstr1.atom):
            atom.expbond = str1.atom[rcatomindice[i]].expbond
            atom.imph =  str1.atom[rcatomindice[i]].imph
            atom.group = str1.group[rcatomindice[i]]
            reactgroup1.append(atom.group)
            if atom.group not in groupdict1.keys():
                groupdict1[atom.group] =[atom.id]
                if atom.id in atomlist0:
                    groupdict01[atom.group] =[atom.id]

                #groupdict1[atom.group].append(atom.id)
            else:
                groupdict1[atom.group].append(atom.id)
                if atom.id in atomlist0:
                    groupdict01[atom.group].append(atom.id)

        reactgroup1.sort()
        reactgroup1 = np.unique(reactgroup1)

        reactgroup2 = []
        groupdict2 = {}
        groupdict02 = {}
        for i,atom in enumerate(_tmpstr2.atom):
            atom.expbond = str2.atom[rcatomindice[i]].expbond
            atom.imph =  str2.atom[rcatomindice[i]].imph
            atom.group = str2.group[rcatomindice[i]]
            reactgroup2.append(atom.group)
            if atom.group not in groupdict2.keys():
                groupdict2[atom.group] =[atom.id]

                if atom.id in atomlist0:
                    groupdict02[atom.group] =[atom.id]

                #groupdict1[atom.group].append(atom.id)
            else:
                groupdict2[atom.group].append(atom.id)
                if atom.id in atomlist0:
                    groupdict02[atom.group].append(atom.id)
        reactgroup2.sort()
        reactgroup2 = np.unique(reactgroup2)

#        reactcenter = AllStr()
#        reactcenter.append(_tmpstr1)
#        reactcenter.append(_tmpstr2)
#
#        reactcenter.printall('reactcenter.arc')
#

#        if len(reactgroup1) >2 or len(reactgroup2) > 2:
#            return  0
#        else:
#            return _tmpstr1,_tmpstr2,groupdict1,groupdict2

        if (LaddfakeZPE):
            return _tmpstr1,_tmpstr2,groupdict1,groupdict2,groupdict01, groupdict02,atomlist0,atomlist, reacttype,ZPE_TS,ZPE_heat

        if (not Ltype):
            return _tmpstr1,_tmpstr2,groupdict1,groupdict2,groupdict01, groupdict02,atomlist0,atomlist
        else:
            return _tmpstr1,_tmpstr2,groupdict1,groupdict2,groupdict01, groupdict02,atomlist0,atomlist,reacttype

#        reactcenter = AllStr()
#        reactcenter.append(_tmpstr1)
#        reactcenter.append(_tmpstr2)
#
#        reactcenter.OutMolFile(printlist= [0], outfile = 'rc1.mol')
#        reactcenter.OutMolFile(printlist= [1], outfile = 'rc2.mol')
#
#        reactcenter.printall('reactcenter.arc')
#
#        return reactcenter
#

    def GetSubstrECFP_list(self,ECFP,groupdict1,groupdict0):
        self.rp = ReactPatternNew()
        for key,value in groupdict0.items():
            self.rp.tmpatomall = []
            startfp =ECFP.genbylayer[1][value[0]]
            atomall = groupdict1[key]

            self.GetSubstrECFPlink(ECFP,value,atomall,[],startfp)

            for atomid in value:
                if set(self.rp.tmpatomall) == set(atomall): break
                if atomid not in self.rp.tmpatomall:
                    startfp =ECFP.genbylayer[1][atomid]
                    self.GetSubstrECFPlink(ECFP,value,atomall,[],startfp)
        #print self.rp.atomall
        #print self.rp.fpmatch
        return


    def GetSubstrECFPlink(self,ECFP,atomlist0,atomall,nowsite,fp):
        #print atomall
        #print atomlist0
        if set(self.rp.tmpatomall)==set(atomall):
            return

        fplist = []
        nowsitel = []
        nextfp = []
        #print fp.coreid,fp.trace
        for iall,atomid in enumerate(fp.trace):
            if atomid in atomlist0 and atomid not in self.rp.contents:
                fp = ECFP.genbylayer[1][atomid]
                self.rp.contents.append(atomid)
                _tmpnowsite = nowsite[:]
                _tmpnowsite.append(iall)
                nowsitel.append(_tmpnowsite)
                self.rp.fpmatch.append((_tmpnowsite,fp.index))
                nextfp.append(fp)
            if atomid in atomall and atomid not in self.rp.tmpatomall:
                self.rp.atomall.append(atomid)
                self.rp.tmpatomall.append(atomid)
                _tmpnowsite = nowsite[:]
                _tmpnowsite.append(iall)
                self.rp.atommatch.append((_tmpnowsite,atomid))


        if set(self.rp.tmpatomall) == set(atomall):
            return

        if len(nextfp) > 0:
            ECFPl= [ECFP]*len(nextfp)
            atomlistl= [atomlist0]*len(nextfp)
            atomalll= [atomall]*len(nextfp)
            #atomnowl= [atomnow]*len(nextfp)
            return map(self.GetSubstrECFPlink,ECFPl,atomlistl,atomalll,nowsitel,nextfp)
        else:
            return


    def CheckReactPattern_new(self,str,rp , stdout,depth =3):
        #rp1, rp2 =self.GenSingleReactPattern()

        plist = self.RunPattern(rp,str,depth=depth)#,CalName =True)

        #for str in plist:
        #    print str.sminame
        rightp= []
        for product in plist:
            lmatch =self.StrCheckSame(product,stdout,depth)
            if lmatch == 1:
                rightp.append(plist)
        if len(rightp) >0:
            return 1
        else:
            print 'wrong product'
            return 0


    def RunPattern(self,rp,reactstr,depth=3, CalName = False):
        FP = ECFP(reactstr,3)
        FP.GenECFP()

        matchcombine =rp.GetMatchFp(FP)
        #print matchcombine
        #print rp.atomall

        newstrlist = AllStr()

        if matchcombine == -1:
            print 'fail to match fp'
            return newstrlist

        for matchfp in matchcombine:
            r= self.GenNewStr_bmx(rp,reactstr,FP,matchfp)
            if r == -1: continue
            newstrlist.append(r)
        if (CalName):
            newstrlist.GetAllsminame_fromreactstr()

        #print newstrlist[0].sminame
        return newstrlist


    def StrCheckSame(self,product, stdout, depth = 3):
#        for i in range(product.natom):
#            for j in range(product.natom):
#                if product.bmx2D[i][j] != stdout.bmx2D[i][j]:
#                    print i,j ,product.bmx2D[i][j],  stdout.bmx2D[i][j]


        fp1 =ECFP(product,depth)
        fp1.GenECFP()

        fp2 =ECFP(stdout,depth)
        fp2.GenECFP()

        lmatch  = 1
        for i,index in enumerate(fp2.allindex):
            if fp1.allindex[i] != index:
                #print fp1.allindex
                #print fp2.allindex
                #for ilayer in range(3):
                #    for iatom in range(stdout.natom):
                #        print ilayer,iatom,fp1.genbylayer[ilayer][iatom].index, fp1.genbylayer[ilayer][iatom].trace
                #        print ilayer,iatom,fp2.genbylayer[ilayer][iatom].index, fp2.genbylayer[ilayer][iatom].trace


                lmatch = 0
                break


        return lmatch

    def GenNewStr_bmx(self,rp,reactstr,FP,matchcombine):

        newstr = BmxStr()
        newstr.Inherit(reactstr)
        isect = 0

        #for fp in matchcombine:
        #    print fp.coreid

        trace = []
        for atominfo in rp.atommatch:
            if len(atominfo[0]) == 1 and (atominfo[0][0]== 0):
                isect = isect+1
            atomid = rp.GenMatch(matchcombine[isect-1].coreid,FP,atominfo[0])
            if atomid in trace:
                #print 'wrong match in GenNewStr'
                return -1
            trace.append(atomid)
            #print trace
        #print 'trace',trace
        #print 'trace',trace


        for i,reali in enumerate(trace):
            for j,realj in enumerate(trace):
                if newstr.bmx2D[reali][realj] != rp.bmxIS[i][j]:
                    print trace
                    print i,j,reali,realj,  newstr.bmx2D[reali][realj] ,rp.bmxIS[i][j]
                    print 'wrong match in bmxIS'

                    print rp.bmxIS
                    print newstr.bmx2D

                    #print 'break',fe
                    return -1
                newstr.bmx2D[reali][realj]= rp.bmxFS[i][j]
        return  newstr

    def OutProduct(self,reactstr,rp,depth =3):
        #plist =self.RunReactPattern(reactstr,rp,depth)
        plist =self.RunPattern(rp,reactstr,depth)
        #plist.GetAllsminame_fromreactstr()

        finalp= []
        allname = []
        for Str in plist:
            Str.GenECFPname()
            if Str.ECFPname not in allname :
                finalp.append(Str)
                allname.append(Str.ECFPname)
            #print Str.ECFPname
        #for str in plist:
        #    print str.sminame

        return finalp

