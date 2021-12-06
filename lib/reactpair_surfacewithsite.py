import os
from ECFP_withsite import ECFP
from ECFP_withsite import SurfaceStr 
from ECFP_withsite import AllStr 
import numpy as np
from reactstr import BmxStr

class ReactionPair(AllStr):
    def __init__(self):
        list.__init__(self)

    def AddStr(self,str1,str2,TS= 0,TSstr=False,sortflag = True):
        self.append(str1)
        self.append(str2)
        #self.append(TSstr)
        self.TS= TS
        self.TSstr =TSstr

        #if sortflag: self.sort(key= lambda X: X.name)
        #self.name = '%s'%(self[0].name+self[1].name)
        return

    def GenECFP_withsite(self,rclist0):

        self[0].RemoveUnrelaventBond(rclist0)
        self[1].RemoveUnrelaventBond(rclist0)

        self[0].GenECFPname(depth=3,dim=10000000,savefp =True)
        self[1].GenECFPname(depth=3,dim=10000000,savefp =True)

        self.name = '%s'%(self[0].ECFPname+self[1].ECFPname)


    def SitePatchtoRP(self, rp1,rp2,rclist0):
        self.rp1 =rp1
        self.rp2 =rp2

        nmetal,metalatom,metalbmx2D1,metalbmx2D2, surfacemetal1, surfacemetal2, supercellinfo1, supercellinfo2 = self.SaveSurfaceMetalBmx()
        self.rp1.nmetal = nmetal
        self.rp1.ISmetalbmx2D= metalbmx2D1
        self.rp1.FSmetalbmx2D= metalbmx2D2
        self.rp1.ISsurfacemetal = surfacemetal1
        self.rp1.FSsurfacemetal = surfacemetal2
        self.rp1.ISsupercellinfo = supercellinfo1
        self.rp1.FSsupercellinfo = supercellinfo2


        self.rp2.nmetal = nmetal
        self.rp2.ISmetalbmx2D= metalbmx2D2
        self.rp2.FSmetalbmx2D= metalbmx2D1
        self.rp2.ISsurfacemetal = surfacemetal2
        self.rp2.FSsurfacemetal = surfacemetal1
        self.rp2.ISsupercellinfo = supercellinfo2
        self.rp2.FSsupercellinfo = supercellinfo1

        self.rp1.ISsitematch = self.AddSiteInfo(self[0],self.rp1,metalatom)
        self.rp1.FSsitematch = self.AddSiteInfo(self[1],self.rp1,metalatom)
        self.rp2.ISsitematch = self.AddSiteInfo(self[1],self.rp2,metalatom)
        self.rp2.FSsitematch = self.AddSiteInfo(self[0],self.rp2,metalatom)

        self.rp1.GenNamewithsite(self[0].FP,self[1].FP,rclist0)
        self.rp2.GenNamewithsite(self[1].FP,self[0].FP,rclist0)

        return self.rp1,self.rp2 


    def SaveSurfaceMetalBmx(self):
        ISsurfacemetal =[]
        FSsurfacemetal =[]
        metalatom = []

        for iatom,atom in enumerate(self[0].atom):
            if atom.ele > 18:
                metalatom.append(iatom)
                ISsurfacemetal.append(atom)
                FSsurfacemetal.append(self[1].atom[iatom])

        ISmetalbmx2D = np.array([[0 for i in range(len(metalatom))] for i in range(len(metalatom))]) 
        FSmetalbmx2D = np.array([[0 for i in range(len(metalatom))] for i in range(len(metalatom))]) 

        ISsupercellinfo = {}
        FSsupercellinfo = {}
        for i,reali in enumerate(metalatom):
            for j,realj in enumerate(metalatom):
                ISmetalbmx2D[i][j]= self[0].bmx2D[reali][realj]
                FSmetalbmx2D[i][j]= self[1].bmx2D[reali][realj]
                if reali == realj: continue
                ISsupercellinfo[i,j] = self[0].supercellinfo[reali,realj]
                FSsupercellinfo[i,j] = self[1].supercellinfo[reali,realj]
        return len(metalatom),metalatom,ISmetalbmx2D,FSmetalbmx2D, ISsurfacemetal, FSsurfacemetal, ISsupercellinfo, FSsupercellinfo

    def AddSiteInfo(self,Str,rp,metalatom):
        newmatch =[]

        for pair in Str.sitematch:
            for i,atominfo in enumerate(rp.atommatch):
                if atominfo[1] == pair[0]:
                    for j,realatomj in enumerate(metalatom):
                        if realatomj == pair[1]:
                            newmatch.append([i,j])
                            break
                    break
        #print "Str.sitematch"
        #print rp.atommatch
        #print Str.sitematch
        #print newmatch
        return newmatch




    def AddSurfaceSiteBond(self,rp ,reactstr,trace):
        for pair in rp.sitematch:
            reali = trace[pair[0]]
            realj = reactstr.norganicatom+ pair[1]
            newstr.bmx2D[reali][realj] = 1
            newstr.bmx2D[realj][reali] = 1

        
    def GenNewSurfaceStr_bmx(self,rp,reactstr,FP,matchcombine):
        newstr = BmxStr()
        newstr.AddMetalBmx(reactstr,rp.nmetal,rp.ISmetalbmx2D,rp.ISsurfacemetal,rp.ISsupercellinfo)
        isect = 0

        #print rp.ISsupercellinfo
        #for fp in matchcombine:
        #    print fp.coreid

        trace = []
        for atominfo in rp.atommatch:
            if len(atominfo[0]) == 1 and (atominfo[0][0]== 0):
                isect = isect+1
            atomid = rp.GenMatch(matchcombine[isect-1].coreid,FP,atominfo[0])
            if atomid in trace:
                #print ('atomid: %d'%atomid)
                #print 'wrong match in GenNewStr'
                return -1
            trace.append(atomid)
        #print 'trace',trace
        #print 'trace',trace
        newstr.addmetal = []
        newstr.sitematch = []
        newstr.adsorpatom = {}

        #print "rp.ISsitematch:"
        #print rp.ISsitematch
        for pair in rp.ISsitematch:
            reali = trace[pair[0]]
            realj = reactstr.norganicatom+ pair[1]
            newstr.bmx2D[reali][realj] = 1
            newstr.bmx2D[realj][reali] = 1
            #newstr.addmetal.append(realj)
            #print reali,realj 
            newstr.sitematch.append([reali,realj])
            if realj not in newstr.addmetal: 
                newstr.addmetal.append(realj)
            if reali not in newstr.adsorpatom.keys():
                newstr.adsorpatom[reali]= [realj]
            elif realj not in newstr.adsorpatom[reali]:
                newstr.adsorpatom[reali].append(realj)

        return newstr


    def StrCheckSame(self,modified,stdout,depth =3):
        #print ('start fp1')
        fp1 =ECFP(modified,depth)
        fp1.GenECFP()
        #print fp1.allindex
        #print modified.bmx2D

        #print ('start fp2')
        fp2 =ECFP(stdout,depth)
        fp2.GenECFP()
        #print fp2.allindex
        #print stdout.bmx2D

        lmatch  = 1
        for i,index in enumerate(fp2.allindex):
            if fp1.allindex[i] != index:
                lmatch = 0
#                print fp1.allindex
#                print fp2.allindex
#                print index,fp2.IndextoFp[index][0].coreid,fp2.IndextoAtomset[index]
#                print fp1.genbylayer[2][fp2.IndextoFp[index][0].coreid].trace
#                print fp1.genbylayer[2][fp2.IndextoFp[index][0].coreid].allatomset
#                print fp1.genbylayer[1][32].trace
#                print fp2.genbylayer[1][32].trace
#                sys.exit()
                break
        return lmatch

