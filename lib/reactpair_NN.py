#from allstr_uncm import AllStr
import os
import numpy as np
import ctypes
from ctypes import pointer
from reactpair_pure import Pair as pair_nosurface
from reactpair_surfacewithsite import ReactionPair as pair_withsite
from ECFP_withsite import ECFP as ECFPwithsite
from ECFP_withsite import AllStr as allstrwithsite
from ECFP_withsite import SurfaceStr as SurfaceStr_ecfpwithsite
from ECFP0 import ECFP as ECFP0

def wrapCalISFSdistance(datain):
    pair = datain
    #print str.sminame
    dis = pair.ISFSdistance()
    return dis

def wrapQuickdoublename(pair, str1_withsitebmx,str2_withsitebmx, Lreturn , parisortflag):
    r = pair.Quickdoublename(str1_withsitebmx,str2_withsitebmx, Lreturn , parisortflag)
    return r

def reversetype(type):
    t = type.split()
    if (t[-1]== 'break'):
        reverse = 'formation'
        r = t[0] + ' '+ reverse
    elif(t[-1]== 'formation'):
        reverse = 'break'
        r = t[0] + ' '+ reverse
    else:
        r = 'unknown'

    return r 

class Pair(pair_nosurface):
    def __init__(self):
        list.__init__(self)

    def AddStr(self,str1,str2,TS= 0, TSstr=False):#, sortflag = True):
        self.append(str1)
        self.append(str2)
        #self.append(TSstr)
        self.TS= TS
        self.TSstr =TSstr

        #if sortflag: self.sort(key= lambda X: X.name)
        #self.name = '%s'%(self[0].name+self[1].name)
        return


    def GenName0(self,str1_withsitebmx,str2_withsitebmx,sortflag=True):

        if sortflag:
            zipped = zip([self[0],self[1]],[str1_withsitebmx,str2_withsitebmx])
            _sort =  zipped.sort(key= lambda X: X[0].ECFPname)
            r = zip(*zipped)
            str1_withsitebmx = r[1][0]
            str2_withsitebmx = r[1][1]
            self.sort(key= lambda X: X.ECFPname)
        self.ECFP0name = '%s'%(self[0].ECFPname+self[1].ECFPname)
        return str1_withsitebmx,str2_withsitebmx

    def ExtractPureReaction(self):
        _tmp= pair_nosurface()
        _tmp.AddStr(self[0],self[1],sortflag= False)
        #r =_tmp.GetReactCenter()
        r = _tmp.GetReactCenter(depth = 1, Ltype = 1, LaddfakeZPE = 1)
        if isinstance(r,int):
            print 'error in extract react center'
            return -1
        else:
            _tmpstr1,_tmpstr2,groupdict1,groupdict2,groupdict01,groupdict02,atomlist0,atomlist, reacttype, ZPE_TS, ZPE_heat =r

        if ZPE_TS < 0:
            ZPE_TS2 = 0
        else:
            ZPE_TS2 = -ZPE_heat
        ZPE_heat2 = -ZPE_heat


        self.ISreacttype = reacttype
        self.FSreacttype = reversetype(reacttype)
        self.ZPE_TS_IS2FS = ZPE_TS
        self.ZPE_heat_IS2FS = ZPE_heat
        self.ZPE_TS_FS2IS = ZPE_TS2
        self.ZPE_heat_FS2IS = ZPE_heat2

        self.rclist0 = atomlist0

        addlist = []
        for iatom,atom in enumerate(self[0].atom):
            if self[0].group[iatom] in groupdict1.keys():
                addlist.append(iatom)
        for iatom,atom in enumerate(self[1].atom):
            if (self[1].group[iatom] in groupdict2.keys()) and (iatom not in addlist):
                print ('ttout')
                addlist.append(iatom)
        for iatom,atom in enumerate(self[0].atom):
            if atom.ele> 18:
                addlist.append(iatom)

        outstr1 = self[0].GenSubStr(addlist)
        outstr2 = self[1].GenSubStr(addlist)
        TSstr = self.TSstr.GenSubStr(addlist)
        return outstr1,TSstr, outstr2,reacttype

    def Quickdoublename(self,str1_withsitebmx,str2_withsitebmx, lreturn = 0,sortflag =True):
        str1_withsitebmx,str2_withsitebmx = self.GenName0(str1_withsitebmx,str2_withsitebmx, sortflag)

        _tmp= pair_nosurface()
        _tmp.AddStr(self[0],self[1],sortflag= False)
        #r =_tmp.GetReactCenter()
        r =_tmp.GenReactPattern()
        if isinstance(r,int):
            print 'fail to extract react pattern '
            return -1
        else:
            self.rp1_0,self.rp2_0,self.rclist0,self.rclist1 =r

        self[0].ECFP0name = _tmp[0].ECFPname
        self[1].ECFP0name = _tmp[1].ECFPname
        self[0].ECFP0 = _tmp[0].ECFP
        self[1].ECFP0 = _tmp[1].ECFP

        self.ECFP0name = '%s'%(_tmp[0].ECFPname+_tmp[1].ECFPname)

        _tmp= pair_withsite()
        _tmp.AddStr(str1_withsitebmx,str2_withsitebmx,sortflag= False)
        _tmp.GenECFP_withsite(self.rclist0)
        self.rp1,self.rp2 =_tmp.SitePatchtoRP(self.rp1_0,self.rp2_0,self.rclist0)

        self[0].ECFPname_withsite= _tmp[0].ECFPname
        self[1].ECFPname_withsite= _tmp[1].ECFPname
        self[0].ECFP_withsite= _tmp[0].FP
        self[1].ECFP_withsite= _tmp[1].FP
        self[0].sitematch = _tmp[0].sitematch
        self[1].sitematch = _tmp[1].sitematch

        self.ECFPname_withsite = '%s'%(_tmp[0].ECFPname+_tmp[1].ECFPname)

        if lreturn:
            return self.rp1_0,self.rp2_0,self.rclist0,self.rclist1, self[0].ECFP0name,self[1].ECFP0name,self[0].ECFP0,self[1].ECFP0,self.ECFP0name,\
                    self[0].ECFPname_withsite, self[1].ECFPname_withsite, self[0].ECFP_withsite, self[1].ECFP_withsite,self[0].sitematch,self[1].sitematch,\
                    self.ECFPname_withsite, self.rp1,self.rp2
        else:
            return  1

    def doublename(self,str1_withsitebmx,str2_withsitebmx, sortflag= True):

        str1_withsitebmx,str2_withsitebmx = self.GenName0(str1_withsitebmx,str2_withsitebmx, sortflag)

        _tmp= pair_nosurface()
        _tmp.AddStr(self[0],self[1],sortflag= False)
        #r =_tmp.GetReactCenter()
        r =_tmp.GenReactPattern()
        if isinstance(r,int):
            print 'fail to extract react pattern '
            return -1
        else:
            self.rp1_0,self.rp2_0,self.rclist0,self.rclist1 =r
    
    
        l1 =_tmp.CheckReactPattern_new(self[0],self.rp1_0,self[1],depth =3)
        l2 =_tmp.CheckReactPattern_new(self[1],self.rp2_0,self[0],depth =3)
    
        if l1 != 1 or l2 != 1:
            print 'wrong in gen rp0'
            return -2

        self[0].ECFP0name = _tmp[0].ECFPname
        self[1].ECFP0name = _tmp[1].ECFPname
        self[0].ECFP0 = _tmp[0].ECFP
        self[1].ECFP0 = _tmp[1].ECFP

        self.ECFP0name = '%s'%(_tmp[0].ECFPname+_tmp[1].ECFPname)
        l = self.AddSiteInfo(str1_withsitebmx,str2_withsitebmx)

        return l

    def ChangeMode(self,mode):
        if mode =='ECFP0':
            self[0].name = self[0].ECFP0name 
            self[1].name = self[1].ECFP0name
            self.name = self.ECFP0name
            self[0].ECFP =self[0].ECFP0
            self[1].ECFP =self[1].ECFP0


        elif mode =='ECFP_withsite':
            self[0].name = self[0].ECFPname_withsite
            self[1].name = self[1].ECFPname_withsite
            self.name = self[0].ECFP0name + self[1].ECFP0name + self.rp1.ECFPname_withsite + '-' + self.surfaceinfo
            self[0].ECFP =self[0].ECFP_withsite
            self[1].ECFP =self[1].ECFP_withsite
        
        elif mode == 'ECFP_surface':
            self[0].name = self[0].ECFPname_surface
            self[1].name = self[1].ECFPname_surface
            self.name = self.ECFPname_surface
            self[0].ECFP =self[0].ECFP0
            self[1].ECFP =self[1].ECFP0


        return


    def AddECFPname_withsurface(self):
        self[0].ECFPname_surface  = self[0].ECFP0name + '-' + self.surfaceinfo
        self[1].ECFPname_surface  = self[1].ECFP0name + '-' + self.surfaceinfo
        self.ECFPname_surface  = self.ECFP0name + '-' + self.surfaceinfo
        self.rp1.surfaceinfo = self.surfaceinfo
        self.rp2.surfaceinfo = self.surfaceinfo

    def AddSiteInfo(self,str1,str2):
        _tmp= pair_withsite() 
        _tmp.AddStr(str1,str2,sortflag= False) 
        _tmp.GenECFP_withsite(self.rclist0)
        self.rp1,self.rp2 =_tmp.SitePatchtoRP(self.rp1_0,self.rp2_0,self.rclist0)

        l1 =self.CheckAddSurface(self[0],str1,self.rp1)
        l2 =self.CheckAddSurface(self[1],str2,self.rp2)

        #if l1 != 1 :
        #    print ("wrong in addsiteinfo 1 %s" %self[0].sminame)
        #if l2 != 1 :
        #    print ("wrong in addsiteinfo 2 %s" %self[1].sminame)

        if l1 != 1 or l2 != 1:
            print 'wrong in addsiteinfo'
            return -2


        self[0].ECFPname_withsite= _tmp[0].ECFPname
        self[1].ECFPname_withsite= _tmp[1].ECFPname
        self[0].ECFP_withsite= _tmp[0].FP
        self[1].ECFP_withsite= _tmp[1].FP
        self[0].sitematch = _tmp[0].sitematch
        self[1].sitematch = _tmp[1].sitematch

        self.ECFPname_withsite = '%s'%(_tmp[0].ECFPname+_tmp[1].ECFPname)

        return 1


    def CheckAddSurface(self,strin,stdout,rp,depth =3):
        plist =self.AddSurfacesite(strin,rp,depth,CalName =False)

        rightp = []
        _tmp =pair_withsite()
        for product in plist:
            lmatch =_tmp.StrCheckSame(product,stdout,depth)
            if lmatch == 1:
                rightp.append(plist)
        if len(rightp) >0:
            return 1
        else:
        #    print 'wrong surface site addition'
            return 0


    def AddSurfacesite(self,strin,rp,depth = 3,CalName =False):
        FP = ECFP0(strin,depth)
        FP.GenECFP()

        matchcombine =rp.GetMatchFp(FP)
        #print matchcombine
        #for fplist in matchcombine:
        #    for fp in fplist:
        #        print fp.coreele
        #        print fp.coreid
        #print rp.atomall

        newstrlist = allstrwithsite()
        
        if matchcombine == -1:
            print 'fail to match fp'
            return newstrlist

        for matchfp in matchcombine:
            _tmp = pair_withsite()
            r=  _tmp.GenNewSurfaceStr_bmx(rp,strin,FP,matchfp)
            if r == -1: 
                #print 'fail to GenStr'
                continue
            newstrlist.append(r)

        #print len(newstrlist)

        if (CalName):
            newstrlist.GetAllsminame_fromreactstr()

        return newstrlist


    def ISFSdistance(self, program = '/home10/kpl/pymodule/reactioncoord/reactioncoord.so'):
        self[0].calCtypes()
        self[1].calCtypes()

        self.factxyz = pointer((ctypes.c_double*9)(*[0,0,0,0,0,0,0,0,0]))
        self.length = pointer((ctypes.c_double)(0))
        len3n3 = self[0].natom * 3 + 9
        self.fracIS = pointer((ctypes.c_double*len3n3)(*range(len3n3)))
        self.fracFS = pointer((ctypes.c_double*len3n3)(*range(len3n3)))


     #   print self[0].c_rv.contents
     #   for item in self[0].c_rv.contents:
     #       print item
     #   print self[0].c_xa.contents
     #   for item in self[0].c_xa.contents:
     #       print item
        self.rcal = ctypes.cdll.LoadLibrary(program)
        self.rcal.isfslength_new_(self[0].c_natm, self[0].c_rv, self[0].c_xa, self[1].c_rv,
                         self[1].c_xa, self.factxyz, self.length, self.fracIS, self.fracFS)
        #print self.length.contents
        #print self.length[0]
        self.totdis = self.length[0]
        return self.length[0]

    def ISFSdistance_crystal(self, program = '/home10/kpl/pymodule/reactioncoord/reactioncoord.so'):
        self[0].calCtypes()
        self[1].calCtypes()

        self.factxyz = pointer((ctypes.c_double*9)(*[0,0,0,0,0,0,0,0,0]))
        self.length = pointer((ctypes.c_double)(0))
        len3n3 = self[0].natom * 3 + 9
        self.fracIS = pointer((ctypes.c_double*len3n3)(*range(len3n3)))
        self.fracFS = pointer((ctypes.c_double*len3n3)(*range(len3n3)))


        self.rcal = ctypes.cdll.LoadLibrary(program)
        self.rcal.isfslength_crystal_(self[0].c_natm, self[0].c_rv, self[0].c_xa, self[1].c_rv,
                         self[1].c_xa, self.factxyz, self.length, self.fracIS, self.fracFS)
        #print self.length.contents
        #print self.length[0]
        self.totdis = self.length[0]
        return self.length[0]


