#! /usr/bin/env python

import sys
import re
from reactpair_NN import Pair
import traceback
import pickle
import os
import numpy as np
from ECFP_withsite import ECFP as ECFPwithsite
from ECFP_withsite import AllStr as allstr_ecfpwithsite
from ECFP0 import ECFP as ECFP0
from ECFP0 import AllStr as allstr_ecfp0
from reactpair_NN import Pair,wrapQuickdoublename
from multiprocessing import Pool
from analyze_surfacemode_NNtmp import AllPair as allpair
from ECFP_allbond import AllStr as allstr_allbond
from simple_func_k import Loadpkl,Dumppkl,Average,CheckFile,Shuffle_listlike
import random
import time
import glob
from dataprint import *

def GetFakek(Ea):
    #return 1000000*np.exp(-Ea*96485/8.314/1000)
    #return 1000000*np.exp(-Ea*10)
    return np.exp(-Ea*10)
    #return 1

class FakeRP(object):
    def __init__(self,name1, name2, barrier, heat, surfaceinfo):
        self.ECFPname_withsite = name1
        self.ECFP0name = name2
        self.barrier = barrier
        self.heat = heat
        self.surfaceinfo = surfaceinfo

class RangeInfo():
    def __init__(self,rp =False):
        if rp:
            self.barriermin = rp.barrier
            self.barriermax = rp.barrier
            self.barrieravg = rp.barrier
            self.heatmin = rp.heat
            self.heatmax = rp.heat
            self.heatavg = rp.heat
            self.count = 1
        else:
            self.barriermin = 999
            self.barriermax = 0
            self.barrieravg = 0
            self.heatmin    = 999
            self.heatmax    = -999
            self.heatavg    = 0
            self.count = 0


    def Update(self, rp):
        if rp.barrier > self.barriermax:
            self.barriermax = rp.barrier
        if rp.barrier < self.barriermin:
            self.barriermin = rp.barrier
        if rp.heat > self.heatmax:
            self.heatmax = rp.heat
        if rp.heat < self.heatmin:
            self.heatmin = rp.heat
        self.barrieravg = (self.barrieravg*self.count + rp.barrier)/(self.count+1)
        self.heatavg =    (self.heatavg*self.count   +  rp.heat)/(self.count+1)
        self.count  =self.count+1

    
    def UpdatefromList(self,rplist):
        self.barriermin = min([rp.barrier for rp in rplist])
        self.barriermax = max([rp.barrier for rp in rplist])
        self.barrieravg = Average([rp.barrier for rp in rplist])
        self.heatmin = min([rp.heat for rp in rplist])
        self.heatmax = max([rp.heat for rp in rplist])
        self.heatavg = Average([rp.heat for rp in rplist])
        self.count = len(rplist)

class AllPair(allpair):
    def __init__(self):
        list.__init__(self)
        self.Lallmin = 0
        
    def PreScreen(self,filename='allgoodsect.arc',filemode = 1, pairsortflag =True):
        _tmp= allstr_ecfp0()
        _tmp.readfile(filename)
        _tmp.GetAllsminame(numproc=6,colorflag = 1)
        _tmp.GetAllECFPname(numproc=6)


        for i in range(0,len(_tmp),3):
            if filemode == 1:
                _tmp[i].screenuppersurf()
                _tmp[i+1].screenuppersurf()
                if (_tmp[i].upper ==1 ) or (_tmp[i+1].upper == 1):  continue
                if (not _tmp1[i].sminame) or (not _tmp[i+1].sminame):  continue
                if (_tmp[i].ECFPname == _tmp[i+1].ECFPname):  continue
                _tmp[i].name   = _tmp[i].ECFPname + '-' + _tmp[i].label 
                _tmp[i+1].name = _tmp[i+1].ECFPname + '-' + _tmp[i+1].label 

                TS= _tmp[i+2].energy
                if (_tmp[i].energy > TS) or (_tmp[i+1].energy > TS): continue
                p =Pair()
                p.AddStr(_tmp[i],_tmp[i+1], TS,_tmp[i+2])
                self.append(p)

            if filemode == 2:
                _tmp[i].screenuppersurf()
                _tmp[i+2].screenuppersurf()
                if (_tmp[i].upper ==1 ) or (_tmp[i+2].upper == 1):  continue
                if (not _tmp[i].sminame) or (not _tmp[i+1].sminame):  continue
                if (_tmp[i].ECFPname == _tmp[i+2].ECFPname):  continue
                _tmp[i].name   = _tmp[i].ECFPname + '-' + _tmp[i].label 
                _tmp[i+2].name = _tmp[i+2].ECFPname + '-' + _tmp[i+1].label 

                TS= _tmp[i+1].energy
                if (_tmp[i].energy > TS) or (_tmp[i+2].energy > TS): continue
                p =Pair()
                p.AddStr(_tmp[i],_tmp[i+2], TS,_tmp[i+1])
                self.append(p)


        self.AllName()
        self.Strsavebyname = {}
        for name in self.allname.keys():
            if len(self.allstrbyname[name]) < 2:
                self.Strsavebyname[name] = 0
            else:
                self.Strsavebyname[name] = 1
        wronglist = allstr_ecfp0()
        printlist = allstr_ecfp0()

        for pair in self:
            if (self.Strsavebyname[pair[0].name] ==0) and (self.Strsavebyname[pair[1].name] ==0):
                wronglist.append(pair[0])
                wronglist.append(pair.TSstr)
                wronglist.append(pair[1])
            else:
                printlist.append(pair[0])
                printlist.append(pair.TSstr)
                printlist.append(pair[1])
           
        printlist.printall('allgoodsect.arc_out')
        wronglist.printall('rarereact.arc')

    def Readfile_Simple(self,filename='allgoodsect.arc',filemode = 1,Lallmininput = 0,allminfile= 'allname.arc', pairsortflag = True):
        _tmp1= allstr_ecfp0()
        _tmp1.readfile(filename)
        _tmp1.GetTmpFakebmx(numproc=6,colorflag = 1)

        outstr = allstr_ecfp0()
        for i in range(0,len(_tmp1),3):
            if filemode ==2:
                TS = _tmp1[i+1].energy
                p= Pair()
                p.AddStr(_tmp1[i],_tmp1[i+2],TS,_tmp1[i+1])
                self.append(p)
                out1,TSstr,out2, rtype = p.ExtractPureReaction( )
                outstr.append(out1)
                outstr.append(TSstr)
                outstr.append(out2)
        #outstr.printall('purereact.arc')
        return

           
    def readfile(self,filename='allgoodsect.arc',filemode = 1,Lallmininput = 0,allminfile= 'allname.arc', pairsortflag = True):
        self.Lallmin = Lallmininput

        t0 = time.time()

        _tmp2= allstr_ecfpwithsite()
        _tmp2.readfile(filename)
        _tmp2.GetTmpFakebmx(numproc=6,colorflag = 1)

        t1 = time.time()
        print ('smi time1: %f' %(t1-t0))


        _tmp1= allstr_ecfp0()
        _tmp1.readfile(filename)
        _tmp1.GetAllsminame(numproc=6,colorflag = 1)

        t2 = time.time()
        print ('smi time2: %f' %(t2-t1))
        _tmp1.GetAllECFPname(numproc=6)

        t3 = time.time()
        print ('ECFP time: %f' %(t3-t2))

        wronglist =allstr_ecfpwithsite()
        for i in range(0,len(_tmp1),3):

            if filemode == 1:
                #_tmp1[i].screenuppersurf()
                #_tmp1[i+1].screenuppersurf()
                #if (_tmp1[i].upper ==1 ) or (_tmp1[i+1].upper == 1):  continue
                #if (not _tmp[i].Lminstr) or (not _tmp[i+1].Lminstr): continue
                #if (not _tmp1[i].sminame) or (not _tmp1[i+1].sminame):  continue
                TS = _tmp1[i+2].energy
                #if (_tmp1[i].energy > TS) or (_tmp1[i+1].energy > TS): continue
                p= Pair()
                p.AddStr(_tmp1[i],_tmp1[i+1],TS,_tmp1[i+2])
                l = p.Quickdoublename(_tmp2[i],_tmp2[i+1],0 ,pairsortflag)

            elif filemode ==2:

                #_tmp1[i].screenuppersurf()
                #_tmp1[i+2].screenuppersurf()
                #if (_tmp1[i].upper ==1 ) or (_tmp1[i+2].upper == 1):  continue
                #if (not _tmp1[i].sminame) or (not _tmp1[i+2].sminame):  continue

                TS = _tmp1[i+1].energy
                #if (_tmp1[i].energy > TS) or (_tmp1[i+2].energy > TS): continue
                p= Pair()
                p.AddStr(_tmp1[i],_tmp1[i+2],TS,_tmp1[i+1])
                l = p.Quickdoublename(_tmp2[i],_tmp2[i+2],0 ,pairsortflag)
                #l = p.doublename(_tmp2[i],_tmp2[i+2])

            if l == 1:
                p.Lextract = True
                self.append(p)
            else:
                print('wrong in Gen FP, pair %d'%(i/3))
                wronglist.append(_tmp1[i])
                wronglist.append(_tmp1[i+1])
                wronglist.append(_tmp1[i+2])

        print ('all pair extract: %d'%(len(self)))
        wronglist.printall('wrongpair.arc')
        self.Arcprint('extract.arc')

        t4 = time.time()
        print ('extract rp time: %f' %(t4-t3))




    def ReadFile_parallel(self,filename='allgoodsect.arc',filemode = 1,Lallmininput = 0,allminfile= 'allname.arc', pairsortflag = True):
        self.Lallmin = Lallmininput

        t0 = time.time()

        _tmp2= allstr_ecfpwithsite()
        _tmp2.readfile(filename)
        _tmp2.GetTmpFakebmx(numproc=6,colorflag = 1)

        t1 = time.time()
        print ('smi time1: %f' %(t1-t0))


        _tmp1= allstr_ecfp0()
        _tmp1.readfile(filename)
        _tmp1.GetAllsminame(numproc=6,colorflag = 1)

        t2 = time.time()
        print ('smi time2: %f' %(t2-t1))
        _tmp1.GetAllECFPname(numproc=6)

        t3 = time.time()
        print ('ECFP time: %f' %(t3-t2))

        wronglist =allstr_ecfpwithsite()
        pool = Pool(processes=6)

        _t = []
        result= []
        for i in range(0,len(_tmp1),3):

            if filemode == 1:
                #_tmp1[i].screenuppersurf()
                #_tmp1[i+1].screenuppersurf()
                #if (_tmp1[i].upper ==1 ) or (_tmp1[i+1].upper == 1):  continue
                #if (not _tmp[i].Lminstr) or (not _tmp[i+1].Lminstr): continue
                #if (not _tmp1[i].sminame) or (not _tmp1[i+1].sminame):  continue

                TS = _tmp1[i+2].energy
                if (_tmp1[i].energy > TS) or (_tmp1[i+1].energy > TS): continue
                p= Pair()
                p.AddStr(_tmp1[i],_tmp1[i+1],TS,_tmp1[i+2])
                #l = p.Quickdoublename(_tmp2[i],_tmp2[i+1])
                result.append(pool.apply_async(wrapQuickdoublename,args= (p, _tmp2[i],_tmp2[i+1],1,pairsortflag)))
                _t.append(p)

            elif filemode ==2:

                #_tmp1[i].screenuppersurf()
                #_tmp1[i+2].screenuppersurf()
                #if (_tmp1[i].upper ==1 ) or (_tmp1[i+2].upper == 1):  continue
                #if (not _tmp1[i].sminame) or (not _tmp1[i+2].sminame):  continue

                TS = _tmp1[i+1].energy
                if (_tmp1[i].energy > TS) or (_tmp1[i+2].energy > TS): continue
                p= Pair()
                p.AddStr(_tmp1[i],_tmp1[i+2],TS,_tmp1[i+1])
                #l = p.Quickdoublename(_tmp2[i],_tmp2[i+2])
                #l = p.doublename(_tmp2[i],_tmp2[i+2])
                result.append(pool.apply_async(wrapQuickdoublename,args= (p, _tmp2[i],_tmp2[i+2],1,pairsortflag)))
                #if pairsortflag: p.sort(key= lambda X: X.ECFPname)
                _t.append(p)
        pool.close()
        pool.join()

        for pair,re in zip(_t,result):
            r = re.get()
            if not isinstance(r,int):
                try:
                    pair.Lextract = True
                    if pairsortflag: pair.sort(key= lambda X: X.ECFPname)
                    pair.rp1_0,pair.rp2_0,pair.rclist0,pair.rclist1, pair[0].ECFP0name,pair[1].ECFP0name,pair[0].ECFP0,pair[1].ECFP0,pair.ECFP0name,\
                    pair[0].ECFPname_withsite, pair[1].ECFPname_withsite, pair[0].ECFP_withsite, pair[1].ECFP_withsite,pair[0].sitematch,pair[1].sitematch,\
                    pair.ECFPname_withsite, pair.rp1,pair.rp2 = r
                    self.append(pair)
                except:
                    print('wrong in Get result')
                    wronglist.append(pair[0])
                    wronglist.append(pair.TSstr)
                    wronglist.append(pair[1])

            else:
                print('wrong in Gen FP')
                wronglist.append(pair[0])
                wronglist.append(pair.TSstr)
                wronglist.append(pair[1])

        print ('all pair extract: %d'%(len(self)))
        wronglist.printall('wrongpair.arc')
        self.Arcprint('extract.arc')

        t4 = time.time()
        print ('extract rp time: %f' %(t4-t3))


    def ChangeMode(self,mode):
        for pair in self:
            pair.ChangeMode(mode)

    def GenECFPwithsitedict(self):
        self.ECFPwithsitedict = {}
        for pair in self:
            if pair[0].ECFPname_withsite not in self.ECFPwithsitedict.keys():
                self.ECFPwithsitedict[pair[0].ECFPname_withsite] = pair[0].ECFP_withsite

            if pair[1].ECFPname_withsite not in self.ECFPwithsitedict.keys():
                self.ECFPwithsitedict[pair[1].ECFPname_withsite] = pair[1].ECFP_withsite
 
    def GenECFP0dict(self):
        self.ECFP0dict = {}
        for pair in self:
            if pair[0].ECFP0name not in self.ECFP0dict.keys():
                self.ECFP0dict[pair[0].ECFP0name] = pair[0].ECFP0

            if pair[1].ECFP0name not in self.ECFP0dict.keys():
                self.ECFP0dict[pair[1].ECFP0name] = pair[1].ECFP0

 
    def AllRPDict(self):
        self.allrpdict_withsite ={}
        self.allrpdict_withsite_pairexample = {}
        self.allrpdict_withsite_rangeinfo = {}

        for pair in self:
            if pair.rp1.ECFPname_withsite not in self.allrpdict_withsite.keys():
                self.allrpdict_withsite[pair.rp1.ECFPname_withsite] = pair.rp1
                self.allrpdict_withsite_pairexample[pair.rp1.ECFPname_withsite] = [pair]
                self.allrpdict_withsite_rangeinfo[pair.rp1.ECFPname_withsite] =RangeInfo(pair.rp1)

            else:
                if pair.rp1.barrier < self.allrpdict_withsite[pair.rp1.ECFPname_withsite].barrier:
                    self.allrpdict_withsite[pair.rp1.ECFPname_withsite] = pair.rp1
                self.allrpdict_withsite_pairexample[pair.rp1.ECFPname_withsite].append(pair)
                self.allrpdict_withsite_rangeinfo[pair.rp1.ECFPname_withsite].Update(pair.rp1)

            if pair.rp2.ECFPname_withsite not in self.allrpdict_withsite.keys():
                self.allrpdict_withsite[pair.rp2.ECFPname_withsite] = pair.rp2
                self.allrpdict_withsite_pairexample[pair.rp2.ECFPname_withsite] = [pair]
                self.allrpdict_withsite_rangeinfo[pair.rp2.ECFPname_withsite] =RangeInfo(pair.rp2)

            else:
                if pair.rp2.barrier < self.allrpdict_withsite[pair.rp2.ECFPname_withsite].barrier:
                    self.allrpdict_withsite[pair.rp2.ECFPname_withsite] = pair.rp2
                self.allrpdict_withsite_pairexample[pair.rp2.ECFPname_withsite].append(pair)
                self.allrpdict_withsite_rangeinfo[pair.rp2.ECFPname_withsite].Update(pair.rp2)

        for key in self.allrpdict_withsite.keys():
            self.allrpdict_withsite_pairexample[key].sort(key = lambda x: x.TS)
            self.allrpdict_withsite[key].rangeinfo = self.allrpdict_withsite_rangeinfo[key]
            self.allrpdict_withsite[key].nsample = len(self.allrpdict_withsite_pairexample[key])

        return

    def AllRPDict0(self):
        self.allrpdict0 ={}
        self.allrpdict0_all ={}
        self.allrpdict0_count = {}
        for rp in self.allrpdict_withsite.values():
            if rp.ECFP0name not in self.allrpdict0_all.keys():
                self.allrpdict0_all[rp.ECFP0name] = [rp]
                self.allrpdict0_count[rp.ECFP0name] = rp.nsample
            else:
                self.allrpdict0_all[rp.ECFP0name].append(rp)
                self.allrpdict0_count[rp.ECFP0name] = self.allrpdict0_count[rp.ECFP0name] + rp.nsample


        #self.allrpdict0_rangeinfo = {}
        self.allrpdict0_ex = {}
        for key in self.allrpdict0_all.keys():
            self.allrpdict0_all[key].sort(key = lambda x: x.barrier)
            self.allrpdict0[key] = [rp.ECFPname_withsite for rp in self.allrpdict0_all[key]]
            self.allrpdict0_ex[key] = self.allrpdict0_all[key][0]
            #_tmprangeinfo = RangeInfo()
            #_tmprangeinfo.UpdatefromList(self.allrpdict0_all[key])
            #self.allrpdict0_rangeinfo[key] = _tmprangeinfo

        allrpdict0_rangeinfo = {}
        for pair in self:
            if pair.rp1.ECFP0name not in allrpdict0_rangeinfo.keys():
                allrpdict0_rangeinfo[pair.rp1.ECFP0name] = RangeInfo(pair.rp1)
            else:
                allrpdict0_rangeinfo[pair.rp1.ECFP0name].Update(pair.rp1)

            if pair.rp2.ECFP0name not in allrpdict0_rangeinfo.keys():
                allrpdict0_rangeinfo[pair.rp2.ECFP0name] = RangeInfo(pair.rp2)
            else:
                allrpdict0_rangeinfo[pair.rp2.ECFP0name].Update(pair.rp2)
        self.allrpdict0_rangeinfo = allrpdict0_rangeinfo


    def RPdict0_Sampled(self):
        allrpdict_withsite_rangeinfo = {}
        print len(self)

        for pair in self:
            if (pair[0].lsampled == 0) or (pair[1].lsampled == 0) :
                print pair[0].lsampled, pair[1].lsampled
                pair.Ldoublesamped = 0
                continue
            pair.Ldoublesamped = 1

            if pair.rp1.ECFPname_withsite not in allrpdict_withsite_rangeinfo.keys():
                allrpdict_withsite_rangeinfo[pair.rp1.ECFPname_withsite] = RangeInfo(pair.rp1)
            else:
                allrpdict_withsite_rangeinfo[pair.rp1.ECFPname_withsite].Update(pair.rp1)

            if pair.rp2.ECFPname_withsite not in allrpdict_withsite_rangeinfo.keys():
                allrpdict_withsite_rangeinfo[pair.rp2.ECFPname_withsite] = RangeInfo(pair.rp2)
            else:
                allrpdict_withsite_rangeinfo[pair.rp2.ECFPname_withsite].Update(pair.rp2)

        self.allrpdict_withsite_rangeinfo_sampled = allrpdict_withsite_rangeinfo

        allrpdict0_rangeinfo = {}
        allrpdict0_ZPEinfo = {}
        allrpdict0_doublesampled = {}

        for pair in self:
            if pair.rp1.ECFP0name not in allrpdict0_rangeinfo.keys():
                allrpdict0_ZPEinfo[pair.rp1.ECFP0name] = [pair.ZPE_TS_IS2FS, pair.ZPE_heat_IS2FS]
                allrpdict0_ZPEinfo[pair.rp2.ECFP0name] = [pair.ZPE_TS_FS2IS, pair.ZPE_heat_FS2IS]
            #else:
            #    print ('addZPE', allrpdict0_ZPEinfo[pair.rp1.ECFP0name], [pair.ZPE_TS_IS2FS, pair.ZPE_heat_IS2FS] )

            if pair.Ldoublesamped == 0 :  continue
            if pair.rp1.ECFP0name not in allrpdict0_rangeinfo.keys():
                allrpdict0_rangeinfo[pair.rp1.ECFP0name] = RangeInfo(pair.rp1)
                allrpdict0_doublesampled[pair.rp1.ECFP0name] = [pair.rp1]
            else:
                allrpdict0_rangeinfo[pair.rp1.ECFP0name].Update(pair.rp1)
                allrpdict0_doublesampled[pair.rp1.ECFP0name].append(pair.rp1)

            if pair.rp2.ECFP0name not in allrpdict0_rangeinfo.keys():
                allrpdict0_rangeinfo[pair.rp2.ECFP0name] = RangeInfo(pair.rp2)
                allrpdict0_doublesampled[pair.rp2.ECFP0name] = [pair.rp2]

            else:
                allrpdict0_rangeinfo[pair.rp2.ECFP0name].Update(pair.rp2)
                allrpdict0_doublesampled[pair.rp2.ECFP0name].append(pair.rp2)


        for key in allrpdict0_doublesampled.keys():
            allrpdict0_doublesampled[key].sort(key = lambda x: x.barrier)
        self.allrpdict0_sampled = allrpdict0_doublesampled
        self.allrpdict0_rangeinfo_sampled = allrpdict0_rangeinfo
        self.allrpdict0_ZPEinfo = allrpdict0_ZPEinfo

        print len(allrpdict_withsite_rangeinfo)
        print len(allrpdict0_rangeinfo)
        print len(allrpdict0_ZPEinfo)

        Dumppkl('allrpdict_withsite_rangeinfo_doublesampled.pkl', self.allrpdict_withsite_rangeinfo_sampled)
        Dumppkl('allrpdict0_rangeinfo_doublesampled.pkl',self.allrpdict0_rangeinfo_sampled)
        Dumppkl('allrpdict0_doublesampled.pkl',self.allrpdict0_sampled)
        Dumppkl('allrpdict0_ZPEinfo.pkl',self.allrpdict0_ZPEinfo)

        return

    def DumpECFP(self,ECFP0file = False,ECFPwithsitefile = False):
        if ECFP0file :
            if not hasattr(self,'ECFP0dict'):
                self.GenECFP0dict()
            #outputfile = open(ECFP0file,'wb')
            #pickle.dump(self.ECFP0dict,outputfile)
            #outputfile.close()


        if ECFPwithsitefile:
            if not hasattr(self,'ECFPwithsitedict'):
                self.GenECFPwithsitedict()
            #outputfile = open(ECFPwithsitefile,'wb')
            #pickle.dump(self.ECFPwithsitedict,outputfile)
            #outputfile.close()



    def DumpAllRPDict(self,RPDict0=False,RPDict0_RangeInfo = False,RPDict0_count= False,RPDict=False,RPDictwithexample=False, RP0Dict_ex =False):
        if RPDict0:
            print ('all rp0 : %d'%len(self.allrpdict0))
            outputfile = open(RPDict0,'wb')
            pickle.dump(self.allrpdict0,outputfile)
            outputfile.close()

        if RPDict0_RangeInfo:
            Dumppkl(RPDict0_RangeInfo,self.allrpdict0_rangeinfo)

        if RPDict0_count:
            Dumppkl(RPDict0_count,self.allrpdict0_count)

        if RPDict:
            print ('all rpwithsite : %d'%len(self.allrpdict_withsite))
            outputfile = open(RPDict,'wb')
            pickle.dump(self.allrpdict_withsite,outputfile)
            outputfile.close()

        if RPDictwithexample:
            outputfile = open(RPDictwithexample,'wb')
            pickle.dump(self.allrpdict_withsite_pairexample,outputfile)
            outputfile.close()
        if RP0Dict_ex:
            Dumppkl(RP0Dict_ex, self.allrpdict0_ex)



    def OutRPData(self,mode,outfile):
        self.ChangeMode(mode)
        self.AllName()
        self.GenPairDict()
        #self.GetAllECFPname_surface()

        _tmp= AllPair()
        _tmp.sampleddict = self.sampleddict 
        
        for pair in self.minpair.values():
            print pair.name
            print pair.TS
            _tmp.append(pair)
       
        for i in range(100): random.shuffle(_tmp)
        _tmp.OutDataSet(outfile)
        _tmp.Arcprint('allpair.arc')
        _tmp.AllRPDict()
        _tmp.AllRPDict0()
        #_tmp.DumpAllRPDict('allrpdict0.pkl','allrpdict0_rangeinfo.pkl','allrpdict0_count.pkl','allrpdict_withsite.pkl','allrpdict_withsite_pairexample.pkl')
        #_tmp.DumpAllRPDict('allrpdict0.pkl','allrpdict0_rangeinfo.pkl','allrpdict0_count.pkl','allrpdict_withsite.pkl')#,'allrpdict_withsite_pairexample.pkl')
        _tmp.DumpAllRPDict('allrpdict0.pkl','allrpdict0_rangeinfo.pkl','allrpdict0_count.pkl',RP0Dict_ex='allrpdict0_ex.pkl')#,'allrpdict_withsite_pairexample.pkl')
        _tmp.DumpECFP('ECFP0.pkl','ECFPwithsite.pkl')

        _tmp.OutAllpairexample()
        return  _tmp

    def OutAllpairexample(self):
        if not glob.glob('Arcsave_all'):
            os.mkdir('Arcsave_all')

        for namelist in self.allrpdict0.values():
            for rpname in namelist:
                pairlist = self.allrpdict_withsite_pairexample[rpname]
                _tmp = allstr_ecfp0()
                for pair in pairlist:
                    pair[0].label = pair[0].ECFPname_surface
                    #pair.TSstr.label = pair.rp1.ECFPname_withsite
                    pair[1].label = pair[1].ECFPname_surface
                    if pair.rp1.ECFPname_withsite == rpname:
                        pair.TSstr.label = pair.rp1.ECFPname_withsite
                        _tmp.append(pair[0])
                        _tmp.append(pair.TSstr)
                        _tmp.append(pair[1])
                    elif pair.rp2.ECFPname_withsite == rpname:
                        pair.TSstr.label = pair.rp2.ECFPname_withsite
                        _tmp.append(pair[1])
                        _tmp.append(pair.TSstr)
                        _tmp.append(pair[0])
            
                _tmp.printall('Arcsave_all/%s.arc'%rpname)

        if not glob.glob('Arcsave_all_rp0'):
            os.mkdir('Arcsave_all_rp0')

        for key,namelist in self.allrpdict0.items():
            _tmp = allstr_ecfp0()
            for rpname in namelist:
                pairlist = self.allrpdict_withsite_pairexample[rpname]
                for pair in pairlist:
                    pair[0].label = pair[0].ECFPname_surface
                    #pair.TSstr.label = pair.rp1.ECFPname_withsite
                    pair[1].label = pair[1].ECFPname_surface
                    if pair.rp1.ECFPname_withsite == rpname:
                        pair.TSstr.label = pair.rp1.ECFPname_withsite
                        _tmp.append(pair[0])
                        _tmp.append(pair.TSstr)
                        _tmp.append(pair[1])
                    elif pair.rp2.ECFPname_withsite == rpname:
                        pair.TSstr.label = pair.rp2.ECFPname_withsite
                        _tmp.append(pair[1])
                        _tmp.append(pair.TSstr)
                        _tmp.append(pair[0])

            _tmp.printall('Arcsave_all_rp0/%s.arc'%key)


        return

    def Arcprint(self,outfile,list =False):
        if not list:
            list = range(len(self))

        _tmp=allstr_ecfp0()
        for ipair in list:
            pair = self[ipair]
            #try:
            #    pair[0].label = pair[0].ECFPname_surface 
            #    pair[1].label = pair[1].ECFPname_surface
            #    pair.TSstr.label = 'RP_ECFP0name-'+'%s'%(pair.rp1.ECFP0name)+'-barrier-'+'%f'%(pair.barrier1)+'-heat-'+'%f'%(pair.heat1)

            #except:
            #    pair[0].label = ''
            #    pair[1].label = ''
            #    pair.TSstr.label = ''
            #    #pair[0].label = 'IS-%d'%ipair
            #    #pair[1].label = 'FS-%d'%ipair
            #    #pair.TSstr.label = 'TS-%d'%ipair
            _tmp.append(pair[0])
            _tmp.append(pair.TSstr)
            _tmp.append(pair[1])

        _tmp.printlist(range(len(_tmp)),outfile)

    def OutDataSet(self,outfile,list =False):
        if not list:
            list = range(len(self))

        f= open(outfile,'w')
        f.write('Raction Pair Dataset\n')
        f.write('tmptest for NN------format: kpl200101\n')
        f.write('Note: Only the chemisorption sites of reaction center are recorded and encoded in ECFPname_withste\n')
        f.write('\n')
        for ipair in list:
            pair = self[ipair]
            pair[0].lsampled = self.sampleddict[pair[0].ECFPname_surface]
            pair[1].lsampled = self.sampleddict[pair[1].ECFPname_surface]
            f.write('ReactPair:       %d\n'%(ipair+1))
            f.write('surfaceinfo:     %s\n'%pair.surfaceinfo)
            f.write('LDoubleEndSampled:     %s     %s\n'%(self.sampleddict[pair[0].ECFPname_surface],self.sampleddict[pair[1].ECFPname_surface]))
            f.write('sminame of IS:   %s\n'%pair[0].sminame)
            f.write('sminame of FS:   %s\n'%pair[1].sminame)
            f.write('ECFPname_nosurface:       %s\n'%pair.ECFP0name)
            f.write('ECFPname_nosurface of IS: %s\n'%pair[0].ECFP0name)
            f.write('ECFPname_nosurface of FS: %s\n'%pair[1].ECFP0name)
            f.write('ECFPname_surface of IS:   %s\n'%pair[0].ECFPname_surface)
            f.write('LsampledIS:     %s\n'%self.sampleddict[pair[0].ECFPname_surface])
            f.write('ECFPname_surface of FS:   %s\n'%pair[1].ECFPname_surface)
            f.write('LsampledFS:     %s\n'%self.sampleddict[pair[1].ECFPname_surface])
            f.write('ECFPname_withsite:        %s\n'%pair.ECFPname_withsite)
            f.write('ECFPname_withsite of IS:  %s\n'%pair[0].ECFPname_withsite)
            f.write('ECFPname_withsite of FS:  %s\n'%pair[1].ECFPname_withsite)
            f.write('Barrier (IStoFS)  %14.6f\n'%pair.barrier1 )
            f.write('Heat    (IStoFS)  %14.6f\n'%pair.heat1 )
            f.write('Barrier (FStoIS)  %14.6f\n'%pair.barrier2 )
            f.write('Heat    (FStoIS)  %14.6f\n'%pair.heat2 )
            f.write('Reaction Pattern(IStoFS) ECFPname_nosurface: %s\n'%pair.rp1.ECFP0name)
            f.write('Reaction Pattern(IStoFS) ECFPname_withsite:  %s\n'%pair.rp1.ECFPname_withsite)
            f.write('Reaction Pattern(IStoFS) fpmatch link:       %s\n'%pair.rp1.fpmatch)
            f.write('Reaction Pattern(IStoFS) fpmatch atommap:    %s\n'%pair.rp1.atommatch)
            f.write('Reaction Pattern(FStoIS) ECFPname_nosurface: %s\n'%pair.rp2.ECFP0name)
            f.write('Reaction Pattern(FStoIS) ECFPname_withsite:  %s\n'%pair.rp2.ECFPname_withsite)
            f.write('Reaction Pattern(FStoIS) fpmatch link:       %s\n'%pair.rp2.fpmatch)
            f.write('Reaction Pattern(FStoIS) fpmatch atommap:    %s\n'%pair.rp2.atommatch)
            if pair.Lextract:
                f.write('Reaction Center list 0 :  %s\n'%pair.rclist0)
                f.write('Reaction Center list 1 :  %s\n'%pair.rclist1)
                f.write('Sitematchpair of IS: %s\n'%pair[0].sitematch)
                f.write('Sitematchpair of FS: %s\n'%pair[1].sitematch)

            str=pair[0]
            #f.write('ISstr\n')
            f.write('     Energy     %8d    IS     %12.6f\n'%((3*ipair),str.energy))
            f.write('!DATE \n')
            f.write('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n'%
                    (str.abc[0],str.abc[1],str.abc[2],str.abc[3],str.abc[4],str.abc[5]))
            #f.write('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n'%(100,100,100,90,90,90))
            for i,atom in enumerate(str.atom):
                f.write('%-3s %15.9f %15.9f %15.9f CORE %4d %2s %2s %8.4f %4d\n' %
                        (atom.elesymbol,atom.xyz[0],atom.xyz[1],atom.xyz[2],i+1,atom.elesymbol,atom.elesymbol,atom.charge,i+1))
            f.write('end\nend\n')


            str=pair[1]
            #f.write('FSstr\n')
            f.write('     Energy     %8d    FS     %12.6f\n'%((3*ipair)+1,str.energy))
            f.write('!DATE \n')
            f.write('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n'%
                    (str.abc[0],str.abc[1],str.abc[2],str.abc[3],str.abc[4],str.abc[5]))
            #f.write('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n'%(100,100,100,90,90,90))
            for i,atom in enumerate(str.atom):
                f.write('%-3s %15.9f %15.9f %15.9f CORE %4d %2s %2s %8.4f %4d\n' %
                        (atom.elesymbol,atom.xyz[0],atom.xyz[1],atom.xyz[2],i+1,atom.elesymbol,atom.elesymbol,atom.charge,i+1))
            f.write('end\nend\n')

            str=pair.TSstr
            #f.write('TSstr\n')
            f.write('     Energy     %8d    TS     %12.6f\n'%((3*ipair)+2,str.energy))
            f.write('!DATE \n')
            f.write('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n'%
                    (str.abc[0],str.abc[1],str.abc[2],str.abc[3],str.abc[4],str.abc[5]))
            #f.write('PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n'%(100,100,100,90,90,90))
            for i,atom in enumerate(str.atom):
                f.write('%-3s %15.9f %15.9f %15.9f CORE %4d %2s %2s %8.4f %4d\n' %
                        (atom.elesymbol,atom.xyz[0],atom.xyz[1],atom.xyz[2],i+1,atom.elesymbol,atom.elesymbol,atom.charge,i+1))
            f.write('end\nend\n')

            f.write('\n') 

        f.close()

    def ReadLabel(self, inputfile):
        _arcread =allstr_allbond()
        _arcread.readfile(inputfile)

        print ('All Str read : %d'%len(_arcread))
        f =  open(inputfile,'r')
        currentPair = len(self)-1
        for line in f:
            if ('ReactPair' in line):
                self.append(Pair())
                currentPair = currentPair+ 1
                TS = _arcread[3*currentPair+2].energy
                self[currentPair].AddStr(_arcread[3*currentPair],_arcread[3*currentPair+1],TS, _arcread[3*currentPair+2])
            elif ('surfaceinfo' in line):
                self[currentPair].surfaceinfo = line.split()[-1]

            elif ('LsampledIS' in line):
                self[currentPair][0].lsampled = int(line.split()[-1])
            elif ('LsampledFS' in line):
                self[currentPair][1].lsampled = int(line.split()[-1])

            elif ('Barrier' in line) and ('IStoFS' in line):
                self[currentPair].barrier1 = float(line.split()[-1])
            elif ('Heat' in line) and ('IStoFS' in line):
                self[currentPair].heat1 = float(line.split()[-1])
            elif ('Barrier' in line) and ('FStoIS' in line):
                self[currentPair].barrier2 = float(line.split()[-1])
            elif ('Heat' in line) and ('FStoIS' in line):
                self[currentPair].heat2 = float(line.split()[-1])

            elif ('ECFPname_nosurface' in line) and ('IS' not in line) and ('FS' not in line):
                self[currentPair].ECFP0name = line.split()[-1]
            elif ('ECFPname_nosurface' in line) and ('IS' in line) and ('FS' not in line):
                self[currentPair][0].ECFP0name = line.split()[-1]
            elif ('ECFPname_nosurface' in line) and ('FS' in line) and ('IS' not in line):
                self[currentPair][1].ECFP0name = line.split()[-1]

            elif ('ECFPname_withsite' in line) and ('IS' not in line) and ('FS' not in line):
                self[currentPair].ECFPname_withsite = line.split()[-1]
            elif ('ECFPname_withsite' in line) and ('IS' in line) and ('FS' not in line):
                self[currentPair][0].ECFPname_withsite = line.split()[-1]
            elif ('ECFPname_withsite' in line) and ('FS' in line) and ('IS' not in line):
                self[currentPair][1].ECFPname_withsite = line.split()[-1]

            elif ('ECFPname_withsite' in line) and ('IStoFS' in line):
                self[currentPair].rp1name = line.split()[-1]
            elif ('ECFPname_withsite' in line) and ('FStoIS' in line):
                self[currentPair].rp2name = line.split()[-1]

            elif ('ECFPname_nosurface' in line) and ('IStoFS' in line):
                self[currentPair].rp1name0 = line.split()[-1]
            elif ('ECFPname_nosurface' in line) and ('FStoIS' in line):
                self[currentPair].rp2name0 = line.split()[-1]


        for pair in self:
            pair.rp1 = FakeRP(pair.rp1name, pair.rp1name0, pair.barrier1, pair.heat1 , pair.surfaceinfo )
            pair.rp2 = FakeRP(pair.rp2name, pair.rp2name0, pair.barrier2, pair.heat2 , pair.surfaceinfo )


    def OutRPInfo(self,outputfile ):
        f =  open(outputfile,'w')
        for ipair,pair in enumerate(self):
            f.write('Array   %d\n' %(2*ipair))
            f.write('RP0:   %s\n'%(pair.rp1name0))
            f.write('RPwithsite:  %s\n' %(pair.rp1name))
            f.write('Array   %d\n' %(2*ipair+1))
            f.write('RP0:   %s\n'%(pair.rp2name0))
            f.write('RPwithsite:  %s\n' %(pair.rp2name))
        f.close()


    def ReadRPData(self,inputfile):
        _arcread =allstr_allbond()
        _arcread.readfile(inputfile)

        print ('All Str read : %d'%len(_arcread))
        f =  open(inputfile,'r')
        currentPair = len(self)-1
        for line in f:
            if ('ReactPair' in line):
                self.append(Pair())
                currentPair = currentPair+ 1
                TS = _arcread[3*currentPair+2].energy
                self[currentPair].AddStr(_arcread[3*currentPair],_arcread[3*currentPair+1],TS, _arcread[3*currentPair+2])
            elif ('surfaceinfo' in line):
                self[currentPair].surfaceinfo = line.split()[-1]

            elif ('Barrier' in line) and ('IStoFS' in line):
                self[currentPair].barrier1 = float(line.split()[-1])
            elif ('Heat' in line) and ('IStoFS' in line):
                self[currentPair].heat1 = float(line.split()[-1])
            elif ('Barrier' in line) and ('FStoIS' in line):
                self[currentPair].barrier2 = float(line.split()[-1])
            elif ('Heat' in line) and ('FStoIS' in line):
                self[currentPair].heat2 = float(line.split()[-1])


            elif ('sminame' in line) and ('IS' in line):
                self[currentPair][0].sminame = line.split()[-1]
            elif ('sminame' in line) and ('FS' in line):
                self[currentPair][1].sminame = line.split()[-1]

            elif ('ECFPname_nosurface' in line) and ('IS' not in line) and ('FS' not in line):
                self[currentPair].ECFP0name = line.split()[-1]
            elif ('ECFPname_nosurface' in line) and ('IS' in line) and ('FS' not in line):
                self[currentPair][0].ECFP0name = line.split()[-1]
                self[currentPair][0].ECFP0 = self.ECFP0dict[self[currentPair][0].ECFP0name]
            elif ('ECFPname_nosurface' in line) and ('FS' in line) and ('IS' not in line):
                self[currentPair][1].ECFP0name = line.split()[-1]
                self[currentPair][1].ECFP0 = self.ECFP0dict[self[currentPair][1].ECFP0name]

            elif ('ECFPname_withsite' in line) and ('IS' not in line) and ('FS' not in line):
                self[currentPair].ECFPname_withsite = line.split()[-1]
            elif ('ECFPname_withsite' in line) and ('IS' in line) and ('FS' not in line):
                self[currentPair][0].ECFPname_withsite = line.split()[-1]
                self[currentPair][0].ECFP_withsite = self.ECFPwithsitedict[self[currentPair][0].ECFPname_withsite]
            elif ('ECFPname_withsite' in line) and ('FS' in line) and ('IS' not in line):
                self[currentPair][1].ECFPname_withsite = line.split()[-1]
                self[currentPair][1].ECFP_withsite = self.ECFPwithsitedict[self[currentPair][1].ECFPname_withsite]

            elif ('ECFPname_withsite' in line) and ('IStoFS' in line):
                self[currentPair].rp1 = self.allrpdict_withsite[line.split()[-1]]
            elif ('ECFPname_withsite' in line) and ('FStoIS' in line):
                self[currentPair].rp2 = self.allrpdict_withsite[line.split()[-1]]

        f.close()

        for pair in self:
            pair.AddECFPname_withsurface()
            pair.Lextract = False

        print ('allpair read : %d'%len(self))

    def Loadpkl_all(self, ECFP0file=False, ECFPwithsitefile=False, RPDict0=False, RPDict=False,\
                    RPDict0_count=False, RPDict0_rangeinfo=False, RPDictwithexample=False):
        if ECFP0file :
            inputfile = open(ECFP0file,'rb')
            self.ECFP0dict = pickle.load(inputfile)
            inputfile.close()

        if ECFPwithsitefile:
            inputfile = open(ECFPwithsitefile,'rb')
            self.ECFPwithsitedict = pickle.load(inputfile)
            inputfile.close()

        if RPDict0:
            inputfile = open(RPDict0,'rb')
            self.allrpdict0 = pickle.load(inputfile)
            inputfile.close()

        if RPDict0_count:
            self.allrpdict0_count = Loadpkl(RPDict0_count)

        if RPDict0_rangeinfo:
            self.allrpdict0_rangeinfo = Loadpkl(RPDict0_rangeinfo)

        if RPDict:
            inputfile = open(RPDict,'rb')
            self.allrpdict_withsite = pickle.load(inputfile)
            inputfile.close()

        if RPDictwithexample:
            inputfile = open(RPDictwithexample,'rb')
            self.allrpdict_withsite_pairexample = pickle.load(inputfile)
            inputfile.close()


    def ECFP2NNinput(self,dict):
        RealIndexToInput = {}
        InputIndexToReal = {}

        fpall =[]
        for ECFP in dict.values():
            fpall.extend(ECFP.allindex)

        fpall.sort()
        fpall,fpcount =np.unique(fpall, return_counts=True)
        #fpall.sort()

        print fpall
        print fpcount
        for inputindex,index  in enumerate(fpall):
            #if fpcount[inputindex] <2 :
            #    continue
            RealIndexToInput[index] = inputindex
            InputIndexToReal[inputindex] = index
            
#        fpall =np.unique(fpall)
#        fpall.sort()
#        for inputindex,index  in enumerate(fpall):
#            RealIndexToInput[index] = inputindex
#            InputIndexToReal[inputindex] = index

        return RealIndexToInput,InputIndexToReal

    def AddECFPname_allbond(self,arcfile):
        _tmp=allstr_allbond()
        _tmp.readfile(arcfile)
        _tmp.GetTmpFakebmx()
        _tmp.GetAllECFPname()
        for i in range(0,len(_tmp),3):
            self[i/3][0].ECFPname_allbond = _tmp[i].ECFPname
            self[i/3][1].ECFPname_allbond = _tmp[i+2].ECFPname

    def AddSurfaceinfo_tmp(self):
        for pair in self:
            #print pair[0].label
            pair[0].surfaceinfo =   pair[0].label
            pair[1].surfaceinfo =   pair[1].label
            pair.surfaceinfo  =     pair[0].label

#        for pair in self:
#            if pair[0].natompe[-1] == 27:
#                pair[0].surfaceinfo = '111-27'
#                pair[1].surfaceinfo = '111-27'
#                pair.surfaceinfo = '111-27'
#
#            elif pair[0].natompe[-1] == 36:
#                pair[0].surfaceinfo = '100-36'
#                pair[1].surfaceinfo = '100-36'
#                pair.surfaceinfo = '100-36'
#
#            elif pair[0].natompe[-1] == 24:
#                pair[0].surfaceinfo = '211-24'
#                pair[1].surfaceinfo = '211-24'
#                pair.surfaceinfo = '211-24'
#
#            elif pair[0].natompe[-1] == 64:
#                if pair[0].atom[-1].xyz[-1] < 6 or \
#                    pair[0].atom[-2].xyz[-1] < 6 or \
#                    pair[0].atom[-3].xyz[-1] < 6 :
#                    pair[0].surfaceinfo = '211-64'
#                    pair[1].surfaceinfo = '211-64'
#                    pair.surfaceinfo = '211-64'
#                else:
#                    pair[0].surfaceinfo = '111-64'
#                    pair[1].surfaceinfo = '111-64'
#                    pair.surfaceinfo = '111-64'



    def GetAllECFPname_surface(self):
        self.AddSurfaceinfo_tmp()
        for pair in self:
            pair.AddECFPname_withsurface()

    def SelectPair(self):
        if not hasattr(self,'allrpdict0_rangeinfo'):
            print ('Gen RP Dict')
            self.AllRPDict()
            self.AllRPDict0()

        new = AllPair()
        possiblewrong = AllPair()
        delete = AllPair()
        for ipair,pair in enumerate(self):
            if pair.barrier1 > 2.5 or pair.barrier2 > 2.5: 
                delete.append(pair)
                continue
            if pair.barrier1 < 0 or pair.barrier2 < 0: 
                delete.append(pair)
                continue
            if ((pair.barrier1 - self.allrpdict0_rangeinfo[pair.rp1.ECFP0name].barrieravg) > 1.5 or \
                (pair.barrier2 - self.allrpdict0_rangeinfo[pair.rp2.ECFP0name].barrieravg) > 1.5 or \
                (pair.barrier1 - self.allrpdict0_rangeinfo[pair.rp1.ECFP0name].barriermin) > 2.0 or \
                (pair.barrier2 - self.allrpdict0_rangeinfo[pair.rp2.ECFP0name].barriermin) > 2.0) :
                print ('too high barrier: Pair %d Rp0: %s Min Barrier: %8.4f  Avg barrier: %8.4f  This Barrrier: %8.4f'\
                    %(ipair,pair.rp1.ECFP0name,self.allrpdict0_rangeinfo[pair.rp1.ECFP0name].barriermin, \
                        self.allrpdict0_rangeinfo[pair.rp1.ECFP0name].barrieravg, pair.barrier1))
                delete.append(pair)
                continue

            #print self.allrpdict0_rangeinfo[pair.rp1.ECFP0name].barrieravg
            if (abs(pair.barrier1 - self.allrpdict0_rangeinfo[pair.rp1.ECFP0name].barrieravg) > 0.5) or \
               (abs(pair.heat1- self.allrpdict0_rangeinfo[pair.rp1.ECFP0name].heatavg) > 0.5):
                possiblewrong.append(pair)
                new.append(pair)
                continue
            #print self.allrpdict0_rangeinfo[pair.rp2.ECFP0name].barrieravg
            if (abs(pair.barrier2 - self.allrpdict0_rangeinfo[pair.rp2.ECFP0name].barrieravg) > 0.5) or \
               (abs(pair.heat2- self.allrpdict0_rangeinfo[pair.rp2.ECFP0name].heatavg) > 0.5):
                possiblewrong.append(pair)
                new.append(pair)
                continue
            if (self.allrpdict0_count[pair.rp1.ECFP0name]) >= 3 or\
                (self.allrpdict_withsite[pair.rp1.ECFPname_withsite].nsample >= 2) or\
                (pair.barrier1 < 1.5 and pair.barrier2 < 1.5):
                new.append(pair)
            else:
                delete.append(pair)
        possiblewrong.Arcprint('possiblewrongTS.arc')
        delete.Arcprint('delete.arc')
        return new

    def SelectPattern(self):
        if not hasattr(self,'allrpdict0_rangeinfo'):
            print ('Gen RP Dict')
            self.AllRPDict()
            self.AllRPDict0()

        new = AllPair()
        AddList = []
        for rp0 in self.allrpdict0.keys():
            for rpwithsite in self.allrpdict0[rp0][:25]:
                for pair in self.allrpdict_withsite_pairexample[rpwithsite]:
                    if pair.ECFPname_withsite not in AddList:
                        new.append(pair)
                        AddList.append(pair.ECFPname_withsite)


        return new

    def Combine_Select(self):
        if not hasattr(self,'selectpair_all'):
            self.AllPair_label()
          
        if not hasattr(self,'allrpdict0_rangeinfo'):
            print ('Gen RP Dict')
            self.AllRPDict()
            self.AllRPDict0()

        new = AllPair()
        AddList = []
        for rp0 in self.allrpdict0.keys():
            for rpwithsite in self.allrpdict0[rp0][:15]:
                for pair in self.allrpdict_withsite_pairexample[rpwithsite]:
                    if pair.ECFPname_withsite not in AddList:
                        new.append(pair)
                        AddList.append(pair.ECFPname_withsite)

        pw = AllPair()
        for label in self.selectpair_all.keys():
            if self.selectpair_all[label][0].ECFPname_withsite not in AddList:
                new.append(self.selectpair_all[label][0])
                pw.append(self.selectpair_all[label][0])
        pw.Arcprint('strangelowpair.arc')

        for i in range(100): random.shuffle(new)
        #outnew= AllPair()
        #outnew= Shuffle_listlike(new,outnew) 
        return new


    def AllPair_label(self):
        self.selectpair_all = {}

        for pair in self:
            pair.label = pair.ECFPname_surface + pair.rp1.ECFP0name
            if pair.label not in self.selectpair_all.keys():
                self.selectpair_all[pair.label]= []
                self.selectpair_all[pair.label].append(pair)
            else:
                self.selectpair_all[pair.label].append(pair)

        new = AllPair()
        for label in self.selectpair_all.keys():
            self.selectpair_all[label].sort(key = lambda x:x.TS)
        return 


    def SelectNew(self):
        self.AllPair_label()
        new = AllPair()
        for label in self.selectpair_all.keys():
            for pair in self.selectpair_all[label][:3]:
                new.append(pair)

        return new



    def GenInputVector(self,RealIndexToInput):
        dataInput= {}
        _tmp =allstr_ecfp0()
        for pair in self:
#            if pair.barrier1 > 2.5 or pair.barrier2 > 2.5: continue
#            if self.allrpdict0_count[pair.rp1.ECFP0name] > 2 or (pair.barrier1 < 1.5 and pair.barrier2 < 1.5):
            _tmp.append(pair[0])
            _tmp.append(pair[1])

        for i in range(len(_tmp)):
            str1 = _tmp[i]
            fp= _tmp[i].ECFP
            inputvector =  np.zeros((len(RealIndexToInput),), dtype = np.int)
            for key in fp.allindex:
                try:
                    inputvector[RealIndexToInput[key]] = len(fp.IndextoFp[key])
                except:
                    print ('warning: contains unsupported index: %d'%key)
            dataInput[str1.name]= inputvector
        
        return dataInput

    def AllMinchange_ECFPsurface(self,minfile):

        _tmp= allstr_ecfp0()
        _tmp.readfile(minfile)
        _tmp.GetTmpFakebmx()
        _tmp.GetAllECFPname()
        for Str in _tmp:
            Str.surfaceinfo = Str.label
#            if Str.natompe[-1] == 27:
#                Str.surfaceinfo = '111-27'
#            elif Str.natompe[-1] == 36:
#                Str.surfaceinfo = '100-36'
#            elif Str.natompe[-1] == 24:
#                Str.surfaceinfo = '211-24'
#            elif Str.natompe[-1] == 64:
#                if Str.atom[-1].xyz[-1] <  6 or \
#                    Str.atom[-2].xyz[-1] < 6 or \
#                    Str.atom[-3].xyz[-1] < 6:
#                    Str.surfaceinfo = '211-64'
#                else:
#                    Str.surfaceinfo = '111-64'
            Str.ECFPname_surface = Str.ECFPname + '-'  + Str.surfaceinfo

        Mindict= {}
        for Str in _tmp:
            if Str.ECFPname_surface not in Mindict.keys():
                Mindict[Str.ECFPname_surface]= Str
            else:
                if Mindict[Str.ECFPname_surface].energy > Str.energy:
                    Mindict[Str.ECFPname_surface] = Str
        

        out = allstr_ecfp0()
        for key in Mindict.keys():
            Mindict[key].label = Mindict[key].ECFPname_surface
            out.append(Mindict[key])
        out.printall('allname.arc_ECFPsurface')
        return

    def GenSampledsict(self, inputfile):
        _tmp= allstr_ecfp0()
        _tmp.readfile(inputfile)
        _tmp.GetTmpFakebmx(numproc = 24)
        _tmp.GetAllECFPname(numproc = 24)
        for Str in _tmp:
            Str.surfaceinfo = Str.label
            Str.ECFPname_surface = Str.ECFPname + '-'  + Str.surfaceinfo

        Sampleddict= {}
        for key in self.mindict_surface.keys():
            Sampleddict[key] = 0
        for Str in _tmp:
            #self.mindict_surface[Str.ECFPname_surface] = Str.energy
            Sampleddict[Str.ECFPname_surface]= 1
        self.sampleddict = Sampleddict

        Dumppkl('sampleddict.pkl',self.sampleddict)
        del _tmp
        return



    def GenMindict(self,minfile):

        _tmp= allstr_ecfp0()
        _tmp.readfile(minfile)
        _tmp.GetTmpFakebmx(numproc = 3)
        _tmp.GetAllECFPname(numproc = 3)
        for Str in _tmp:
            Str.surfaceinfo = Str.label
#            if Str.natompe[-1] == 27:
#                Str.surfaceinfo = '111-27'
#            elif Str.natompe[-1] == 36:
#                Str.surfaceinfo = '100-36'
#            elif Str.natompe[-1] == 24:
#                Str.surfaceinfo = '211-24'
#            elif Str.natompe[-1] == 64:
#                if Str.atom[-1].xyz[-1] <  6 or \
#                    Str.atom[-2].xyz[-1] < 6 or \
#                    Str.atom[-3].xyz[-1] < 6:
#                    Str.surfaceinfo = '211-64'
#                else:
#                    Str.surfaceinfo = '111-64'
#
            Str.ECFPname_surface = Str.ECFPname + '-'  + Str.surfaceinfo 

#        self.GetAllECFPname_surface()


        Mindict= {}
        for Str in _tmp:
            if Str.ECFPname_surface not in Mindict.keys():
                Mindict[Str.ECFPname_surface]= Str.energy
            else:
                Mindict[Str.ECFPname_surface] = min(Mindict[Str.ECFPname_surface],Str.energy)
        self.mindict_surface = Mindict
        del _tmp
        return 

    def ReadDict(self,filename):
        dict = {}
        f = open(filename, 'r')
        for l in f.readlines():
            line = l.strip("\n").split()
            dict[line[0]] = float(line[1])
        f.close()
        return dict


    def GetMinNamebyCombine(self,minfile):
        surfacedict = self.ReadDict('surfacedict')
        fragdict = self.ReadDict('fragdict')
        print surfacedict
        _tmp= allstr_ecfp0()
        _tmp.readfile(minfile)
        _tmp.GetTmpFakebmx(numproc = 12)
        _tmp.GetAllECFPname(numproc = 12)
        _tmp.GetAllFragECFPname()

        out = allstr_ecfp0()
        out2 = allstr_ecfp0()
        for Str in _tmp:
            Str.surfaceinfo = Str.label.split('-')[-2] + '-' + Str.label.split('-')[-1]
            #print Str.surfaceinfo
            combineenergy = 0
            try:
                for frag in Str.AllFragECFPname:
                    fragname = frag+ '-' + Str.surfaceinfo
                    #print fragname    
                    combineenergy = combineenergy + fragdict[fragname] - surfacedict[Str.surfaceinfo]
                combineenergy = combineenergy + surfacedict[Str.surfaceinfo]
            except:
                out.append(Str)
            print Str.energy, combineenergy
            if combineenergy < Str.energy :
                Str.energy = combineenergy
            else:
                print ('possible inter-mol interact')
                out2.append(Str)
            Str.ECFPname_surface = Str.ECFPname + '-'  + Str.surfaceinfo

        out.printall('pwmin.arc')
        out2.printall('possibleinteract.arc')


        Mindict= {}
        for Str in _tmp:
            if Str.ECFPname_surface not in Mindict.keys():
                Mindict[Str.ECFPname_surface]= Str.energy
            else:
                Mindict[Str.ECFPname_surface] = min(Mindict[Str.ECFPname_surface],Str.energy)
        self.mindict_surface = Mindict

        print 'self.mindict_surface'
        print len(self.mindict_surface)
        _tmp.printall('allname.arc_finalused')
        del _tmp
        return


    def AddRPenergyData(self):
        Mindict = self.mindict_surface
        for pair in self:
            try:
                barrier= pair.TSstr.energy - Mindict[pair[0].ECFPname_surface]
                if barrier <0 or barrier > 10:
                    print 'wrongbarrier', pair.TSstr.energy, Mindict[pair[0].ECFPname_surface], pair[0].energy
                try:
                    heat = Mindict[pair[1].ECFPname_surface] -Mindict[pair[0].ECFPname_surface]
                except:
                    print pair[1].ECFPname_surface,pair[0].ECFPname_surface
                pair.barrier1 = barrier
                pair.heat1 = heat
                pair.rp1.barrier = barrier
                pair.rp1.heat = heat
#                if pair.rp1.barriermax < barrier :
#                    pair.rp1.barriermax == barrier
#                if pair.rp1.barriermin > barrier:
#                    pair.rp1.barriermin == barrier


                barrier= pair.TSstr.energy - Mindict[pair[1].ECFPname_surface]
                if barrier <0 or barrier > 10:
                    print 'wrongbarrier', pair.TSstr.energy, Mindict[pair[1].ECFPname_surface], pair[1].energy
                heat = Mindict[pair[0].ECFPname_surface] -Mindict[pair[1].ECFPname_surface]
                pair.barrier2 = barrier
                pair.heat2 = heat
                pair.rp2.barrier = barrier
                pair.rp2.heat = heat

            except:
                print ('problem pair: %s %s'%(pair[1].ECFPname_surface,pair[0].ECFPname_surface))
                barrier = 999
                heat = 0
                pair.barrier1 = barrier
                pair.heat1 = heat
                pair.rp1.barrier = barrier
                pair.rp1.heat = heat
                pair.barrier2 = barrier
                pair.heat2 = heat
                pair.rp2.barrier = barrier
                pair.rp2.heat = heat




    def GetAllrpECFPdict(self):

        fpall= []
        for rp in self.allrpdict_withsite.values():
            fpall.extend(rp.ECFP_withsite)
            #print "rp.ECFP_withsite"
            #print rp.ECFP_withsite


        RealIndexToInput = {}
        InputIndexToReal = {}

        fpall.sort()
        fpall,fpcount =np.unique(fpall, return_counts=True)
        #fpall.sort()

        print fpall
        print fpcount
        for inputindex,index  in enumerate(fpall):
            #if fpcount[inputindex] <2 :
            #    continue
            RealIndexToInput[index] = inputindex
            InputIndexToReal[inputindex] = index


        return RealIndexToInput,InputIndexToReal


    def GetAllrpECFPdict_ECFP0(self):

        fpall= []
        for rp in self.allrpdict_withsite.values():
            fpall.extend(rp.ECFP0)
            #print "rp.ECFP_withsite"
            #print rp.ECFP_withsite


        RealIndexToInput = {}
        InputIndexToReal = {}

        fpall.sort()
        fpall,fpcount =np.unique(fpall, return_counts=True)
        #fpall.sort()

        print fpall
        print fpcount
        for inputindex,index  in enumerate(fpall):
            #if fpcount[inputindex] <2 :
            #    continue
            RealIndexToInput[index] = inputindex
            InputIndexToReal[inputindex] = index


        return RealIndexToInput,InputIndexToReal


    def GenRPInputVector(self,RealIndexToInput):

        dataInput= {}
        for rp in self.allrpdict_withsite.values():
            inputvector =  np.zeros((len(RealIndexToInput),), dtype = np.int)
            for key in rp.ECFP_withsite:
                try:
                    inputvector[RealIndexToInput[key]] = inputvector[RealIndexToInput[key]] + 1
                except:
                    print ('warning: contains unsupported index: %d'%key)
            dataInput[rp.ECFPname_withsite]= inputvector
        return dataInput

    def GenRPInputVector_ECFP0(self,RealIndexToInput):

        dataInput= {}
        for rpl in self.allrpdict0_all.values():
            rp = rpl[0]
            inputvector =  np.zeros((len(RealIndexToInput),), dtype = np.int)
            for key in rp.ECFP0:
                try:
                    inputvector[RealIndexToInput[key]] = inputvector[RealIndexToInput[key]] + 1
                except:
                    print ('warning: contains unsupported index: %d'%key)
            dataInput[rp.ECFP0name]= inputvector
        return dataInput


    def ScreenUpperPair(self):
        for pair in self:
            pair[0].screenuppersurf()
            pair[1].screenuppersurf()
            if pair[0].upper ==1 or pair[1].upper ==1:
                continue

    def GenInfoNetout_ECFP0(self, selectsurface = 0,  Laddsurface = 0):
        RealIndexToInput_rp,InputIndexToReal_rp= self.GetAllrpECFPdict_ECFP0()
        rpdataInput = self.GenRPInputVector_ECFP0(RealIndexToInput_rp)

        RealIndexToInput_str,InputIndexToReal_str = self.ECFP2NNinput(self.ECFP0dict)
        strdataInput= self.GenInputVector(RealIndexToInput_str)


        dataInput2 = []
        dataInput3 = []
        dataOutput2 = []
        dataLabel = []

        for pair in self:
            s = [0,0,0]
            _surf = int(pair.surfaceinfo.split('-')[0])
            if _surf == 111:
                s[0] = 1
            elif _surf == 100:
                s[1] = 1
            elif _surf == 211:
                s[2] = 1

            #print pair.surfaceinfo.split('-')[0]
            if ((selectsurface != 0) and (selectsurface != int(pair.surfaceinfo.split('-')[0]))):  continue

            val1 = strdataInput[pair[0].ECFP0name]
            val2 = rpdataInput[pair.rp1.ECFP0name]
            val = np.append(val1,val2)
            dataInput2.append(val)

            _rangeinfo = self.allrpdict0_rangeinfo[pair.rp1.ECFP0name]
            if Laddsurface ==1:
                _val= [_rangeinfo.barriermin, _rangeinfo.barriermax, _rangeinfo.barrieravg, _rangeinfo.heatmax, _rangeinfo.heatmin, _rangeinfo.heatavg,\
                         s[0],s[1],s[2]]
            else:
                _val= [_rangeinfo.barriermin, _rangeinfo.barriermax, _rangeinfo.barrieravg, _rangeinfo.heatmax, _rangeinfo.heatmin, _rangeinfo.heatavg]
            #_val= [_rangeinfo.barriermin, _rangeinfo.barriermax, _rangeinfo.barrieravg, _rangeinfo.heatmax, _rangeinfo.heatmin, _rangeinfo.heatavg]
            dataInput3.append(_val)
            dataOutput2.append([pair.barrier1, pair.heat1])

            if pair[0].lsampled == 1  and pair[1].lsampled == 1:
                labeled = [10]
            else :
                labeled = [1]
            dataLabel.append(labeled)

            val1 = strdataInput[pair[1].ECFP0name]
            val2 = rpdataInput[pair.rp2.ECFP0name]
            val = np.append(val1,val2)
            dataInput2.append(val)
            _rangeinfo = self.allrpdict0_rangeinfo[pair.rp2.ECFP0name]
            if Laddsurface ==1:
                _val= [_rangeinfo.barriermin, _rangeinfo.barriermax, _rangeinfo.barrieravg, _rangeinfo.heatmax, _rangeinfo.heatmin, _rangeinfo.heatavg,\
                         s[0],s[1],s[2]]
            else:
                _val= [_rangeinfo.barriermin, _rangeinfo.barriermax, _rangeinfo.barrieravg, _rangeinfo.heatmax, _rangeinfo.heatmin, _rangeinfo.heatavg]
            #_val= [_rangeinfo.barriermin, _rangeinfo.barriermax, _rangeinfo.barrieravg, _rangeinfo.heatmax, _rangeinfo.heatmin, _rangeinfo.heatavg]
            dataInput3.append(_val)
            dataOutput2.append([pair.barrier2, pair.heat2])

            dataLabel.append(labeled)
        Dumppkl('RealIndextoInput_str_%s.pkl'%(selectsurface),RealIndexToInput_str)
        Dumppkl('RealIndextoInput_rp_%s.pkl'%(selectsurface),RealIndexToInput_rp)

        return dataInput2, dataOutput2,dataInput3, dataLabel


    def GenInfoNetout_ECFP0_TwoStr(self, selectsurface = 0, Laddsurface = 0):
        #RealIndexToInput_rp,InputIndexToReal_rp= self.GetAllrpECFPdict_ECFP0()
        #rpdataInput = self.GenRPInputVector_ECFP0(RealIndexToInput_rp)

        RealIndexToInput_str,InputIndexToReal_str = self.ECFP2NNinput(self.ECFP0dict)
        strdataInput= self.GenInputVector(RealIndexToInput_str)
    
        


        dataInput2 = []
        dataInput3 = []
        dataOutput2 = []
        dataLabel = []

        for pair in self:

            s = [0,0,0]
            _surf = int(pair.surfaceinfo.split('-')[0])
            if _surf == 111:
                s[0] = 1
            elif _surf == 100:
                s[1] = 1
            elif _surf == 211:
                s[2] = 1

            #print pair.surfaceinfo.split('-')[0]
            if ((selectsurface != 0) and (selectsurface != _surf)):  continue

            val1 = strdataInput[pair[0].ECFP0name]
            val2 = strdataInput[pair[1].ECFP0name]
            val = np.append(val1,val2)
            dataInput2.append(val)

            _rangeinfo = self.allrpdict0_rangeinfo[pair.rp1.ECFP0name]
            if Laddsurface ==1:
                _val= [_rangeinfo.barriermin, _rangeinfo.barriermax, _rangeinfo.barrieravg, _rangeinfo.heatmax, _rangeinfo.heatmin, _rangeinfo.heatavg,\
                         s[0],s[1],s[2]]
            else:
                _val= [_rangeinfo.barriermin, _rangeinfo.barriermax, _rangeinfo.barrieravg, _rangeinfo.heatmax, _rangeinfo.heatmin, _rangeinfo.heatavg]
            dataInput3.append(_val)
            dataOutput2.append([pair.barrier1, pair.heat1])

            if pair[0].lsampled == 1  and pair[1].lsampled == 1:
                labeled = [10]
            else :
                labeled = [1]
            dataLabel.append(labeled)

            val1 = strdataInput[pair[1].ECFP0name]
            val2 = strdataInput[pair[0].ECFP0name]
            val = np.append(val1,val2)
            dataInput2.append(val)
            _rangeinfo = self.allrpdict0_rangeinfo[pair.rp2.ECFP0name]
            if Laddsurface ==1:
                _val= [_rangeinfo.barriermin, _rangeinfo.barriermax, _rangeinfo.barrieravg, _rangeinfo.heatmax, _rangeinfo.heatmin, _rangeinfo.heatavg,\
                         s[0],s[1],s[2]]
            else:
                _val= [_rangeinfo.barriermin, _rangeinfo.barriermax, _rangeinfo.barrieravg, _rangeinfo.heatmax, _rangeinfo.heatmin, _rangeinfo.heatavg]
            dataInput3.append(_val)
            dataOutput2.append([pair.barrier2, pair.heat2])

            dataLabel.append(labeled)
        Dumppkl('RealIndextoInput_str_%s.pkl'%(selectsurface),RealIndexToInput_str)
        #Dumppkl('RealIndextoInput_rp_%s.pkl'%(selectsurface),RealIndexToInput_rp)

        return dataInput2, dataOutput2,dataInput3, dataLabel






    def GenInfoNetout(self):
        RealIndexToInput_rp,InputIndexToReal_rp= self.GetAllrpECFPdict()
        rpdataInput = self.GenRPInputVector(RealIndexToInput_rp)


        RealIndexToInput_str,InputIndexToReal_str = self.ECFP2NNinput(self.ECFPwithsitedict)
        strdataInput= self.GenInputVector(RealIndexToInput_str)


        dataInput2 = []
        dataInput3 = []
        dataOutput2 = []
        dataLabel = []

        for pair in self:
            val1 = strdataInput[pair[0].ECFPname_withsite]
            val2 = rpdataInput[pair.rp1.ECFPname_withsite]
            val = np.append(val1,val2)
            dataInput2.append(val)

            _rangeinfo = self.allrpdict_withsite_rangeinfo[pair.rp1.ECFPname_withsite]
            _val= [_rangeinfo.barriermin, _rangeinfo.barriermax, _rangeinfo.barrieravg, _rangeinfo.heatmax, _rangeinfo.heatmin, _rangeinfo.heatavg]
            dataInput3.append(_val)
            dataOutput2.append([pair.barrier1, pair.heat1])

            if pair[0].lsampled == 1  and pair[1].lsampled == 1:
                labeled = [1]
            else :
                labeled = [0.1]
            dataLabel.append(labeled)

            val1 = strdataInput[pair[1].ECFPname_withsite]
            val2 = rpdataInput[pair.rp2.ECFPname_withsite]
            val = np.append(val1,val2)
            dataInput2.append(val)
            _rangeinfo = self.allrpdict_withsite_rangeinfo[pair.rp2.ECFPname_withsite]
            _val= [_rangeinfo.barriermin, _rangeinfo.barriermax, _rangeinfo.barrieravg, _rangeinfo.heatmax, _rangeinfo.heatmin, _rangeinfo.heatavg]
            dataInput3.append(_val)
            dataOutput2.append([pair.barrier2, pair.heat2])

            dataLabel.append(labeled)

            #del pair.rp1
            #del pair.rp2

        #input2 = np.array(dataInput2)
        #result2 = np.array(dataOutput2)

        Dumppkl('RealIndextoInput_strwithsite.pkl',RealIndexToInput_str)
        Dumppkl('RealIndextoInput_rp.pkl',RealIndexToInput_rp)

        return dataInput2, dataOutput2,dataInput3, dataLabel




    def CheckAllRP(self):
        PattDict = {}
        for pair in self:
            if pair.rp1.ECFP0name not in PattDict.keys():
                iPattern= len(PattDict)
                PattDict[pair.rp1.ECFP0name] = iPattern
            if pair.rp2.ECFP0name not in PattDict.keys():
                iPattern= len(PattDict)
                PattDict[pair.rp2.ECFP0name] = iPattern
        print len(PattDict)

        _tmp= AllPair()
        for name in self.allrpdict0.keys():
            if name not in PattDict.keys():
                _nopair =self.allrpdict_withsite_pairexample[self.allrpdict0[name][0]][0]
                _tmp.append(_nopair)

        _tmp.Arcprint('pattern.arc')
        return


    def GenPolicyNetout(self):

        RealIndexToInput,InputIndexToReal = self.ECFP2NNinput(self.ECFP0dict)
        dataInput= self.GenInputVector(RealIndexToInput)


        PattDict = {}
        PattIndexDict = {}
        PattDict_Str ={}

        for pair in self:
            if pair.barrier1 < 2.5:
                fakek = GetFakek(pair.barrier1)
                if pair.rp1.ECFP0name not in PattDict.keys():
                    iPattern= len(PattDict)
                    PattDict[pair.rp1.ECFP0name] = iPattern
                    PattIndexDict[iPattern] = pair.rp1
    
                iPattern = PattDict[pair.rp1.ECFP0name] 
                if pair[0].ECFP0name in PattDict_Str.keys():
                    if iPattern not in PattDict_Str[pair[0].ECFP0name].keys():
                        PattDict_Str[pair[0].ECFP0name][iPattern] = fakek
                    elif fakek > PattDict_Str[pair[0].ECFP0name][iPattern]: 
                        PattDict_Str[pair[0].ECFP0name][iPattern]= fakek
                else:
                    PattDict_Str[pair[0].ECFP0name] = {}
                    PattDict_Str[pair[0].ECFP0name][iPattern] =fakek


            if pair.barrier2 < 2.5:            
                fakek = GetFakek(pair.barrier2)
                if pair.rp2.ECFP0name not in PattDict.keys():
                    iPattern= len(PattDict)
                    PattDict[pair.rp2.ECFP0name] = iPattern
                    PattIndexDict[iPattern] = pair.rp2
    
                iPattern = PattDict[pair.rp2.ECFP0name]
                if pair[1].ECFP0name in PattDict_Str.keys():
                    if iPattern not in PattDict_Str[pair[1].ECFP0name].keys():
                        PattDict_Str[pair[1].ECFP0name][iPattern] = fakek
                    elif fakek > PattDict_Str[pair[1].ECFP0name][iPattern] :
                        PattDict_Str[pair[1].ECFP0name][iPattern]= fakek
                else:
                    PattDict_Str[pair[1].ECFP0name] = {}
                    PattDict_Str[pair[1].ECFP0name][iPattern] =fakek

        allpatt=len(PattDict)
        dataOutput={}
            
        for name in PattDict_Str.keys():
            dataOutput[name]= np.zeros(allpatt)
            all = 0
            for key,value in  PattDict_Str[name].items() :
                #dataOutput[name][key] = value
                #all= all +np.exp(value)
                all= all + value
            for key,value in  PattDict_Str[name].items() :
                #dataOutput[name][key]= np.exp(value)/all
                inputp =  value/all
                if np.isnan(inputp):
                    print 'wrongpossibility',name,value,all
                    inputp =np.nan_to_num(inputp)
                dataOutput[name][key]= inputp

        datain = []
        dataresult = []
        
        namelist = dataOutput.keys()
        for i in range(200): random.shuffle(namelist)
        
        for name in namelist:
            #if name in wrongkeylist: continue
            datain.append(dataInput[name])
            dataresult.append(dataOutput[name])
        
        #print datain
        #print dataresult
        
        #input = np.array(datain)
        #result = np.array(dataresult)

        inputfile = open('stdoutput.pkl','wb')
        pickle.dump(dataOutput,inputfile)
        inputfile.close()


        #inputfile=open('datain.pkl','wb')
        #pickle.dump(input,inputfile)
        #inputfile.close()
        #
        #outputfile=open('dataresult.pkl','wb')
        #pickle.dump(result,outputfile)
        #outputfile.close()

        Dumppkl('PattIndexDict.pkl',PattIndexDict)
        Dumppkl('RealIndextoInput_0.pkl',RealIndexToInput)

        return  datain,dataresult 



if __name__ == "__main__":

    """init mode"""

    lready = CheckFile(['allgoodsect.arc','allname.arc'])
    if not lready:
        sys.exit()

    _pre = AllPair()
    _pre.PreScreen('allgoodsect.arc',filemode =2)
    del _pre

    test= AllPair()
    test.GenMindict('allname.arc')
    #test.GetMinNamebyCombine('allname.arc')
    test.GenSampledsict('sampled.arc')
    #test.readfile('allgoodsect.arc_out',filemode =2)
    test.ReadFile_parallel('allgoodsect.arc_out',filemode =2)
    test.GetAllECFPname_surface()
    #test.GenMindict('allname.arc')
    test.AddRPenergyData()



    """read data mode"""
#    lready = CheckFile(['ECFP0.pkl_save','ECFPwithsite.pkl_save','allrpdict0.pkl_save','allrpdict_withsite.pkl_save'\
#                      ,'allrpdict0_count.pkl_save','allrpdict0_rangeinfo.pkl_save','testdata_save'])
#    if not lready:
#        sys.exit()
#
#    new = AllPair()
#    #test.Loadpkl_all('ECFP0.pkl','ECFPwithsite.pkl','allrpdict0.pkl','allrpdict_withsite.pkl',RPDictwithexample='allrpdict_withsite_pairexample.pkl')
#    new.Loadpkl_all('ECFP0.pkl_save','ECFPwithsite.pkl_save','allrpdict0.pkl_save','allrpdict_withsite.pkl_save'\
#                      ,'allrpdict0_count.pkl_save','allrpdict0_rangeinfo.pkl_save',RPDictwithexample=False)
#    new.ReadRPData('testdata_save')
#




    """Pre Select"""
    new1 =test.SelectPair()
    new1.sampleddict = test.sampleddict
    del test    
    print ('SelectPair : %d' %(len(new1)))

    """Select low barrier pair"""
    #new = new1.SelectNew()
    new = new1.Combine_Select()
    new.sampleddict = new1.sampleddict
    del new1
    print ('SelectPair_lowbarrier : %d' %(len(new)))

    
    new = new.OutRPData('ECFP_surface','testdata')
    new.ChangeMode('ECFP0')
    input1, result1 = new.GenPolicyNetout()
    OutArray_Int(input1,'datain')
    OutArray_Float(result1,'dataresult')


#    new.ChangeMode('ECFP_withsite')
#    input2, result2, input_add , datalabel = new.GenInfoNetout()
#    OutArray_Int(input2,'datain_info')
#    OutArray_Float(result2,'dataresult_info')
#    OutArray_Float(input_add,'datain_addinfo')
#    #print datalabel
#    OutArray_Float(datalabel ,'datalabel')
    #Dumppkl('datainput_info.pkl',input2)
    #Dumppkl('dataout_info.pkl',result2)


    new.ChangeMode('ECFP0')
    input2, result2, input_add , datalabel = new.GenInfoNetout_ECFP0(Laddsurface = 1)
    #input2, result2, input_add , datalabel = new.GenInfoNetout_ECFP0_TwoStr(Laddsurface = 1)
    #input2, result2, input_add , datalabel = new.GenInfoNetout_ECFP0( selectsurface= 100)
    OutArray_Int(input2,'datain_info')
    OutArray_Float(result2,'dataresult_info')
    OutArray_Float(input_add,'datain_addinfo')
    OutArray_Float(datalabel ,'datalabel')


#    new.ChangeMode('ECFP0')
#    #input2, result2, input_add , datalabel = new.GenInfoNetout_ECFP0()
#    input2, result2, input_add , datalabel = new.GenInfoNetout_ECFP0_TwoStr(211)
#    #input2, result2, input_add , datalabel = new.GenInfoNetout_ECFP0( selectsurface= 100)
#    OutArray_Int(input2,'datain_info_211')
#    OutArray_Float(result2,'dataresult_info_211')
#    OutArray_Float(input_add,'datain_addinfo_211')
#    #print datalabel
#    OutArray_Float(datalabel ,'datalabel_211')
#
#
#
#    input2, result2, input_add , datalabel = new.GenInfoNetout_ECFP0_TwoStr(111)
#    OutArray_Int(input2,'datain_info_111')
#    OutArray_Float(result2,'dataresult_info_111')
#    OutArray_Float(input_add,'datain_addinfo_111')
#    OutArray_Float(datalabel ,'datalabel_111')
#
#
#    input2, result2, input_add , datalabel = new.GenInfoNetout_ECFP0_TwoStr(100)
#    OutArray_Int(input2,'datain_info_100')
#    OutArray_Float(result2,'dataresult_info_100')
#    OutArray_Float(input_add,'datain_addinfo_100')
#    OutArray_Float(datalabel ,'datalabel_100')

