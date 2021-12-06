import sys
sys.path.append('/home9/kpl/ReactNN/surface-new/newmode/pymodule_tmp/')

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
from reactpair_NN import Pair
from analyze_surfacemode_NNtmp import AllPair as allpair
from ECFP_allbond import AllStr as allstr_allbond
from simple_func_k import Loadpkl,Dumppkl,Average
import random

def GetFakek(Ea):
    #return 10000000000*np.exp(-Ea*96485/8.314/500)
    #return 10000000000*np.exp(-Ea*10)
    return np.exp(-Ea*10)
    #return 1

class RangeInfo():
    def __init__(self,rp):
        self.barriermin = rp.barrier
        self.barriermax = rp.barrier
        self.heatmin = rp.heat
        self.heatmax = rp.heat

    def Update(self, rp):
        if rp.barrier > self.barriermax:
            self.barriermax = rp.barrier
        if rp.barrier < self.barriermin:
            self.barriermin = rp.barrier
        if rp.heat > self.heatmax:
            self.heatmax = rp.heat
        if rp.heat < self.heatmin:
            self.heatmin = rp.heat


class AllPair(allpair):
    def __init__(self):
        list.__init__(self)
        self.Lallmin = 0

    def readfile(self,filename='allgoodsect.arc',filemode = 1,Lallmininput = 0,allminfile= 'allname.arc', pairsortflag = True):
        self.Lallmin = Lallmininput

        _tmp1= allstr_ecfp0()
        _tmp1.readfile(filename)
        _tmp1.GetAllsminame(numproc=12,colorflag = 1)
        _tmp1.GetAllECFPname(numproc=12)


        _tmp2= allstr_ecfpwithsite()
        _tmp2.readfile(filename)
        _tmp2.GetAllsminame(numproc=12,colorflag = 1)


        wronglist =allstr_ecfpwithsite()
        for i in range(0,len(_tmp1),3):

            if filemode == 1:
                _tmp1[i].screenuppersurf()
                _tmp1[i+1].screenuppersurf()
                if (_tmp1[i].upper ==1 ) or (_tmp1[i+1].upper == 1):  continue
                #if (not _tmp[i].Lminstr) or (not _tmp[i+1].Lminstr): continue
                if (not _tmp1[i].sminame) or (not _tmp1[i+1].sminame):  continue

                TS = _tmp1[i+2].energy
                p= Pair()
                p.AddStr(_tmp1[i],_tmp1[i+1],TS,_tmp1[i+2],pairsortflag)
                l = p.doublename(_tmp2[i],_tmp2[i+1])

            elif filemode ==2:

                _tmp1[i].screenuppersurf()
                _tmp1[i+2].screenuppersurf()
                if (_tmp1[i].upper ==1 ) or (_tmp1[i+2].upper == 1):  continue
                if (not _tmp1[i].sminame) or (not _tmp1[i+2].sminame):  continue

                TS = _tmp1[i+1].energy
                p= Pair()
                p.AddStr(_tmp1[i],_tmp1[i+2],TS,_tmp1[i+1],pairsortflag)
                l = p.doublename(_tmp2[i],_tmp2[i+2])

            if l == 1:
                self.append(p)
            else:
                print('wrong in Gen FP, pair %d'%(i/3))
                wronglist.append(_tmp1[i])
                wronglist.append(_tmp1[i+1])
                wronglist.append(_tmp1[i+2])

        wronglist.printall('wrongpair.arc')

    def ChangeMode(self,mode):
        for pair in self:
            pair.ChangeMode(mode)

    def ECFPwithsitedict(self):
        self.ECFPwithsitedict = {}
        for pair in self:
            if pair[0].ECFPname_withsite not in self.ECFPwithsitedict.keys():
                self.ECFPwithsitedict[pair[0].ECFPname_withsite] = pair[0].ECFP_withsite

            if pair[1].ECFPname_withsite not in self.ECFPwithsitedict.keys():
                self.ECFPwithsitedict[pair[1].ECFPname_withsite] = pair[1].ECFP_withsite
 
    def ECFP0dict(self):
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


        self.allrpdict0_avgbarrier = {}
        self.allrpdict0_avgheat = {}
        for key in self.allrpdict0_all.keys():
            self.allrpdict0_all[key].sort(key = lambda x: x.barrier)
            self.allrpdict0[key] = [rp.ECFPname_withsite for rp in self.allrpdict0_all[key]]
            self.allrpdict0_avgbarrier[key] = Average([rp.barrier for rp in self.allrpdict0_all[key]] )
            self.allrpdict0_avgheat[key] = Average([rp.heat for rp in self.allrpdict0_all[key]])


    def DumpECFP(self,ECFP0file = False,ECFPwithsitefile = False):
        if ECFP0file :
           self.ECFP0dict()
           outputfile = open(ECFP0file,'wb')
           pickle.dump(self.ECFP0dict,outputfile)
           outputfile.close()


        if ECFPwithsitefile:
           self.ECFPwithsitedict()
           outputfile = open(ECFPwithsitefile,'wb')
           pickle.dump(self.ECFPwithsitedict,outputfile)
           outputfile.close()



    def DumpAllRPDict(self,RPDict0=False,RPDict=False,RPDictwithexample=False):
        if RPDict0:
            print len(self.allrpdict0)
            outputfile = open(RPDict0,'wb')
            pickle.dump(self.allrpdict0,outputfile)
            outputfile.close()

        if RPDict:
            print len(self.allrpdict_withsite)
            outputfile = open(RPDict,'wb')
            pickle.dump(self.allrpdict_withsite,outputfile)
            outputfile.close()

        if RPDictwithexample:
            outputfile = open(RPDictwithexample,'wb')
            pickle.dump(self.allrpdict_withsite_pairexample,outputfile)
            outputfile.close()


    def OutRPData(self,mode,outfile):
        self.ChangeMode(mode)
        self.AllName()
        self.GenPairDict()
        #self.GetAllECFPname_surface()

        _tmp= AllPair()
        
        for pair in self.minpair.values():
            _tmp.append(pair)
        
        _tmp.OutDataSet(outfile)
        _tmp.Arcprint('allpair.arc')
        #self.AllRPDict()
        #self.AllRPDict0()
        #self.DumpAllRPDict('allrpdict0.pkl','allrpdict_withsite.pkl','allrpdict_withsite_pairexample.pkl')
        #self.DumpECFP('ECFP0.pkl','ECFPwithsite.pkl')
        return

    def Arcprint(self,outfile,list =False,Llabel =False):
        if not list:
            list = range(len(self))

        _tmp=allstr_ecfp0()
        for ipair in list:
            pair = self[ipair]
            pair[0].label = pair[0].ECFPname_surface 
            pair[1].label = pair[1].ECFPname_surface
            #pair.TSstr.label = 'RP_ECFP0name-'+'%s'%(pair.rp1.ECFP0name)+'-barrier-'+'%f'%(pair.barrier1)+'-heat-'+'%f'%(pair.heat1)
            _tmp.append(pair[0])
            _tmp.append(pair.TSstr)
            _tmp.append(pair[1])

        _tmp.printlist(range(len(_tmp)),outfile)

    def OutDataSet(self,outfile,list =False):
        if not list:
            list = range(len(self))

        f= open(outfile,'w')
        f.write('Raction Pair Dataset\n')
        f.write('tmptest for NN------format: kpl191222\n')
        f.write('Note: Only the chemisorption sites of reaction center are recorded and encoded in ECFPname_withste\n')
        f.write('\n')
        for ipair in list:
            pair = self[ipair]
            f.write('ReactPair:       %d\n'%(ipair+1))
            f.write('surfaceinfo:     %s\n'%pair.surfaceinfo)
            f.write('sminame of IS:   %s\n'%pair[0].sminame)
            f.write('sminame of FS:   %s\n'%pair[1].sminame)
            f.write('ECFPname_nosurface:       %s\n'%pair.ECFP0name)
            f.write('ECFPname_nosurface of IS: %s\n'%pair[0].ECFP0name)
            f.write('ECFPname_nosurface of FS: %s\n'%pair[1].ECFP0name)
            f.write('ECFPname_withsite:        %s\n'%pair.ECFPname_withsite)
            f.write('ECFPname_withsite of IS:  %s\n'%pair[0].ECFPname_withsite)
            f.write('ECFPname_withsite of FS:  %s\n'%pair[1].ECFPname_withsite)
            f.write('Reaction Center list 0 :  %s\n'%pair.rclist0)
            f.write('Reaction Center list 1 :  %s\n'%pair.rclist1)
            f.write('Reaction Pattern(IStoFS) ECFPname_nosurface: %s\n'%pair.rp1.ECFP0name)
            f.write('Reaction Pattern(IStoFS) ECFPname_withsite:  %s\n'%pair.rp1.ECFPname_withsite)
            f.write('Reaction Pattern(IStoFS) fpmatch link:       %s\n'%pair.rp1.fpmatch)
            f.write('Reaction Pattern(IStoFS) fpmatch atommap:    %s\n'%pair.rp1.atommatch)
            f.write('Reaction Pattern(FStoIS) ECFPname_nosurface: %s\n'%pair.rp2.ECFP0name)
            f.write('Reaction Pattern(FStoIS) ECFPname_withsite:  %s\n'%pair.rp2.ECFPname_withsite)
            f.write('Reaction Pattern(FStoIS) fpmatch link:       %s\n'%pair.rp2.fpmatch)
            f.write('Reaction Pattern(FStoIS) fpmatch atommap:    %s\n'%pair.rp2.atommatch)
            f.write('Sitematchpair of IS: %s\n'%pair[0].sitematch)
            f.write('Sitematchpair of FS: %s\n'%pair[1].sitematch)

            str=pair[0]
            #f.write('ISstr\n')
            f.write('     Energy     %8d    IS     %12.6f\n'%(ipair,str.energy))
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
            f.write('     Energy     %8d    FS     %12.6f\n'%(ipair,str.energy))
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
            f.write('     Energy     %8d    TS     %12.6f\n'%(ipair,str.energy))
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

    def ReadRPData(self,inputfile):
        _arcread =allstr_allbond()
        _arcread.readfile(inputfile)

        print len(_arcread)
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

    def Loadpkl(self, ECFP0file=False, ECFPwithsitefile=False, RPDict0=False, RPDict=False, RPDictwithexample=False):
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
            
        fpall =np.unique(fpall)
        fpall.sort()
        for inputindex,index  in enumerate(fpall):
            RealIndexToInput[index] = inputindex
            InputIndexToReal[inputindex] = index

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
            if pair[0].natompe[-1] == 27:
                pair[0].surfaceinfo = '111-27'
                pair[1].surfaceinfo = '111-27'
                pair.surfaceinfo = '111-27'

            elif pair[0].natompe[-1] == 36:
                pair[0].surfaceinfo = '110-36'
                pair[1].surfaceinfo = '110-36'
                pair.surfaceinfo = '110-36'

            elif pair[0].natompe[-1] == 24:
                pair[0].surfaceinfo = '211-24'
                pair[1].surfaceinfo = '211-24'
                pair.surfaceinfo = '211-24'

            elif pair[0].natompe[-1] == 64:
                pair[0].surfaceinfo = '211-64'
                pair[1].surfaceinfo = '211-64'
                pair.surfaceinfo = '211-64'

    def GetAllECFPname_surface(self):
        self.AddSurfaceinfo_tmp()
        for pair in self:
            pair.AddECFPname_withsurface()

    def SelectPair(self):
        self.AllRPDict()
        self.AllRPDict0()

        new = AllPair()
        possiblewrong = AllPair()
        for pair in self:
            if pair.barrier1 > 2.5 or pair.barrier2 > 2.5: continue
            if (abs(pair.barrier1 - self.allrpdict0_avgbarrier[pair.rp1.ECFP0name]) > 0.5) or (abs(pair.heat1- self.allrpdict0_avgheat[pair.rp1.ECFP0name]) > 0.5):
                possiblewrong.append(pair)
                continue
            if (abs(pair.barrier2 - self.allrpdict0_avgbarrier[pair.rp2.ECFP0name]) > 0.5) or (abs(pair.heat2- self.allrpdict0_avgheat[pair.rp2.ECFP0name]) > 0.5):
                possiblewrong.append(pair)
                continue
            if self.allrpdict0_count[pair.rp1.ECFP0name] > 2 or (pair.barrier1 < 1.5 and pair.barrier2 < 1.5):
                new.append(pair)
        possiblewrong.Arcprint('possiblewrongTS.arc')
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
                inputvector[RealIndexToInput[key]] = len(fp.IndextoFp[key])
            dataInput[str1.name]= inputvector
        
        return dataInput

    def AllMinchange_ECFPsurface(self,minfile):

        _tmp= allstr_ecfp0()
        _tmp.readfile(minfile)
        _tmp.GetTmpFakebmx()
        _tmp.GetAllECFPname()
        for Str in _tmp:
            if Str.natompe[-1] == 27:
                Str.surfaceinfo = '111-27'
            elif Str.natompe[-1] == 36:
                Str.surfaceinfo = '110-36'
            elif Str.natompe[-1] == 24:
                Str.surfaceinfo = '211-24'
            elif Str.natompe[-1] == 64:
                Str.surfaceinfo = '211-64'
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


    def GenMindict(self,minfile):
#        self.AddECFPname_allbond('allpair.arc')

#        _tmp = allstr_allbond()
#        _tmp.readfile(minfile)
#        _tmp.GetTmpFakebmx()
#        _tmp.GetAllECFPname()

        _tmp= allstr_ecfp0()
        _tmp.readfile(minfile)
        _tmp.GetTmpFakebmx()
        _tmp.GetAllECFPname()
        for Str in _tmp:
            if Str.natompe[-1] == 27:
                Str.surfaceinfo = '111-27'
            elif Str.natompe[-1] == 36:
                Str.surfaceinfo = '110-36'
            elif Str.natompe[-1] == 24:
                Str.surfaceinfo = '211-24'
            elif Str.natompe[-1] == 64:
                Str.surfaceinfo = '211-64'
            Str.ECFPname_surface = Str.ECFPname + '-'  + Str.surfaceinfo 

#        self.GetAllECFPname_surface()


        Mindict= {}
        for Str in _tmp:
            if Str.ECFPname_surface not in Mindict.keys():
                Mindict[Str.ECFPname_surface]= Str.energy
            else:
                Mindict[Str.ECFPname_surface] = min(Mindict[Str.ECFPname_surface],Str.energy)
        self.mindict_surface = Mindict
        return 

    def AddRPenergyData(self):
        Mindict = self.mindict_surface
        for pair in self:
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
#            if pair.rp1.barriermax < barrier :
#                pair.rp1.barriermax == barrier
#            if pair.rp1.barriermin > barrier:
#                pair.rp1.barriermin == barrier


            barrier= pair.TSstr.energy - Mindict[pair[1].ECFPname_surface]
            if barrier <0 or barrier > 10:
                print 'wrongbarrier', pair.TSstr.energy, Mindict[pair[1].ECFPname_surface], pair[1].energy
            heat = Mindict[pair[0].ECFPname_surface] -Mindict[pair[1].ECFPname_surface]
            pair.barrier2 = barrier
            pair.heat2 = heat
            pair.rp2.barrier = barrier
            pair.rp2.heat = heat




    def GetAllrpECFPdict(self):

        fpall= []
        for rp in self.allrpdict_withsite.values():
            fpall.extend(rp.ECFP_withsite)

        RealIndexToInput = {}
        InputIndexToReal = {}


        fpall =np.unique(fpall)
        fpall.sort()

        for inputindex,index  in enumerate(fpall):
            RealIndexToInput[index] = inputindex
            InputIndexToReal[inputindex] = index


        return RealIndexToInput,InputIndexToReal

    def GenRPInputVector(self,RealIndexToInput):

        dataInput= {}
        for rp in self.allrpdict_withsite.values():
            inputvector =  np.zeros((len(RealIndexToInput),), dtype = np.int)
            for key in rp.ECFP_withsite:
                inputvector[RealIndexToInput[key]] = inputvector[RealIndexToInput[key]] + 1
            dataInput[rp.ECFPname_withsite]= inputvector
        return dataInput

    def ScreenUpperPair(self):
        for pair in self:
            pair[0].screenuppersurf()
            pair[1].screenuppersurf()
            if pair[0].upper ==1 or pair[1].upper ==1:
                continue

    def GenInfoNetout(self):
        RealIndexToInput_rp,InputIndexToReal_rp= self.GetAllrpECFPdict()
        rpdataInput = self.GenRPInputVector(RealIndexToInput_rp)


        RealIndexToInput_str,InputIndexToReal_str = self.ECFP2NNinput(self.ECFPwithsitedict)
        strdataInput= self.GenInputVector(RealIndexToInput_str)


        dataInput2 = []
        dataOutput2 = []
        for pair in self:
            val1 = strdataInput[pair[0].ECFPname_withsite]
            val2 = rpdataInput[pair.rp1.ECFPname_withsite]
            val = np.append(val1,val2)
            dataInput2.append(val)
            dataOutput2.append([pair.barrier1, pair.heat1])

            val1 = strdataInput[pair[1].ECFPname_withsite]
            val2 = rpdataInput[pair.rp2.ECFPname_withsite]
            val = np.append(val1,val2)
            dataInput2.append(val)
            dataOutput2.append([pair.barrier2, pair.heat2])

        input2 = np.array(dataInput2)
        result2 = np.array(dataOutput2)

        Dumppkl('RealIndextoInput_strwithsite.pkl',RealIndexToInput_str)
        Dumppkl('RealIndextoInput_rp.pkl',RealIndexToInput_rp)
        Dumppkl('datainput_info.pkl',input2)
        Dumppkl('dataout_info.pkl',result2)




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
        
        input = np.array(datain)
        result = np.array(dataresult)

        inputfile = open('stdoutput.pkl','wb')
        pickle.dump(dataOutput,inputfile)
        inputfile.close()


        inputfile=open('datain.pkl','wb')
        pickle.dump(input,inputfile)
        inputfile.close()
        
        outputfile=open('dataresult.pkl','wb')
        pickle.dump(result,outputfile)
        outputfile.close()

        Dumppkl('PattIndexDict.pkl',PattIndexDict)
        Dumppkl('RealIndextoInput_0.pkl',RealIndexToInput)

        return 



if __name__ == "__main__":
    test= AllPair()
    test.readfile('allgoodsect.arc',filemode =2)
    test.GetAllECFPname_surface()
    test.OutRPData('ECFP_withsite','testdata')
