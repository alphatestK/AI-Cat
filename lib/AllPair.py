import os
#from allstr_new import allstr
import numpy as np
#import matplotlib.pyplot as plt
#from reactpair import Pair as allpair
import  pickle
from ECFP import ECFP
from ECFP import AllStr as allstr

class Pair(list):
    def __init__(self,str1,str2,TS,TSstr):
        list.__init__(self)
        self.append(str1)
        self.append(str2)
        #self.append(TSstr)
        self.TS= TS
        self.sort(key= lambda X: X.name)
        self.name = '%s'%(self[0].name+self[1].name)
        self.TSstr =TSstr
        return


class AllPair(Pair):
    def __init__(self):
        list.__init__(self)

    def readfile(self,filename='allgoodsect.arc',filemode = 1):
        _tmp =allstr()

        _tmp.readfile(filename)
        #_tmp.GetAllsminame(colorflag = 1)
        _tmp.calAllFakeBondMaxtrix()
        _tmp.GetAllECFPname()
        for i in range(0,len(_tmp),3):
            if filemode == 1:
                #if (not _tmp[i].Lminstr) or (not _tmp[i+1].Lminstr): continue
                #if (not _tmp[i].sminame) or (not _tmp[i+1].sminame):  continue

                #_tmp[i].GenECFPname()
                #_tmp[i+1].GenECFPname()
                _tmp[i].name = _tmp[i].ECFPname
                _tmp[i+1].name = _tmp[i+1].ECFPname
                TS= _tmp[i+2].energy
                self.append(Pair(_tmp[i],_tmp[i+1], TS,_tmp[i+2]))
            elif filemode == 2:
                #if (not _tmp[i].Lminstr) or (not _tmp[i+2].Lminstr): continue
                #if (not _tmp[i].sminame) or (not _tmp[i+2].sminame):  continue
                #_tmp[i].GenECFPname()
                #_tmp[i+2].GenECFPname()
                _tmp[i].name = _tmp[i].ECFPname
                _tmp[i+2].name = _tmp[i+2].ECFPname

                TS= _tmp[i+1].energy
                self.append(Pair(_tmp[i],_tmp[i+2], TS,_tmp[i+1]))


        self.npair =len(self)

    def GenPairDict(self):
        self.allpair = {}
        for pair in self:
            if pair.name not in self.allpair.keys():
                self.allpair[pair.name] = []
                self.allpair[pair.name].append(pair)
            else:
                self.allpair[pair.name].append(pair)

        self.minpair={}
        for name in self.allpair.keys():
            self.minpair[name] =min(self.allpair[name],key=lambda x: x.TS)

