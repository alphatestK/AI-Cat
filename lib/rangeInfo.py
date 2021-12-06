from Arc2RPdata import RangeInfo
from simple_func_k import Loadpkl,Dumppkl
import numpy as np
import time
import glob
import os
from allstr_new import allstr
import datetime
def PrintTargetrpexample2(rpname,workdir):
    os.system('cp %s/%s.arc .'%(workdir,rpname))

class RP(object):
    def __init__(self,barrier, heat):
        self.barrier = barrier
        self.heat = heat

def ReadDict(puredict):
    allrpdict_withsite_pure = {}
    f = open(puredict,'r')
    for line in f:
        rpname = (line.split()[0]).split('-')[2]
        barrier = float(line.split()[2])
        heat = float(line.split()[4])
        rp = RP(barrier,heat)
        if rpname not in allrpdict_withsite_pure.keys():
            allrpdict_withsite_pure[rpname] = RangeInfo(rp)
        else:
            allrpdict_withsite_pure[rpname].update(rp)
    return allrpdict_withsite_pure


allrpdict_withsite = Loadpkl('allrpdict_withsite.pkl')

maxrange = 0

rangelist = []
savekey = 0
for key in allrpdict_withsite:
    barrierrange = (allrpdict_withsite[key].rangeinfo.barriermax - allrpdict_withsite[key].rangeinfo.barriermin)
    rangelist.append([key, barrierrange])
    if barrierrange > maxrange:
        savekey = key
        maxrange = barrierrange 

print savekey
print maxrange

rangelist.sort(key = lambda x: x[1], reverse= True)

for item in rangelist[:10]:
    key = item[0]
    print item[0],item[1]
    PrintTargetrpexample2(key, 'Arcsave_all/')

allrpdict_withsite_pure = ReadDict('purereactdict')
for key in allrpdict_withsite_pure.keys():
    if allrpdict_withsite[key].rangeinfo.barriermin - allrpdict_withsite_pure[key].barriermin  < -0.5:
        print key


