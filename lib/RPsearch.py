from Arc2RPdata import RangeInfo
from simple_func_k import Loadpkl,Dumppkl
import numpy as np
import time
import glob
import os
from allstr_new import allstr
import datetime

def PrintTargetrpexample(rpname,rpdict):
    try:
        pairlist =rpdict[rpname]
    except:
        print ('key error: no rpname_withsite %s'%rpname)
        return
    _tmp = allstr()
    for pair in pairlist:
        pair[0].label = pair[0].ECFPname_surface
        pair.TSstr.label = pair.rp1.ECFPname_withsite
        pair[1].label = pair[1].ECFPname_surface
        _tmp.append(pair[0])
        _tmp.append(pair.TSstr)
        _tmp.append(pair[1])
    _tmp.printall('%s.arc'%rpname)


    

pairall =Loadpkl('allrpdict_withsite_pairexample.pkl')
rp0all =Loadpkl('allrpdict0.pkl')

print ('Data search service strat: %s'%(datetime.datetime.now()))


while not glob.glob('Endsearch.kpl'):
    if glob.glob('rpsearch'):
        print ('find search target')
        searchlist = []
        f =  open('rpsearch','r')
        for line in f:
            if ('RPwithsite' in line):
                searchlist.append(line.split()[-1])
            if ('RP0' in line):
                rp0 = line.split()[-1]
                try:
                    namelist = rp0all[rp0]
                except:
                    print ('key error: no rpname0 %s'%rp0)
                for rpname in namelist:
                    searchlist.append(rpname)
        f.close()

        for name in searchlist:
            PrintTargetrpexample(name,pairall)
        print ('search done %s'%(datetime.datetime.now()))
        os.system('rm -rf rpsearch')    

    if glob.glob('printall.kpl'):
        for namelist in rp0all.values():
            for rpname in namelist:
                PrintTargetrpexample(rpname,pairall)
        print ('print all done %s'%(datetime.datetime.now()))
        os.system('rm -rf printall.kpl')


    time.sleep(10)

print ('Data search service stop: %s'%(datetime.datetime.now()))
