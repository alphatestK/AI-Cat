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

    _tmp.printall('%s.arc'%rpname)


def PrintTargetrpexample2(rpname,workdir):
    os.system('cp %s/%s.arc .'%(workdir,rpname))

#def PrintTargetrpexample2(rpname,workdir,name1,name2):
#    os.system('cp %s/%s.arc .'%(workdir,rpname))


def GenReverseDict(rp0all):
    RPwithsite2RO0dict = {}
    for rp in rp0all.keys():
        rpwithsitelist = rp0all[rp]
        for item in rpwithsitelist:
            RPwithsite2RO0dict[item] = rp


    return RPwithsite2RO0dict
        
    

#pairall =Loadpkl('allrpdict_withsite_pairexample.pkl')
rp0all =Loadpkl('allrpdict0.pkl')
#RPwithsite2RO0dict = GenReverseDict(rp0all)
#Dumppkl('rpwithsiteDict.pkl',RPwithsite2RO0dict)

print ('Data search service strat: %s'%(datetime.datetime.now()))


while not glob.glob('Endsearch.kpl'):
    if glob.glob('rpsearch'):
        print ('find search target')
        searchlist = []
        f =  open('rpsearch','r')
        for line in f:
            if ('RPwithsite' in line):
                searchlist.append(line.split()[3])
            if ('RP0' in line):
                rp0 = line.split()[1]
                try:
                    namelist = rp0all[rp0]
                except:
                    print ('key error: no rpname0 %s'%rp0)
                for rpname in namelist:
                    searchlist.append(rpname)
        f.close()

        for name in searchlist:
            PrintTargetrpexample2(name,'Arcsave_all/')
            #PrintTargetrpexample(name,pairall)
        print ('search done %s'%(datetime.datetime.now()))
        os.system('rm -rf rpsearch')    

    if glob.glob('pathprint'):
        print ('path connection')
        searchlist = []
        f =  open('pathprint','r')
        for line in f:
            if ('RPwithsite' in line):
                searchlist.append(line.split()[3])
            if ('RP0' in line):
                #rp0 = line.split()[1]
                rp0 = line.split()[3]
                try:
                    namelist = rp0all[rp0]
                except:
                    print ('key error: no rpname0 %s'%rp0)
                for rpname in namelist:
                    searchlist.append(rpname)
        f.close()

        out = allstr()
        for name in searchlist:
            file = 'Arcsave_all/'+'%s'%(name)+'.arc'
            test = allstr()
            test.readfile(file)
            #lp =min(len(test),3)
            lp =3
            for i in range(0,lp,3):
                out.append(test[i])
                out.append(test[i+1])
                out.append(test[i+2])
        out.printall('path.arc')
        os.system('rm -rf pathprint')

    if glob.glob('reproduce'):
        searchlist = []
        searchlist2 = []
        f =  open('reproduce','r')
        for line in f:
            if('RPwithsite' in line):
                searchlist.append([line.split()[3], line.split()[4],line.split()[5]])
            if('RP0' in line):
                searchlist2.append([line.split()[3], line.split()[4],line.split()[5]])
        f.close()

        out = allstr()
        if (len(searchlist) > 0):
            for info in searchlist:
                name = info[0]
                file = 'Arcsave_all/'+'%s'%(name)+'.arc'
                test = allstr()
                test.readfile(file)
                test.GetAllsminame()
                for i in range(0,len(test),3):
                    if test[i].sminame == info[1] and test[i+2].sminame == info[2]:
                        out.append(test[i])
                        out.append(test[i+1])
                        out.append(test[i+2])

        if (len(searchlist2) > 0)
            for info in searchlist2:
                name = info[0]
                file = 'Arcsave_all_rp0/'+'%s'%(name)+'.arc'
                test = allstr()
                test.readfile(file)
                test.GetAllsminame()
                for i in range(0,len(test),3):
                    if test[i].sminame == info[1] and test[i+2].sminame == info[2]:
                        out.append(test[i])
                        out.append(test[i+1])
                        out.append(test[i+2])


        out.printall('reproduce.arc')
        os.system('rm -rf reproduce')


    if glob.glob('printall.kpl'):
        for namelist in rp0all.values():
            for rpname in namelist:
                PrintTargetrpexample(rpname,pairall)
        print ('print all done %s'%(datetime.datetime.now()))
        os.system('rm -rf printall.kpl')


    time.sleep(10)

print ('Data search service stop: %s'%(datetime.datetime.now()))
