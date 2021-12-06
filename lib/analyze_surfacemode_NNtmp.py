#! /usr/bin/env python

import os
#from allstr_new import allstr
import numpy as np
import matplotlib.pyplot as plt
#from reactpair import Pair as allpair
import  pickle
from ECFP_withsite import ECFP as ECFPwithsite
from ECFP_withsite import AllStr as allstr_ecfpwithsite

from ECFP0 import ECFP as ECFP0
from ECFP0 import AllStr as allstr_ecfp0
class Pair(list):
    def __init__(self,str1,str2,TS,TSstr,sortflag = True):
        list.__init__(self)
        self.append(str1)
        self.append(str2)
        #self.append(TSstr)
        self.TS= TS
        if sortflag: self.sort(key= lambda X: X.name)
        self.name = '%s'%(self[0].name+self[1].name)
        self.TSstr =TSstr
        return

class Path(list):
    def __init__(self):
        list.__init__(self)
        self.maxTS =999

    def copy(self,cpath):
        for pair in cpath:
            self.append(pair)
        self.maxTS = cpath.maxTS
        self.namechain = cpath.namechain[:]

    def printpath(self,fname):
        _path= allstr()
        for pair in self:
            _path.append(pair[0])
            _path.append(pair[1])
        _path.printall(fname)

    def printpath_TS(self,fname):
        _path= allstr()
        for pair in self:
            _path.append(pair[0])
            _path.append(pair.TSstr)
            _path.append(pair[1])
        _path.printall(fname)

    def GetBarrier(self):
        Ecurrent = 0
        bmax = 0
        for pair in self:
            if pair[0].energy < Ecurrent:
                Ecurrent = pair[0].energy
            bcurrent = pair.TS- Ecurrent
            if bmax < bcurrent:
                bmax = bcurrent
            if pair[1].energy < Ecurrent:
                Ecurrent = pair[1].energy
        self.barrier= bmax
        self.score = bmax*1000+ len(self)
        return
            

    def GetBarrier_new(self,allmin):
        Ecurrent = 0
        bmax = 0
        barrierpair =0 
        for ipair,pair in enumerate(self):
            if allmin[pair[0].name].energy < Ecurrent:
                Ecurrent = allmin[pair[0].name].energy
            bcurrent = pair.TS- Ecurrent
            if bmax < bcurrent:
                bmax = bcurrent
                barrierpair = ipair
            if allmin[pair[1].name].energy < Ecurrent:
                Ecurrent = allmin[pair[1].name].energy
        self.barrier= bmax
        self.score = bmax*1000+ len(self)
        self.barrierpair = barrierpair
        return

    def GetBarrier_list(self,allmin):
        Ecurrent = 0
        bmax = 0

        self.GetBarrier_new(allmin)
        barrierpair = self.barrierpair
        self.barrierlist= [self.barrier]


        while (barrierpair+1) != len(self):
            _tmppath = Path()
            for ipair in range(barrierpair+1,len(self)):
                _tmppath.append(self[ipair])
            _tmppath.GetBarrier_new(allmin)
            barrierpair = barrierpair+1+_tmppath.barrierpair
            self.barrierlist.append(_tmppath.barrier)

        #print self.barrierlist
        self.score = 0
        if len(self.barrierlist)> 4:
            print('maybe too long to cal score')
        for i in range(min(4,len(self.barrierlist))):
            self.score= self.score+self.barrierlist[i]*(1000**(4-i))
        self.score = self.score +len(self)

        return

    def EndtoendSort(self, name1):
        for i,pair in enumerate(self): 
            pair.sort(key= lambda X: X.name)
            if pair[0].name == name1:
                name1 = pair[1].name
            else :
                #print 'reverse', i
                #print pair[0].name,pair[1].name
                pair.sort(key= lambda X: X.name, reverse = True)
                #if pair[0].name != name1:
                #    pair.sort(key= lambda X: X.name)
                #print pair[0].name,pair[1].name
                name1 = pair[1].name
                self.pop(i)
                self.insert(i,pair)

    def plot(self,fname):
        IS,TS,FS =[],[],[]
        line,linex =[],[]
        for i,pair in enumerate(self):
            IS.append(pair[0].energy)
            TS.append(pair.TSstr.energy)
            FS.append(pair[1].energy)
            line.append([pair[0].energy,pair.TSstr.energy,pair[1].energy])
            linex.append([i+1,i+1.5,i+2])
        plotISx = np.linspace(1, len(self), len(self))
        plotTSx = np.linspace(1.5, len(self) + 0.5, len(self))
        plotFSx = np.linspace(2, len(self) + 1, len(self))
        fig = plt.figure(figsize = (8, 6), dpi = 400)
        plt.scatter(plotISx, IS, s = 400, color = 'black', marker = '_', lw = 5)
        plt.scatter(plotTSx, TS, s = 400, color = 'red', marker = '_', lw = 5)
        plt.scatter(plotFSx, FS, s = 400, color = 'black', marker = '_', lw = 5)
        for pairlinex,pairline in zip(linex,line):
            plt.plot(pairlinex, pairline,'b-')

        #plt.xticks([])
        #plt.yticks()
        plt.xlabel('Reaction Coordinate', fontsize = 24, fontweight = 'bold')
        plt.ylabel('$\Delta$E [eV per atom]', fontsize = 24, fontweight = 'bold')
        plt.savefig(fname, dpi = 400)

        return

    def plot_new(self,fname,allmin):
        IS,TS,FS, ISmin,FSmin =[],[],[],[],[]
        line,linex =[],[]
        for i,pair in enumerate(self):
            IS.append(pair[0].energy)
            ISmin.append(allmin[pair[0].name].energy)
            TS.append(pair.TSstr.energy)
            FS.append(pair[1].energy)
            FSmin.append(allmin[pair[1].name].energy)
            line.append([pair[0].energy,pair.TSstr.energy,pair[1].energy])
            linex.append([i+1,i+1.5,i+2])
        plotISx = np.linspace(1, len(self), len(self))
        plotTSx = np.linspace(1.5, len(self) + 0.5, len(self))
        plotFSx = np.linspace(2, len(self) + 1, len(self))
        fig = plt.figure(figsize = (8, 6), dpi = 400)
        plt.scatter(plotISx, IS, s = 400, color = 'black', marker = '_', lw = 5)
        plt.scatter(plotTSx, TS, s = 400, color = 'red', marker = '_', lw = 5)
        plt.scatter(plotFSx, FS, s = 400, color = 'black', marker = '_', lw = 5)
        plt.scatter(plotISx, ISmin, s = 400, color = 'green', marker = '_', lw = 5)
        plt.scatter(plotFSx, FSmin, s = 400, color = 'green', marker = '_', lw = 5)
        for pairlinex,pairline in zip(linex,line):
            plt.plot(pairlinex, pairline,'b-')

        #plt.xticks([])
        #plt.yticks()
        plt.xlabel('Reaction Coordinate', fontsize = 24, fontweight = 'bold')
        plt.ylabel('$\Delta$E [eV per atom]', fontsize = 24, fontweight = 'bold')
        plt.savefig(fname, dpi = 400)

        return


    def plot_paper3(self,allmin):
        point = []
        zeroe = allmin[self[0][0].name].energy
        print (zeroe)
        point = [0,0]
        for i,pair in enumerate(self):
            #if pair.TSstr.energy-zeroe > 3000 : zeroe = allmin[pair[0].name].energy-point[-1]
            point.append(pair.TSstr.energy-zeroe)
            point.append(pair.TSstr.energy-zeroe)
            point.append(allmin[pair[1].name].energy-zeroe)
            point.append(allmin[pair[1].name].energy-zeroe)
            print (allmin[pair[1].name].energy)
            print (pair[1].name)
            print (pair.TSstr.energy - allmin[pair[0].name].energy)
            print (allmin[pair[1].name].energy- allmin[pair[0].name].energy)

        plotx = np.linspace(1-0.125, len(self)+1+0.125, 4*len(self)+2)
        return plotx,point


class AllPair(Pair):
    def __init__(self):
        list.__init__(self)

    def readfile(self,filename='allgoodsect.arc',filemode = 1, Lallmininput = 0,allminfile= 'allname.arc', pairsortflag = True):
        _tmp =allstr()
        
        _tmp.readfile(filename)
        _tmp.GetAllsminame(colorflag = 1)
        _tmp.GetAllECFPname()
        for i in range(0,len(_tmp),3):
            if filemode == 1:
                #if (not _tmp[i].Lminstr) or (not _tmp[i+1].Lminstr): continue
                if (not _tmp[i].sminame) or (not _tmp[i+1].sminame):  continue

                #_tmp[i].GenECFPname()
                #_tmp[i+1].GenECFPname()
                _tmp[i].name = _tmp[i].ECFPname
                _tmp[i+1].name = _tmp[i+1].ECFPname
                TS= _tmp[i+2].energy
                self.append(Pair(_tmp[i],_tmp[i+1], TS,_tmp[i+2], pairsortflag))
            elif filemode == 2:
                #if (not _tmp[i].Lminstr) or (not _tmp[i+2].Lminstr): continue
                if (not _tmp[i].sminame) or (not _tmp[i+2].sminame):  continue
                #_tmp[i].GenECFPname()
                #_tmp[i+2].GenECFPname()
                _tmp[i].name = _tmp[i].ECFPname
                _tmp[i+2].name = _tmp[i+2].ECFPname

                TS= _tmp[i+1].energy
                self.append(Pair(_tmp[i],_tmp[i+2], TS,_tmp[i+1], pairsortflag))


        self.npair =len(self)

        self.Lallmin = Lallmininput

        if Lallmininput == 1:
            _allmin = allstr()
            _allmin.readfile(allminfile)
            _allmin.GetAllsminame(colorflag = 1)
            _allmin.GetAllECFPname()

            self.allname = {}
            self.allnameid = {}
            self.allstrbyname= {}

            for Str in _allmin:
                Str.name = Str.ECFPname
                if Str.name not in self.allname.keys():
                    nameid = len(self.allname)
                    self.allname[Str.name]= nameid
                    self.allnameid[nameid]= Str.name
                    self.allstrbyname[Str.name]= [Str]
                else:
                    self.allstrbyname[Str.name].append(Str)
            self.allstrmin = {}
            for name in self.allname.keys():
                self.allstrmin[name]= min(self.allstrbyname[name],key=lambda x: x.energy)

        if Lallmininput == 2:
            allminfile= 'addmin.arc'
            self.addallmin = allminfile


    def ReadPath(self, filename, filemode =2 ,Lallmininput = 0 ,allminfile= 'allname.arc'):
        self.readfile(filename,filemode = filemode,Lallmininput= Lallmininput, allminfile=allminfile, pairsortflag=False)
        outpath = Path()
        for pair in self:
            outpath.append(pair)
        return outpath


    def AllName(self):
        self.allname = {}
        self.allnameid = {}
        self.allstrbyname= {}
        for pair in self:
            for Str in pair:
                if Str.name not in self.allname.keys():
                    nameid = len(self.allname)
                    self.allname[Str.name]= nameid
                    self.allnameid[nameid]= Str.name
                    self.allstrbyname[Str.name]= [Str]
                else:
                    self.allstrbyname[Str.name].append(Str)

        if self.Lallmin == 2:
            _allmin = allstr()
            _allmin.readfile(self.addallmin)
            _allmin.GetAllsminame(colorflag = 1)
            _allmin.GetAllECFPname()
            for Str in _allmin:
                Str.name = Str.ECFPname
                if Str.name not in self.allname.keys():
                    nameid = len(self.allname)
                    self.allname[Str.name]= nameid
                    self.allnameid[nameid]= Str.name
                    self.allstrbyname[Str.name]= [Str]
                else:
                    self.allstrbyname[Str.name].append(Str)

        self.allstrmin = {}
        for name in self.allname.keys():
            self.allstrmin[name]= min(self.allstrbyname[name],key=lambda x: x.energy)
       
        self.Lallmin =1


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

    def PrintPair(self,name1,name2,  printTSstr =True):
        if self.Lallmin != 1: self.AllName()
        self.GenPairDict()
        namelist = [name1,name2]
        namelist.sort()
        targetname = '%s'%(namelist[0]+namelist[1])

        self.allpair[targetname].sort(key = lambda x: x.TS)
        _tmp = allstr()
        for pair in self.allpair[targetname]:
            print (pair.TS)
            _tmp.append(pair[0])
            if printTSstr: _tmp.append(pair.TSstr)
            _tmp.append(pair[1])
        _tmp.printall('selectpair.arc')


    def PrintPair_list(self,filename, printTSstr =True):
        if self.Lallmin != 1: self.AllName()

        self.GenPairDict()

        _tmp =allstr()
        _tmp.readfile(filename)
        _tmp.GetAllsminame(colorflag =1)
        _tmp.GetAllECFPname()

        out= allstr()
        for i in xrange(0,len(_tmp),2):
            name1 = _tmp[i].ECFPname
            name2 = _tmp[i+1].ECFPname

            namelist = [name1,name2]
            namelist.sort()
            targetname = '%s'%(namelist[0]+namelist[1])
   
            try: 
                self.allpair[targetname].sort(key = lambda x: x.TS)
                for pair in self.allpair[targetname]:
                    print (pair.TS)
                    out.append(pair[0])
                    if printTSstr: out.append(pair.TSstr)
                    out.append(pair[1])
            except:
                print ('no pair %d'%(i/2))


        out.printall('selectpair_list.arc')


    def GenPairdict_byname(self):
        if self.Lallmin != 1: self.AllName()

        self.GenPairDict()
        self.rpdict={}
        for name in self.allname.keys():
            self.rpdict[name]=[]
            for pair in self.minpair.values():
                if pair[0].name ==name or pair[1].name ==name:
                    self.rpdict[name].append(pair)
            self.rpdict[name].sort(key = lambda X: X.TS)
        #print self.rpdict

    def OutPairdict_byname(self,name, printTSstr =True):
        print (name)
        _tmp = allstr()
        for pair in self.rpdict[name]:
            if pair[0].name == name:
                print (self.allname[pair[1].name],pair[1].name,pair.TS)
                _tmp.append(pair[0])
                if printTSstr: _tmp.append(pair.TSstr)
                _tmp.append(pair[1])
            else:
                print (self.allname[pair[0].name],pair[0].name,pair.TS)
                _tmp.append(pair[1])
                if printTSstr: _tmp.append(pair.TSstr)
                _tmp.append(pair[0])
        _tmp.printall('spepair.arc')

    def Outallpair(self,pts= 1):
        _tmp= allstr()
        for pair in self.minpair.values():
            _tmp.append(pair[0])
            if pts:
                _tmp.append(pair.TSstr)
            _tmp.append(pair[1])
        _tmp.printall('allpair.arc')

    def connect_path_all(self, name1,ppath,depth, maxdepth =10):
        if depth > 10 or len(name1) == 0:
            return

        pmid=[]
        pmidpath =[]
        for i,pair in enumerate(self.rpdict[name1]):
            #if pair[0].name == name2 or pair[1].name== name2:
            if i < 10:
                #if pair[0].name == name1 and (pair[1].name not in ppath.namechain ) and (pair.TS - pair[0].energy) < 2.5:
                if (name1 in pair[0].name) and (pair[1].name not in ppath.namechain) and (pair.TS - pair[0].energy) < 2.5 and (pair.TS -self.zeroenergy) < 3:
                    pmid.append(pair[1].name) 
                    #print (pair.TS -self.zeroenergy)
                    _ppath = Path()
                    _ppath.copy(ppath)
                    _ppath.append(pair)
                    _ppath.namechain.append(pair[1].name)
                    pmidpath.append(_ppath)
                    self.searchtrace.append(pair[1].name)
                    name2 = pair[1].name
                    if name2 not in self.AllPath_fromspeStr.keys():
                        self.AllPath_fromspeStr[name2]= [_ppath]
                    else:
                        self.AllPath_fromspeStr[name2].append(_ppath)


                #elif pair[1].name == name1 and (pair[0].name not in ppath.namechain ) and (pair.TS - pair[1].energy) < 2.5:
                elif (name1 in pair[1].name) and (pair[0].name not in ppath.namechain) and (pair.TS - pair[1].energy) < 2.5 and (pair.TS -self.zeroenergy) <3:
                    pmid.append(pair[0].name)
                    #print (pair.TS -self.zeroenergy)
                    _ppath = Path()
                    _ppath.copy(ppath)
                    _ppath.append(pair)
                    _ppath.namechain.append(pair[0].name)
                    pmidpath.append(_ppath)
                    self.searchtrace.append(pair[0].name)
                    name2 = pair[0].name
                    if name2 not in self.AllPath_fromspeStr.keys():
                        self.AllPath_fromspeStr[name2]= [_ppath]
                    else:
                        self.AllPath_fromspeStr[name2].append(_ppath)


        #print 'cycle',depth
        #for name in pmid:
        #    print name
        depth =depth +1
        if len(pmid) > 0:
            pdepth = [depth]*len(pmid)
            return map(self.connect_path_all, pmid,pmidpath,pdepth)
        else:
            return


    def connect_path(self, name1,name2,ppath,depth,maxdepth = 10):
        if depth > 10 or len(name1) ==0:
            return
#        if len(name1)!= 1:
#            r_name1=[name1[0]]
#            name1.pop(0)
#            return self.connect_path(r_name1,name2,ppath,depth) ,self.connect_path(name1,name2,ppath,depth)
#        name1 =name1[0]
        pmid=[]
        pmidpath =[]
        for i,pair in enumerate(self.rpdict[name1]):
            #if pair[0].name == name2 or pair[1].name== name2:
            if  name2 in pair[0].name or name2 in pair[1].name:
                _ppath = Path()
                _ppath.copy(ppath)
                _ppath.append(pair)
                _ppath.namechain.append(pair[1].name)
                self.allpath.append(_ppath) 
            elif i < 10:
                #if pair[0].name == name1 and (pair[1].name not in ppath.namechain ) and (pair.TS - pair[0].energy) < 2.5:
                if (name1 in pair[0].name) and (pair[1].name not in ppath.namechain ) and (pair.TS - pair[0].energy) < 2.5 and (pair.TS -self.zeroenergy) < 3 :
                    #print (pair.TS -self.zeroenergy)
                    pmid.append(pair[1].name) 
                    _ppath = Path()
                    _ppath.copy(ppath)
                    _ppath.append(pair)
                    _ppath.namechain.append(pair[1].name)
                    pmidpath.append(_ppath)
                    self.searchtrace.append(pair[1].name)
                    
                #elif pair[1].name == name1 and (pair[0].name not in ppath.namechain ) and (pair.TS - pair[1].energy) < 2.5:
                elif (name1 in pair[1].name) and (pair[0].name not in ppath.namechain) and (pair.TS - pair[1].energy) < 2.5 and (pair.TS -self.zeroenergy) <3:
                    #print (pair.TS -self.zeroenergy)
                    pmid.append(pair[0].name)
                    _ppath = Path()
                    _ppath.copy(ppath)
                    _ppath.append(pair)
                    _ppath.namechain.append(pair[0].name)
                    pmidpath.append(_ppath)
                    self.searchtrace.append(pair[0].name)
        #print 'cycle',depth
        #for name in pmid:
        #    print name
        depth =depth +1
        if len(pmid) > 0:
            pname2 = [name2]*len(pmid)
            pdepth = [depth]*len(pmid)
            return map(self.connect_path, pmid,pname2,pmidpath,pdepth)
        else:
            return

    def GetTarget(self,file):
        _tmp =allstr()
        _tmp.readfile(file)
        _tmp.GetAllsminame()
        #_tmp.GenCanoicalSmiles()
        _tmp[0].GenECFPname()
        _tmp[1].GenECFPname()

        name1 = _tmp[0].ECFPname
        name2 = _tmp[1].ECFPname
        

        return name1,name2

    def GetSingleTarget(self,file):
        _tmp =allstr()
        _tmp.readfile(file)
        _tmp.GetAllsminame()
        _tmp[0].GenECFPname()
        #_tmp.GenCanoicalSmiles()
        name1 = _tmp[0].ECFPname
        return name1


    def OutInfo(self):
        f= open('allname','w')
        f.write('allname appear in pair\n')
        f.write('  Id   Energy      Npair    Name \n')
        _tmp=allstr()
        for i in range(len(self.allname)):
            f.write('%4d  %12.6f  %4d %-50s\n'%(i,self.allstrmin[self.allnameid[i]].energy,\
                                            len(self.rpdict[self.allnameid[i]]), self.allnameid[i] ))
            _tmp.append(self.allstrmin[self.allnameid[i]])
        _tmp.printall('allname.arc')

        f.close()
        
    
    def OutInfo_script(self):
        self.OutInfo()
        os.system('cat allname')
        #print ('  Id  Energy   Name      Npairlink')
        #for i in range(len(self.allname)):
        #    print ('%4d    %50s     %d'%(i,self.allnameid[i], len(self.rpdict[self.allnameid[i]])))

    def GetAllPathfromSpestr(self,name1):
        barrierdict ={}
        self.lowpathall = []
        for name in self.allname.keys():
            if name != name1:
                self.findpath(name1,name, output = False)
                if len(self.allpath) == 0: continue
                barrierdict[name]= self.allpath[0][-1][1]
                self.lowpathall.append([name,self.allpath[0].barrier,self.allpath[0].score])
        self.lowpathall.sort(key = lambda X : X[2])
        return self.lowpathall

    def OutAllPathfromSpestr(self,name1):
        #self.GetAllPathfromSpestr(name1)
        self.findpath_all(name1)
        _tmp = allstr()
        for item in self.lowpathall:
            name = item[0]
            print (item[0],item[1])
            _tmp.append(self.allstrmin[name])
        _tmp.printall('alllow.arc')


    def findpath_all(self,name1,output =True ,printstring = False):
        self.zeroenergy = self.allstrmin[name1].energy
        self.AllPath_fromspeStr= {}
        self.searchtrace =[name1]
        ppath =Path()
        ppath.namechain=[name1]
        depth =0

        barrierdict ={}
        self.lowpathall = []
        self.connect_path_all(name1,ppath,depth)
        for name in self.AllPath_fromspeStr.keys():
            allpath = self.AllPath_fromspeStr[name]
            for path in allpath:
                path.EndtoendSort(name1)
                path.GetBarrier_list(self.allstrmin)
            allpath.sort(key= lambda x: x.score)
            barrierdict[name]= allpath[0][-1][1]
            self.lowpathall.append([name,allpath[0].barrier,allpath[0].score])

        self.lowpathall.sort(key = lambda X : X[2])
        return self.lowpathall


    def findpath(self,name1,name2,output =True ,printstring = False):
        self.zeroenergy = self.allstrmin[name1].energy
        self.allpath =[]
        self.searchtrace =[name1]
        ppath =Path()
        ppath.namechain=[name1]
        depth =0
        self.connect_path(name1,name2,ppath,depth)
        #print self.allpath

        
        #print 'start sort path'
        for path in self.allpath:
            path.EndtoendSort(name1)
            #path.GetBarrier()
            path.GetBarrier_list(self.allstrmin)
            #path.maxTS= (max(path,key=lambda x:x.TS)).TS
            #print path.maxTS
        self.allpath.sort(key= lambda x: x.score)

        if not output:
            return

        inputfile = open('allpath.pkl','wb')
        pickle.dump(self.allpath,inputfile)
        inputfile.close()


        #print 'end sort path'
        os.system('rm -rf lowestpath-*')
        f= open('lowestpath','w')
        f.write('lowestpath %s  %s\n'%(name1,name2))
        for i,path in enumerate(self.allpath[:10]):
            path.EndtoendSort(name1)
            path.printpath_TS('lowestpath-%d'%(i+1))
            if printstring : path.plot('lowestpath-%d.png'%(i+1))
            f.write('------lowestpath-%d-barrier-%f-score-%f------\n'%(i+1,path.barrier,path.score))
            for pair in path:
                f.write('IS/TS/FS  %14.8f  %14.8f  %14.8f   %s   %s\n' %(pair[0].energy,pair.TS,pair[1].energy,pair[0].sminame,pair[1].sminame))
        f.close()
        #print r[:5]
        #    _path
        #    for pair in path:
        #        pair[0] = 1
        #        pair[1] 


if __name__ == "__main__":
    test= AllPair()
    test.readfile()
    test.GenPairdict_byname()
    
    name1,name2 =test.GetTarget('target.arc')
    #name1= test[0][1].name
    #name2= test[-21][0].name
    
    test.OutInfo()
    test.OutPairdict_byname(test.allnameid[18])
    name1 =test.allnameid[32]
    name2 =test.allnameid[40]
   # print name1 
   # print name2
    #test.findpath(name1,name2)
