import numpy as np
from simple_func_k import hashlist


class ReactPattern(object):
    def __init__(self,fplist,bmxIS,bmxFS,FSfplist):
        self.fplist =fplist
        self.bmxIS = bmxIS
        self.bmxFS = bmxFS
        self.FSfplist = FSfplist
        self.GenName()

    def GenName(self):
        tmp =[]
        tmp= self.fplist[:]
        tmp.sort()

        tmp2 =[]
        tmp2 =self.FSfplist[:]
        tmp2.sort()

        tmp.extend(tmp2[:])

        string = ''
        for index in tmp:
            string = string +str(index)
        self.name = string



class ReactPatternNew(object):
    def __init__(self):
        self.contents =[]
        self.fpmatch= []
        self.atomall = []
        self.atommatch = []
    
    def GenMatch(self,atomid0,ECFP,matchlink):
        atomid = atomid0
        for index in matchlink:
            fp = ECFP.genbylayer[1][atomid]
            atomid = fp.trace[index]
        return atomid

#    def GenName(self,ISECFP,FSECFP,atomlist0):
#        tmpIS = []
#        tmpFS = []
#        for iatom in atomlist0:
#            tmpIS.append(ISECFP.genbylayer[1][iatom].index)
#            tmpFS.append(FSECFP.genbylayer[1][iatom].index)
#
#        tmpIS.sort()
#        tmpFS.sort()
#        tmpIS.extend(tmpFS[:])
#
#
#        self.ECFP= tmpIS
#
#        string = ''
#        for index in tmpIS:
#            string = string +str(index)
#        self.name = string
#
#        return

    def GenName(self,ISECFP,FSECFP,atomlist0):
        tmpIS = []
        tmpFS = []
        for iatom in atomlist0:
            tmpIS.append(hashlist(['IS',ISECFP.genbylayer[1][iatom].index]))
            tmpFS.append(hashlist(['FS',FSECFP.genbylayer[1][iatom].index]))

        tmpIS.sort()
        tmpFS.sort()
        tmpIS.extend(tmpFS[:])


        self.ECFP= tmpIS

        string = ''
        for index in tmpIS:
            string = string +str(index)
        self.name = string

        return



    def BmxSort(self,ISbmx,FSbmx):
        rpbmxIS = np.zeros((len(self.atommatch),len(self.atommatch)),dtype =np.int)
        rpbmxFS = np.zeros((len(self.atommatch),len(self.atommatch)),dtype =np.int)

        for i,atominfoi in enumerate(self.atommatch):
            atomi = atominfoi[1]
            for j,atominfoj in enumerate(self.atommatch):
                atomj = atominfoj[1]
                rpbmxIS[i][j] = ISbmx[atomi][atomj]
                rpbmxFS[i][j] = FSbmx[atomi][atomj]

        self.bmxIS = rpbmxIS
        self.bmxFS = rpbmxFS
        return


    def CheckSect(self,ECFP,sect,fp):
        atomid0 = fp.coreid
        lmatch = 1
        for joint in sect:
            atomid =self.GenMatch(atomid0,ECFP,joint[0])
            if ECFP.genbylayer[1][atomid].index != joint[1]:
                lmatch = 0
                break
        return lmatch

    def GetMatchFp(self,ECFP):
        self.sect =[]
        #print self.fpmatch

        for joint in self.fpmatch:
            if joint[0][0] == 0 and len(joint[0]) == 1:
                try:
                    self.sect.append(newsect)
                    newsect=[]
                except:
                    newsect=[]
            newsect.append(joint)
        self.sect.append(newsect)
       
        #print 'self.sect'
        #print self.sect
        matchfp = {}
        for i,sect in enumerate(self.sect):
            targetfp = sect[0][1]
            #print targetfp
            for fp in ECFP.allfp:
                #print fp.index
                if fp.index == targetfp:
                    lmatch = self.CheckSect(ECFP,sect,fp)
                    if lmatch == 1:
                        if i not in matchfp.keys():
                            matchfp[i] =[fp]
                        else:
                            matchfp[i].append(fp)
        #print matchfp

        if len(matchfp) != len(self.sect):
            return -1

        self.matchfp = matchfp
        self.matchcombine = []
        #print 'matchfp',matchfp[1]
        for fp in matchfp[0]:
            self.Dict2AllCombine(0,[],fp)

        #print 'self.matchcombine'
        #print self.matchcombine
        return self.matchcombine
    
    def Dict2AllCombine(self,idepth,matchlink,current):
        newlink = matchlink[:]
        newlink.append(current)
        if idepth == (len(self.matchfp)-1):
            self.matchcombine.append(newlink)
            return
        idepth = idepth +1

        fplist= []
        for fp in self.matchfp[idepth]:
            if fp not in newlink:
                fplist.append(fp)

        if len(fplist) > 0:
            newlinkl = [newlink]*len(fplist)
            idepthl = [idepth]*len(fplist)
            return map(self.Dict2AllCombine,idepthl,newlinkl,fplist)
        else:
            return

