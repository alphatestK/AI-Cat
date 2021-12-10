import tensorflow as tf
from tensorflow import keras
import numpy as np
from Arc2RPdata import RangeInfo,FakeRP
import pickle
#from reactpair import Pair
from allreactstr import AllReactStr,Path
from allstr_new import allstr
from ECFP0 import AllStr
from simple_func_k import Loadpkl,Dumppkl
import os
import glob
import time
import datetime
from Arc2RPdata import AllPair
from metrics import  MaxError_barrier,MaxError_heat 


def Loadinfo():
    
    datain=open('RealIndextoInput_rp.pkl','rb')
    rpRealIndexToInput=pickle.load(datain)
    datain.close()
    
    datain=open('stdoutput.pkl','rb')
    stdout =pickle.load(datain)
    datain.close()
    
    dataout=open('PattIndexDict.pkl','rb')
    PattIndexDict=pickle.load(dataout)
    dataout.close()
    
    RealIndexToInput_0 =Loadpkl('RealIndextoInput_0.pkl')
    RealIndexToInput_withsite =Loadpkl('RealIndextoInput_strwithsite.pkl')
    allrpdict_withsite =Loadpkl('allrpdict_withsite.pkl')
    allrpdict0 =Loadpkl('allrpdict0.pkl')

    NNset = [rpRealIndexToInput, stdout, PattIndexDict, RealIndexToInput_0,\
             RealIndexToInput_withsite,allrpdict_withsite,allrpdict0]
    return NNset


def Loadinfo_Simple():

    dataout=open('PattIndexDict.pkl','rb')
    PattIndexDict=pickle.load(dataout)
    dataout.close()

    RealIndexToInput_0 =Loadpkl('RealIndextoInput_0.pkl')
    allrpdict0 =Loadpkl('allrpdict0.pkl')

    allrpdict0_rangeinfo = Loadpkl('allrpdict0_rangeinfo.pkl')
    allrpdict0_doublesampled = Loadpkl('allrpdict0_doublesampled.pkl')
    NNset = [PattIndexDict, RealIndexToInput_0, allrpdict0, allrpdict0_rangeinfo, allrpdict0_doublesampled]

    return NNset

def Loadinfo_mode5():
    dataout=open('PattIndexDict.pkl','rb')
    PattIndexDict=pickle.load(dataout)
    dataout.close()


    RealIndexToInput_0 =Loadpkl('RealIndextoInput_0.pkl')
    allrpdict0 =Loadpkl('allrpdict0.pkl')

    allrpdict0_rangeinfo = Loadpkl('allrpdict0_rangeinfo.pkl')
    allrpdict0_doublesampled = Loadpkl('allrpdict0_doublesampled.pkl')

    rpRealIndexToInput0= Loadpkl('RealIndextoInput_rp_0.pkl')
    strRealIndexToInput0= Loadpkl('RealIndextoInput_str_0.pkl')
    allrpdict0_ZPEinfo = Loadpkl('allrpdict0_ZPEinfo.pkl')
    NNset = [PattIndexDict, RealIndexToInput_0, allrpdict0, allrpdict0_rangeinfo, allrpdict0_doublesampled, rpRealIndexToInput0, strRealIndexToInput0, allrpdict0_ZPEinfo]
    return NNset

def Loadinfo_mode7():
    dataout=open('PattIndexDict.pkl','rb')
    PattIndexDict=pickle.load(dataout)
    dataout.close()


    RealIndexToInput_0 =Loadpkl('RealIndextoInput_0.pkl')
    allrpdict0 =Loadpkl('allrpdict0.pkl')

    allrpdict0_rangeinfo = Loadpkl('allrpdict0_rangeinfo.pkl')

    strRealIndexToInput0= Loadpkl('RealIndextoInput_str_0.pkl')
    NNset = [PattIndexDict, RealIndexToInput_0, allrpdict0, allrpdict0_rangeinfo, strRealIndexToInput0]
    return NNset



def LoadModel(mode):
    if mode == 1:
        model = keras.models.load_model('policymodel.h5')

        #model2= keras.models.load_model('net2.h5')
        model2= keras.models.load_model('infonet.h5', custom_objects={'MaxError_barrier':MaxError_barrier, 'MaxError_heat':MaxError_heat}, compile = False)
        #model2= keras.models.load_model('infonet.h5')
        return [mode,model, model2]

    elif mode == 2:
        model = keras.models.load_model('policymodel.h5')
    
        #model2= keras.models.load_model('net2.h5')
        barriermodel = keras.models.load_model('barriermodel.h5')
        heatmodel = keras.models.load_model('heatmodel.h5')
    
        return [mode,model,barriermodel, heatmodel]

    elif mode ==3:
        model = keras.models.load_model('policymodel.h5')

        #model2= keras.models.load_model('net2.h5')
        model2= keras.models.load_model('infonet.h5', custom_objects={'MaxError_barrier':MaxError_barrier, 'MaxError_heat':MaxError_heat}, compile = False)
        return [mode,model, model2]

    elif (mode%10 ==5) or (mode%10 == 7) :
        model = keras.models.load_model('policymodel.h5')
        model2= keras.models.load_model('infonet.h5', custom_objects={'MaxError_barrier':MaxError_barrier, 'MaxError_heat':MaxError_heat}, compile = False)
        return [mode,model, model2]


    elif mode == 0:
        model = keras.models.load_model('policymodel.h5')
        return [mode,model]

class Predictor(AllReactStr):
    def LoadInfo(self,NNset,modelset, ECFPdepth):
        if len(modelset) == 2:
            self.ECFPdepth = ECFPdepth
            self.PattIndexDict, self.RealIndexToInput0, self.allrpdict0, self.allrpdict0_rangeinfo, self.allrpdict0_doublesampled = NNset
            #self.PattIndexDict, self.RealIndexToInput0, self.allrpdict0, self.allrpdict0_rangeinfo, self.allrpdict0_rangeinfo_doublesampled = NNset
            #self.PattIndexDict, self.RealIndexToInput0, self.allrpdict0, self.purerp0dict = NNset
            #print self.allrpdict0_rangeinfo
            self.mode = modelset[0]
            self.model = modelset[1]
            return

        if ( modelset[0]%10 == 5):
            self.ECFPdepth = ECFPdepth
            self.PattIndexDict, self.RealIndexToInput0, self.allrpdict0, self.allrpdict0_rangeinfo, self.allrpdict0_doublesampled, self.rpRealIndexToInput_0, self.strRealIndexToInput_0,self.allrpdict0_ZPEinfo = NNset
            self.mode = modelset[0]
            self.model = modelset[1]
            self.rpmodel = modelset[2]
            return

        if ( modelset[0]%10 == 7):
            self.ECFPdepth = ECFPdepth
            self.PattIndexDict, self.RealIndexToInput0, self.allrpdict0, self.allrpdict0_rangeinfo, self.strRealIndexToInput_0 = NNset
            self.mode = modelset[0]
            self.model = modelset[1]
            self.rpmodel = modelset[2]
            return



        self.ECFPdepth = ECFPdepth
        self.rpRealIndexToInput, self.stdout, self.PattIndexDict,self.RealIndexToInput0,\
            self.RealIndexToInput_withsite,self.allrpdict_withsite,self.allrpdict0 =NNset

        self.mode= modelset[0]
        if ((self.mode)%2 == 1):
            self.model,self.rpmodel = modelset[1:]
            self.LoadNNInfo(self.model,self.RealIndexToInput0,self.PattIndexDict,self.ECFPdepth)
            self.LoadRPinfoNet(self.mode, self.rpmodel,self.rpRealIndexToInput, self.RealIndexToInput_withsite, self.allrpdict0, self.allrpdict_withsite)

        elif ((self.mode)%2 == 0):
            self.model,self.barriermodel, self.heatmodel = modelset[1:]
            self.LoadNNInfo(self.model,self.RealIndexToInput0,self.PattIndexDict,self.ECFPdepth)
            self.Loadtwomodel(self.mode, self.barriermodel,self.heatmodel,self.rpRealIndexToInput, self.RealIndexToInput_withsite,\
                                 self.allrpdict0, self.allrpdict_withsite)

        return 

    def SimplePredict(self,Str):
        if len(self) == 1:
            self.pop(0)
        self.append(Str)
        self.PredictProduct_Simple()
        return self.rpdict[Str.ECFPname]

#    def PredictOnemodel(self,Str):
#        if len(self) == 1:
#            self.pop(0)
#        self.append(Str)
#        self.PredictProduct_Singlemodel()
#        return self.rpdict[Str.ECFPname]

    def Predict(self,Str,selectsurface = 0):
        if len(self) == 1:
            self.pop(0)
        self.append(Str)

        if (selectsurface):
            self.PredictProduct_fixedsurface(selectsurface, LaddZPE= self.LaddZPE)
        else:
            if self.mode ==0:
                self.PredictProduct_Singlemodel()
            else:
                self.PredictProduct()
        return  self.rpdict[Str.ECFPname]

    def Predict_fixedsurface(self,Str,selectsurface):
        if len(self) == 1:
            self.pop(0)
        self.append(Str)
        self.PredictProduct_fixedsurface(selectsurface, LaddZPE= self.LaddZPE)
        return self.rpdict[Str.ECFPname]

    def Predict_selectrp(self,pair,rp,LGibbs, selectsurface = False):
        if len(self) == 1:
            self.pop(0)
        self.append(pair[0])
        if ((self.mode)%2 == 1):
            barrier, heat =self.RunRpinfoNet(rp,pair[0],self.rpmodel,pair[1], selectsurface = selectsurface)
        elif ((self.mode)%2 == 0):
            barrier, heat =self.RunRpinfoNet_twomodel(rp,pair[0],self.barriermodel, self.heatmodel)
        if LGibbs == 1:
            deltaTS, deltaTS_ts = self.AddGibbsEnergy(pair[0],pair[1])
            barrier = barrier + deltaTS_ts
            heat = heat + deltaTS
        return barrier, heat


class MCTS(object):
    def __init__(self, NNset,modelset, LaddZPE = 1):

#        self.model,self.rpmodel, self.rpRealIndexToInput, self.stdout, self.PattIndexDict,self.RealIndexToInput0,\
#                            self.RealIndexToInput_withsite,self.allrpdict_withsite,self.allrpdict0 =NNset
        self.predictor = Predictor()
        self.ECFPdepth = 3
        self.NNmode = modelset[0]
        self.predictor.LoadInfo(NNset,modelset,self.ECFPdepth)
        self.predictor.LaddZPE = LaddZPE
        self.rpdict = {}
        self.expandrpdict = {}
        self.expandrpdict_fixedsurface = {}
        self.allfragdict = {}
        self.terminaldict = {}


    def SetRootNode(self):
        self.allpath = []
        self.allstep = []
        self.allnode = {}
        self.rootnode.parent =0 
        self.rootnode.state = 'root'
        self.rootnode.nowpath = [self.rootnode.Str.ECFPname]
        self.rootnode.path = Path()
        self.rootnode.path.score = 99900000000000
        self.rootnode.bscore = 99900000000000

        self.rootnode.nodeid = 1
        self.allnode[1] = self.rootnode
        #self.rootnode.path.append(self.rootnode.Str)
        self.rootnode.path.namechain = [self.rootnode.Str.ECFPname]

    def SetTarget(self,rootstr,targetstr,selectsurface = 0, maxpathdepth =12):
        self.mode = 1
        rootstr.GenECFPname()
        targetstr.GenECFPname()
        self.rootnode = Node(rootstr)
        self.targetstr = targetstr
        self.SetRootNode()
        self.pathresult = []

        self.selectsurface = selectsurface
        self.maxdepth = maxpathdepth
        self.nowbestscore = 2500000000000

    def SetFuzzyTarget(self,rootstr, targetlist, selectsurface = 0, maxpathdepth =12):

        self.mode = 2
        rootstr.GenECFPname()
        #targetlist.GetAllECFPname()
        self.multitarget = []
        for Str in targetlist:
            Str.GenECFPname()
            self.multitarget.append(Str.ECFPname)
        print (self.multitarget)
        self.rootnode = Node(rootstr)
        #self.targetstr = targetstr
        self.SetRootNode()
        self.pathresult = []

        self.selectsurface = selectsurface
        self.maxdepth = maxpathdepth
        self.nowbestscore = 2500000000000

    def SetFuzzyTarget2(self,rootstr, targetlist, selectsurface = 0, maxpathdepth =12):

        self.mode = 3
        rootstr.GenECFPname()
        #targetlist.GetAllECFPname()
        self.multitarget = []
        for formula in targetlist:
            self.multitarget.append(formula)
        print (self.multitarget)
        self.rootnode = Node(rootstr)
        #self.targetstr = targetstr
        self.SetRootNode()
        self.pathresult = []

        self.selectsurface = selectsurface
        self.maxdepth = maxpathdepth
        self.nowbestscore = 2500000000000



    def Search(self,times):
        #self.rpdict = {}

        l = self.CheckNodeTerminal(self.rootnode.Str)
        if (l ==1 ):
            f= open('result','w')
            f.write('input str is already target\n')
            f.close()
            return

        for i in range(times):
            _ta = time.time()
            Lout =0
            self.allpath = []
            r = self.SelectRollout()
            self.Update(self.simulatenode,r)
            _tb = time.time()
            print ('cycle: %d' %(self.rootnode.N))
            print ('cycle time: %f'%(_tb-_ta))

            #if glob.glob('OutResult'):
            #    os.system('rm -f OutResult')
            #    Lout = 1
            
            if (i%100 ==0 or Lout ): self.OutResult()


            #if glob.glob('EndSearch'):
            #    self.DumpTree()
            #    break



    def Connect_Path(self, name1,ppath,depth, maxdepth = 10):
        if depth > maxdepth :
            return
        pmid=[]
        pmidpath =[]

        try :
            rppair = self.rpdict[name1.ECFPname]
        except:
            #_tmp = AllReactStr()
            #_tmp.append(name1)
            #_tmp.LoadNNInfo(self.model,self.RealIndexToInput0,self.PattIndexDict,self.ECFPdepth)
            ##_tmp.LoadRPinfoNet(self.rpmodel,self.rpRealIndexToInput, self.RealIndexToInput_withsite, self.allrpdict0, self.allrpdict_withsite)
            #_tmp.PredictProduct_Simple()
            #self.rpdict[name1.ECFPname] =  _tmp.rpdict[name1.ECFPname]
            if self.NNmode == 0:
                self.rpdict[name1.ECFPname] = self.predictor.Predict(name1,self.selectsurface )
            else:
                self.rpdict[name1.ECFPname] = self.predictor.SimplePredict(name1)
            rppair = self.rpdict[name1.ECFPname]


        for i,pairinfo in enumerate(rppair):
            pair = pairinfo[0]
            #pair.score = pairinfo[1]
           
            lend = self.CheckNodeTerminal(pair[1])
            #print (pair[1].ECFPname)
            #print ("lend: %s"%(lend))
            if lend ==1:
                _ppath = Path()
                _ppath.copy(ppath)
                _ppath.append(pair)
                _ppath.namechain.append(pair[1].ECFPname)
                #_ppath.score = ppath.score* pairinfo[1]
                #_ppath.score =min( ppath.score,pairinfo[1])
                self.allpath.append(_ppath)
            elif i < 2:
                if (pair[1].ECFPname not in ppath.namechain ):
                    pmid.append(pair[1])
                    _ppath = Path()
                    _ppath.copy(ppath)
                    _ppath.append(pair)
                    _ppath.namechain.append(pair[1].ECFPname)
                    #_ppath.score = ppath.score* pairinfo[1]
                    #_ppath.score =min( ppath.score,pairinfo[1])
                    pmidpath.append(_ppath)

        depth =depth +1
        if len(pmid) > 0:
            #pname2 = [name2]*len(pmid)
            pdepth = [depth]*len(pmid)
            pmaxdepth = [maxdepth]*len(pmid)
            #return map(self.Connect_Path, pmid,pname2,pmidpath,pdepth,pmaxdepth)
            return map(self.Connect_Path, pmid,pmidpath,pdepth,pmaxdepth)
        else :
            return

    def ExpandNodes(self,node,selectsurface):
        if selectsurface == 0:
            if (node.Str.ECFPname in self.expandrpdict.keys()):
                plist = self.expandrpdict[node.Str.ECFPname]
            else:
                self.expandrpdict[node.Str.ECFPname] = self.predictor.Predict(node.Str)
                plist = self.expandrpdict[node.Str.ECFPname]

        else:
            if selectsurface not in self.expandrpdict_fixedsurface.keys():
                self.expandrpdict_fixedsurface[selectsurface] = {}

            if node.Str.ECFPname in self.expandrpdict_fixedsurface[selectsurface].keys():
                plist = self.expandrpdict_fixedsurface[selectsurface][node.Str.ECFPname]
            else:
                self.expandrpdict_fixedsurface[selectsurface][node.Str.ECFPname]= self.predictor.Predict_fixedsurface(node.Str,selectsurface)
                plist = self.expandrpdict_fixedsurface[selectsurface][node.Str.ECFPname]

        return plist

    def CheckNodeTerminal(self, Str):
        if self.mode  == 1:
            if Str.ECFPname == self.targetstr.ECFPname:
                return 1
            else:
                return 0

        elif self.mode == 2:
            #Str= node.Str
            if Str.ECFPname in self.allfragdict.keys():
                _allfrag =self.allfragdict[Str.ECFPname]
            else:
                Str.bmx1D = []
                for line in Str.bmx2D:
                    Str.bmx1D.extend(line)
                Str.group=Str.Segmolecular()

                outstr = AllStr()
                fragments = [[] for i in np.unique(Str.group)]
                for id,atom in enumerate(Str.atom):
                    if (atom.ele> 18): continue
                    atom.id = id
                    fragments[Str.group[id]-1].append(id)
                for frag in fragments:
                    if (len(frag) == 0): continue
                    sub = Str.GenSubStr(frag)
                    sub.CalMass()
                    outstr.append(sub)

                
                substr = sorted(outstr, key=lambda x: x.mass, reverse=True)
                Str.AllFragECFPname = []
                for frag in substr:
                    frag.ECFPname = (frag.GenECFPname(savefp =0 ))[0]
                    Str.AllFragECFPname.append(frag.ECFPname)
                #Str.allfrag = substr
                #print ('Str.AllFragECFPname')
                #print Str.AllFragECFPname
                self.allfragdict[Str.ECFPname] = Str.AllFragECFPname
                _allfrag = self.allfragdict[Str.ECFPname]

            r = 0
            for target in self.multitarget:
                if target in _allfrag:
                    r = 1
                    break
            return r            

        elif self.mode == 3:
            if Str.ECFPname in self.allfragdict.keys():
               # _allfrag =self.allfragdict[Str.ECFPname]
                r= self.terminaldict[Str.ECFPname]
                return r
            else:
                Str.bmx1D = []
                for line in Str.bmx2D:
                    Str.bmx1D.extend(line)
                Str.group=Str.Segmolecular()

                Str.Get_AllFragFormula()
                self.allfragdict[Str.ECFPname] = Str.allfragchemformula
                _allfrag = self.allfragdict[Str.ECFPname]


            r = 0
            for target in self.multitarget:
                for frag in _allfrag:
                    if (target in frag):
                        r= 1
                        break
            self.terminaldict[Str.ECFPname] = r
            return r 

    def OutPathInfo(self):
        f= open('allstep','w')
        for line in self.allstep:
            f.write('%s\n'%(line))
        f.close()

        nodeinfo = []
        for nodeid in self.allnode.keys():
            nodeinfo.append([nodeid,self.allnode[nodeid].bscore])

        nodeinfo.sort(key = lambda X: X[1])

        print('number of node: %d'%(len(nodeinfo)))

        _pathresult = []
        allfinalstate = []

        for inode in range(0,len(nodeinfo),100):
            tmpall= AllStr()
            for info in nodeinfo[inode:inode+100]:
                tmpall.append(self.allnode[info[0]].Str)
            tmpall.GetAllsminame_fromreactstr_serial()

            for i,info in enumerate(nodeinfo[inode:inode+100]):
                if ( ('94m' not in tmpall[i].sminame) and  (tmpall[i].sminame not in allfinalstate)):
                    allfinalstate.append(tmpall[i].sminame)
                    _pathresult.append(self.allnode[info[0]].path)
            if (len(_pathresult) >= 50):
                break

        tmpall= AllStr()
        for i,path in enumerate(_pathresult[:50]):
            for pair in path:
                tmpall.append(pair[0])
                tmpall.append(pair[1])
        tmpall.GetAllsminame_fromreactstr_serial()
             

        f= open('goodresult','w')
        iStr =0
        for i,path in enumerate(_pathresult[:50]):
            if (len(path) == 0): continue
            f.write('-------path-%d-------\n'%(i+1))
            f.write('score-%f----barrier--%f--heat--%f\n'%(path.score,path.barrier,path.heat))
            for ipair,pair in enumerate(path):
                f.write('%d  RPwithsite %s   %s   %s   %s   %s   %s\n' %((ipair+1),tmpall[iStr+1].surfaceinfo,tmpall[iStr+1].rpname,tmpall[iStr].sminame,tmpall[iStr+1].sminame,pair[1].barrier,pair[1].heat))
                iStr =iStr+2
        f.close()

            


    def SelectRollout(self):
        _node =self.rootnode
        #ppath =Path()
        #ppath.namechain=[_node.Str.ECFPname]
        #ppath.score = 1
        depth =0


        while (_node.state != 'unvisited'):
            if not hasattr(_node, 'subnode'):
                _t01 = time.time()

                sublist = self.ExpandNodes(_node,self.selectsurface)
                _node.ExpandSubNodes(sublist)
                for tmpnode in _node.subnode:
                    self.allstep.append('%s -> %s  bscore: %f  barrier: %f  heat: %f' %(_node.nodeid,tmpnode.nodeid,\
                                                tmpnode.bscore, tmpnode.Str.barrier, tmpnode.Str.heat))
                    self.allnode[tmpnode.nodeid] = tmpnode
                #print ('expansion')
                #print (len(_node.subnode))
                _t02 = time.time()
                #print ('expansion time: %f'%(_t02-_t01))
            newnode = _node.SelectSubNode(self.nowbestscore)

            Lright =self.CheckNodeTerminal(newnode.Str)
            if (newnode.state == 'terminal') and (Lright == 0):
                #print ('fail: terminal for no legal action')
                self.simulatenode = newnode
                return -1

#            ppath.append([_node.Str,newnode.Str])
#            ppath.namechain.append(newnode.Str.ECFPname)
#            ppath.score = ppath.score*(newnode.P)
#            depth = depth + 1
            if Lright == 1:
                newnode.state = 'terminal'
                self.pathresult.append(newnode.path)
                if newnode.bscore < self.nowbestscore:
                    self.nowbestscore = newnode.bscore
                self.simulatenode = newnode
                return 1
            _node = newnode

        self.simulatenode = _node
        #print self.simulatenode.Str.ECFPname
        
        _t03 = time.time()
        if self.mode == 1:
            r = self.Simulate(_node.Str,self.targetstr,_node.path,depth,self.maxdepth)
        elif self.mode == 2:
            r = self.Simulate(_node.Str,self.multitarget,_node.path,depth,self.maxdepth)
        elif self.mode == 3:
            r = self.Simulate(_node.Str,self.multitarget,_node.path,depth,self.maxdepth)
        _t04 = time.time()
        #print ('rollout time: %f'%(_t04-_t03))

        return r


    def Simulate(self,name1,targetstr,path,depth,maxdepth):
        # not need targetsrt now
        #print ('depth: %d' %depth)
        self.Connect_Path(name1,path,depth,maxdepth)
        #self.Connect_Path(name1,targetstr,path,depth,maxdepth)
        if len(self.allpath) > 0:
            if (self.NNmode == 0) and len(self.pathresult) == 0:
                for path in self.allpath:
                    path.GetBarrier_list()
                self.allpath.sort(key = lambda X:X.score)
                if self.allpath[0].score < self.nowbestscore:
                    self.nowbestscore = self.allpath[0].score
                    self.pathresult.append(self.allpath[0])
                    print ('add path from rollout')
            r = 1 
        else:
            r = -1

        return r


    def Update(self,node, r):
        if node.state == 'unvisited':
            node.state = 'visited'
        _node = node
        _node.N = _node.N +1
        _node.Q = _node.Q +r


        while _node.parent != 0:
            _node = _node.parent 
            _node.N = _node.N +1 
            _node.Q = _node.Q +r


    def OutResult(self):
#        for path in self.allpath:
#            path.GetBarrier_list()
#
#        self.allpath.sort(key = lambda X:X.score)
#
#
#        for path in self.allpath[:15]:
#            self.pathresult.append(path)

        if len(self.pathresult) == 0:
            return

#        for path in self.pathresult:
#            path.GetBarrier_list()

        print ('allpath in:  %d'%(len(self.pathresult)))
        self.pathresult.sort(key = lambda X:X.score)

        barrierpairlist = []
        sameRScount = 0
        printcount = 0
        _pathresult = []
        for path in self.pathresult:
            if (printcount > 21):  break
            repeatKeyRS = 0
            for pair in path:
                pairname = '%s'%(pair[0].ECFPname)+'-%s'%(pair[1].ECFPname)
                if pairname in barrierpairlist:
                    repeatKeyRS= 1
                    break

            if (repeatKeyRS):
                if (sameRScount > 2):
                    continue
                else:
                    sameRScount += 1
                    printcount += 1
            else:
                pairname = '%s'%(path[path.barrierpair][0].ECFPname)+'-%s'%(path[path.barrierpair][1].ECFPname)
                #barrierpairlist.append(path[path.barrierpair].name)
                barrierpairlist.append(pairname)
                sameRScount = 1
                printcount += 1

            _pathresult.append(path)         

        #_pathresult = self.pathresult[:15]

        print ('allpath:  %d'%(len(self.pathresult)))
        del self.pathresult

        #inputfile = open('pathresult.pkl','wb')
        #pickle.dump(_pathresult,inputfile)
        #inputfile.close()

       
        tmpall= AllStr()
        for i,path in enumerate(_pathresult[:20]):
            for pair in path:
                tmpall.append(pair[0])
                tmpall.append(pair[1])
        tmpall.GetAllsminame_fromreactstr_serial()

        f= open('result','w')
        iStr =0
        for i,path in enumerate(_pathresult[:20]):
            f.write('-------path-%d-------\n'%(i+1))
            f.write('score-%f----barrier--%f--heat--%f\n'%(path.score,path.barrier,path.heat))
            for ipair,pair in enumerate(path):
                f.write('%d  RPwithsite %s   %s   %s   %s   %s   %s\n' %((ipair+1),tmpall[iStr+1].surfaceinfo,tmpall[iStr+1].rpname,tmpall[iStr].sminame,tmpall[iStr+1].sminame,pair[1].barrier,pair[1].heat))
                iStr =iStr+2
        f.close()
        self.pathresult= _pathresult
        self.nowbestscore = self.pathresult[0].score

        #print ('allpath:  %d'%(len(self.pathresult)))

    def DumpTree(self):

        """
        very slow for dump/load, some problem not solved.
        """


        Dumppkl('MCTS.pkl',self)

    def PathReproduce(self,filename, LGibbs= 0, selectsurface =False, LZPE = 1):
        f = open('result','w')
        test = AllPair()
        test.readfile(filename, filemode = 2,pairsortflag = False)
        f.write('path\n')
        for ipair,pair in enumerate(test):
            ECFPname0 = pair.rp1.ECFP0name
            #ECFPname_withsite = pair.rp1.ECFPname_withsite
            print ECFPname0
            #try:
            #    rp = self.predictor.allrpdict_withsite[ECFPname_withsite]
            #except:
            #    print ('ERROR: no rp %s in dataset'%(ECFPname_withsite))
            #    continue
            rp = pair.rp1
            barrier,heat = self.predictor.Predict_selectrp(pair,rp, LGibbs,selectsurface)

            if (LZPE == 1):
                ZPEinfo =  self.predictor.allrpdict0_ZPEinfo[ECFPname0]
                barrier = barrier + ZPEinfo[0]
                heat = heat + ZPEinfo[1]

            if (selectsurface):
                surface = selectsurface
            else:
                surface = 0
            f.write('%d  RP0  %s   %s   %s   %s   %s   %s\n' %((ipair+1),surface ,ECFPname0,pair[0].sminame,pair[1].sminame,barrier,heat))
        f.close()




class Node(object):
    def __init__(self,Str):
        self.Str = Str
        self.state = 'unvisited'
        self.N = 0
        self.Q = 0
        self.statename =  Str.ECFPname
        self.ECFPdepth = 3
        #self.NNset = NNset
        #self.model,self.rpmodel, self.rpRealIndexToInput, self.stdout, self.PattIndexDict,self.RealIndexToInput0,\
        #                    self.RealIndexToInput_withsite,self.allrpdict_withsite,self.allrpdict0 =NNset


#    def expand(self):
#        _tmp = AllReactStr()
#        _tmp.append(self.Str)
#        _tmp.LoadNNInfo(self.model,self.RealIndexToInput0,self.PattIndexDict,self.ECFPdepth)
#        _tmp.LoadRPinfoNet(self.rpmodel,self.rpRealIndexToInput, self.RealIndexToInput_withsite, self.allrpdict0, self.allrpdict_withsite)
#        _tmp.PredictProduct()
#        plist =  _tmp.rpdict[self.statename]
#
#        return plist
#
#
#    def expand_fixedsurface(self,selectsurface):
#        _tmp = AllReactStr()
#        _tmp.append(self.Str)
#        _tmp.LoadNNInfo(self.model,self.RealIndexToInput0,self.PattIndexDict,self.ECFPdepth)
#        _tmp.LoadRPinfoNet(self.rpmodel,self.rpRealIndexToInput, self.RealIndexToInput_withsite, self.allrpdict0, self.allrpdict_withsite)
#        _tmp.PredictProduct_fixedsurface(selectsurface)
#        plist =  _tmp.rpdict[self.statename]
#
#        return plist

    def ExpandSubNodes(self,sublist):
        self.subnode = []
#        if not selectsurface :
#            sublist = self.expand()
#        else:
#            sublist = self.expand_fixedsurface(selectsurface)

        k = 0
        for subpairinfo in sublist:
            if subpairinfo[0][1].ECFPname in self.path.namechain:
                continue
            #_sub = Node(subpairinfo[0][1],self.NNset)
            _sub = Node(subpairinfo[0][1])
            _sub.P = subpairinfo[1]
            #print (_sub.P)
            _sub.parent = self
            _sub.path = Path()
            _sub.path.copy(self.path)
            _sub.path.append([self.Str,_sub.Str])
            _sub.path.namechain.append(_sub.Str.ECFPname)
            _sub.path.GetBarrier_list()
            _sub.bscore = _sub.path.score
            k = k + 1
            _sub.nodeid = self.nodeid*10 + k
            self.subnode.append(_sub)
            #print (_sub.Str.rpname)

        sumK = 0
        for node in self.subnode: 
            #node.K = np.exp(-node.Str.barrier*96485/8.314/1000)
            node.K = np.exp(-node.Str.barrier*10)
            sumK = sumK + node.K

        for node in self.subnode: 
            node.P = node.K/sumK
            #print node.P

    def SelectSubNode(self,nowbestscore):
        if (len(self.subnode) == 0):
            #print ('no legal action:  terminal')
            self.state = 'terminal'
            selectnode = self
            return selectnode
        else:
            for node in self.subnode:
                if (node.state == 'terminal'):
                    node.value = -999999
                    continue
                if (node.bscore > (nowbestscore+300000000000)):
                    node.state = 'terminal'
                    node.value = -999999
                    continue
                #if (node.state == 'unvisited'):
                #    selectnode = node
                #    return selectnode
                node.value = float(node.Q)/float(node.N+1)+(np.sqrt(np.sqrt(node.P))*np.sqrt(self.N)/float(node.N+1))
  
            #print self.state
            #print self.subnode 
            #print len(self.subnode)
            #print self.N
            #print self.Q
            selectnode = max(self.subnode, key= lambda x: x.value)
            #print selectnode.value
            if selectnode.state == 'terminal':
                self.state = 'terminal'

        return selectnode

def ReadPara(file):
    f= open(file,'r')
    line=f.readline()
    para ={}
    while line:
        if len(line.split()) == 0:
            line=f.readline()
            continue
        elif line.split()[0] =='selectsurface':
            para['selectsurface']= int(line.split()[1])
        elif line.split()[0] =='reloadmodel':
            para['reloadmodel']= int(line.split()[1])
        elif line.split()[0] =='NNmode':
            para['NNmode']= int(line.split()[1])
        elif line.split()[0] =='mode':
            para['mode']= int(line.split()[1])
        elif line.split()[0] =='maxpathdepth':
            para['maxpathdepth']= int(line.split()[1])
        elif line.split()[0] =='inputfile':
            para['inputfile']= line.split()[1]
        elif line.split()[0] =='startfile':
            para['startfile']= line.split()[1]
        elif line.split()[0] =='targetfile':
            para['targetfile']= line.split()[1]
        elif line.split()[0] =='maxsearchtime':
            para['maxsearchtime']= int(line.split()[1])
        elif line.split()[0] =='LGibbs':
            para['LGibbs']= int(line.split()[1])
        elif line.split()[0] =='LaddZPE':
            para['LaddZPE']= int(line.split()[1])
#        elif line.split()[0] =='%block':
#            blockname = line.split()[-1]
#            lines= []
#            blockline=f.readline().split()
#            while blockline[0] != '%endblock':
#                lines.append(blockline)
#                blockline=f.readline().split()
#            if blockname == 'ReactPattern':
#                para[blockname]= PatternBlock(lines)
#            else:
#                para[blockname]=readblock(lines)
        line=f.readline()
    f.close()
    return para


if __name__ == "__main__":
   # inputfile = 'simple.arc'

    _t1 = time.time()
   
    para = ReadPara('console')
    NNmode = para['NNmode']
    if (NNmode ==0):
        NNset = Loadinfo_Simple()
    elif (NNmode%10 ==5) :
        NNset = Loadinfo_mode5()
    elif (NNmode%10 ==7) :
        NNset = Loadinfo_mode7()
    else:
        NNset = Loadinfo()
    modelset = LoadModel(NNmode) 
    
    #NNset = Loadinfo2()
    #test = MCTS(rootstr,targetstr, NNset, selectsurface = '111')
    MCtree = MCTS(NNset,modelset)
    _t2 = time.time()

    #MCtree.PathReproduce('lowestpath-CO-100', LGibbs = 1)

    print ('Load   time: %f'%(_t2-_t1))
    while not glob.glob('Endsearch.kpl'):
        if glob.glob('searchjob'):
            #try:
            print ('Start Search %s'%(datetime.datetime.now()))
            para = ReadPara('searchjob')
            if not glob.glob(para['inputfile']):
                print ('Error: not found inputfile %s' %(para['inputfile']))
                os.system('mv searchjob searchjob_done')
                continue
            if para['reloadmodel'] == 1:
                modelset = LoadModel(para['NNmode'])
                MCtree = MCTS(NNset, modelset, LaddZPE = para['LaddZPE'])


            if para['mode'] == 1:
                test= AllReactStr()
                test.readfile(para['inputfile'])
                test.GetAllsminame()

                rootstr = test[0]
                targetstr = test[1]
                MCtree.SetTarget(rootstr,targetstr, para['selectsurface'],para['maxpathdepth'] )

            elif para['mode'] == 2:
                rootstr =  AllReactStr()
                rootstr.readfile(para['startfile'])
                rootstr.GetAllsminame()
                targetlist =  AllReactStr()
                targetlist.readfile(para['targetfile'])
                targetlist.GetAllsminame()
                MCtree.SetFuzzyTarget(rootstr[0],targetlist, para['selectsurface'],para['maxpathdepth'])
            #while (not glob.glob('Endnowjob.kpl')) and (not (_searchtime > para['maxsearchtime'])):

            elif para['mode'] == 3:
                rootstr =  AllReactStr()
                rootstr.readfile(para['inputfile'])
                rootstr.GetAllsminame()
                targetlist = ['C2']
                MCtree.SetFuzzyTarget2(rootstr[0],targetlist, para['selectsurface'],para['maxpathdepth'])


            for i in range(0, para['maxsearchtime'],100):
                if (glob.glob('Endnowjob.kpl')):
                    os.system('rm -f Endnowjob.kpl')
                    break
                MCtree.Search(100)
                #_searchtime = _searchtime + 1000
            MCtree.OutPathInfo()
            os.system('mv searchjob searchjob_done')
            print ('End Search %s'%(datetime.datetime.now()))
            time.sleep(30)
            #except:
            #    os.system('mv searchjob searchjob_fail')
            #    print ('Fail Search %s'%(datetime.datetime.now()))

        if glob.glob('reproducejob'):
            para = ReadPara('reproducejob')
            if para['reloadmodel'] == 1:
                modelset = LoadModel(para['NNmode'])
                MCtree = MCTS(NNset, modelset)
            if not glob.glob(para['inputfile']):
                print ('Error: not found inputfile %s' %(para['inputfile']))
                os.system('mv reproducejob reproducejob_done')
                continue

            MCtree.PathReproduce(para['inputfile'],para['LGibbs'], para['selectsurface'], LZPE= para['LaddZPE'])
            os.system('mv reproducejob reproducejob_done')
            time.sleep(30)

        time.sleep(60)

        
#    _t3 = time.time()
#    test.Search(30000)
#    _t4 = time.time()
#    print ('Load   time: %f'%(_t2-_t1))
#    print ('Search time: %f'%(_t3-_t2))
#    print ('wall   time: %f'%(_t3-_t1))
#    print ('final   time: %f'%(_t4-_t1))
#
