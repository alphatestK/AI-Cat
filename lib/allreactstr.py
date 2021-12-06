
from ECFP0 import ECFP as ECFP0
from ECFP0 import AllStr
from ECFP_withsite import ECFP as ECFP_withsite
from reactpair_NN import Pair
from analyze_surfacemode_NNtmp import Path as pathstd
import numpy as np
import  pickle
from allstr_new import allstr

class Path(pathstd):
    def GetBarrier(self, Ecurrentin = 0 ):
        Ecurrent = Ecurrentin
        bmax = 0
        Emin = 0
        babs = 0
        barrierpair =0
        Ehere = Ecurrent
        for ipair,pair in enumerate(self):
            bcurrent = pair[1].barrier + Ecurrent
            babs = bcurrent - Emin
            Ecurrent = Ecurrent+pair[1].heat
            if bmax < babs:
                bmax = babs
                barrierpair = ipair
                Ehere = Ecurrent
            if Emin > Ecurrent:
                Emin = Ecurrent
        #print bmax
        self.barrier= bmax
        self.barrierpair = barrierpair
        self.Epair = Ehere
        self.score = bmax*1000+ len(self)
        self.Efinal = Ecurrent
        return



    def GetBarrier_list(self):

        self.GetBarrier()
        barrierpair = self.barrierpair
        Ein = self.Epair
        self.barrierlist= [self.barrier]
        self.heat = self.Efinal

        while (barrierpair+1) != len(self):
            _tmppath = Path()
            for ipair in range(barrierpair+1,len(self)):
                _tmppath.append(self[ipair])
            _tmppath.GetBarrier(Ein)
            barrierpair = barrierpair+1+_tmppath.barrierpair
            Ein = _tmppath.Epair
            self.barrierlist.append(_tmppath.barrier)

        #print self.barrierlist
        self.score = 0
        if len(self.barrierlist)> 4:
            print('maybe too long to cal score')
        for i in range(min(4,len(self.barrierlist))):
            self.score= self.score+np.round(self.barrierlist[i],2)*(1000**(4-i))
        self.score = self.score +len(self)

        return


    def plot_paper3(self):
        point = []
        Ecurrent = 0
        point = [0,0]
        for i,pair in enumerate(self):
            #if pair.TSstr.energy-zeroe > 3000 : zeroe = allmin[pair[0].name].energy-point[-1]
            bcurrent = pair[1].barrier + Ecurrent
            point.append(bcurrent)
            point.append(bcurrent)
            Ecurrent = Ecurrent+pair[1].heat
            point.append(Ecurrent)
            point.append(Ecurrent)


        plotx = np.linspace(1-0.125, len(self)+1+0.125, 4*len(self)+2)
        return plotx,point


class AllReactStr(Pair):
    def GenNNInput0(self,RealIndexToInput0,depth =3):
        testin= []
        for str in self:
            fp =ECFP0(str,depth)
            fp.GenECFP()
            str.ECFP0= fp
            inputvector,lwarn = self.Fp2InputVector(fp,RealIndexToInput0)
            testin.append(inputvector)
        self.datain =np.array(testin)

        return

    def GenECFP0(self,rp,depth = 0):
        input = []
        for Str in self:
            fp = rp.ECFP0
            inputvector2 = np.zeros(len(self.rpRealIndexToInput_0),dtype= np.int)
            for index in fp:
                try:
                    inputvector2[self.rpRealIndexToInput_0[index]] = inputvector2[self.rpRealIndexToInput_0[index]] + 1
                except:
                    print ('warning: contains unsupported index: %d'%index)

            fp =ECFP0(Str,depth)
            fp.GenECFP()
            Str.ECFP0= fp
            inputvector1,lwarn = self.Fp2InputVector(fp,self.strRealIndexToInput_0)

            inputvector =np.append(inputvector1,inputvector2)
            if (self.mode == 5) or (self.mode == 25):

                rangeinfo = self.allrpdict0_rangeinfo[rp.ECFP0name]
                addvector = [rangeinfo.barriermin,rangeinfo.barriermax,rangeinfo.barrieravg,\
                             rangeinfo.heatmax, rangeinfo.heatmin, rangeinfo.heatavg]
                #addvector = [rp.barriermax,rp.rangeinfo.barriermin,rp.rangeinfo.barrieravg,\
                #             rp.rangeinfo.heatmax, rp.rangeinfo.heatmin, rp.rangeinfo.heatavg]
                inputvector = np.append(inputvector,addvector)

            input.append(inputvector)

        datain = np.array(input)
        return datain

    def GenECFP0_twostr(self,rp,Str1,Str2,depth = 0):
        input = []
        fp =ECFP0(Str1,depth)
        fp.GenECFP()
        Str1.ECFP0= fp
        inputvector1,lwarn = self.Fp2InputVector(fp,self.strRealIndexToInput_0)

        fp =ECFP0(Str2,depth)
        fp.GenECFP()
        Str2.ECFP0= fp
        inputvector2,lwarn = self.Fp2InputVector(fp,self.strRealIndexToInput_0)

        inputvector =np.append(inputvector1,inputvector2)
        if (self.mode == 7) or (self.mode == 27):

            rangeinfo = self.allrpdict0_rangeinfo[rp.ECFP0name]
            addvector = [rangeinfo.barriermin,rangeinfo.barriermax,rangeinfo.barrieravg,\
                         rangeinfo.heatmax, rangeinfo.heatmin, rangeinfo.heatavg]
            #addvector = [rp.barriermax,rp.rangeinfo.barriermin,rp.rangeinfo.barrieravg,\
            #             rp.rangeinfo.heatmax, rp.rangeinfo.heatmin, rp.rangeinfo.heatavg]
            inputvector = np.append(inputvector,addvector)

        input.append(inputvector)

        datain = np.array(input)
        return datain


    def GenECFPwithsite(self,rp,RealIndexToInput_withsite,depth =3):
        input = []
        for Str in self:
            strlist =self.AddSurfacesite(Str,rp,depth,CalName = False)
            
            fp = rp.ECFP_withsite
            inputvector2 = np.zeros(len(self.rpRealIndexToInput),dtype= np.int)
            for index in fp:
                try:
                    #inputvector[RealIndexToInput[index]] =  1
                    inputvector2[self.rpRealIndexToInput[index]] = inputvector2[self.rpRealIndexToInput[index]] + 1
                except:
                    print ('warning: contains unsupported index: %d'%index)
    
            addstr = []
            for str_withsite in strlist:
                fp = ECFP_withsite(str_withsite,depth)
                fp.GenECFP()
    
                string = ''
                for index in fp.allindex:
                    #print index
                    string = string + str(len(fp.IndextoFp[index]))
                    string = string+ str(index)
                ECFPname= string

                if ECFPname not in addstr:
                    addstr.append(ECFPname)
                    inputvector1,lwarn= self.Fp2InputVector(fp, self.RealIndexToInput_withsite)
                    inputvector =np.append(inputvector1,inputvector2)
                    if self.mode == 3 or self.mode == 4:
                        addvector = [rp.rangeinfo.barriermin,rp.rangeinfo.barriermax,rp.rangeinfo.barrieravg,\
                                     rp.rangeinfo.heatmax, rp.rangeinfo.heatmin, rp.rangeinfo.heatavg]
                        inputvector = np.append(inputvector,addvector)

                    input.append(inputvector)
            if len(addstr)> 1:
                print ('multi match in add site info')
        datain = np.array(input)    
        return datain

    def RunPredictNN(self, model,PattIndexDict):
        out =model.predict(self.datain, batch_size=1)
        for i,vector in enumerate(out):
            predictpatt =[]
            #test = []
            for ikey,e in enumerate(vector):
                if e > 1e-10:
                    #print ikey,e
                    predictpatt.append([PattIndexDict[ikey],e])
                    #test.append([ikey,e])


            predictpatt.sort(key= lambda X : X[1],reverse= True)
            #test.sort(key= lambda X : X[1],reverse= True)
            #print test[:5]
            self[i].predictpatt= predictpatt
        return

    def Fp2InputVector(self,fp,RealIndexToInput):
        inputvector = np.zeros(len(RealIndexToInput),dtype= np.int)

        lwarn = 0
        for index in fp.allindex:
            try:
                #inputvector[RealIndexToInput[index]] =  1
                inputvector[RealIndexToInput[index]] =  len(fp.IndextoFp[index])
            except:
                print ('warning: contains unsupported index: %d'%index)
                lwarn = 1
        return inputvector,lwarn

    def RunRpinfoNet(self,rp,str,model ,FSstr =False, selectsurface= False):
        if ((self.mode)%10 == 5) :
            datain =self.GenECFP0(rp,depth =3)
        elif ((self.mode)%10 == 7) :
            datain =self.GenECFP0_twostr(rp,str,FSstr,depth =3)
        else:
            datain =self.GenECFPwithsite(rp,self.RealIndexToInput_withsite,depth =3)
        if len(datain) == 0:
            return 10,10

        if (self.mode >= 20):
            tmp = []
            if (not selectsurface):
                for array in datain:
                    v1 =np.append(array, np.array([1,0,0]))
                    tmp.append(v1)
                    v2 =np.append(array, np.array([0,1,0]))
                    tmp.append(v2)
                    v3 =np.append(array, np.array([0,0,1]))
                    tmp.append(v3)
            elif (selectsurface == 111):
                for array in datain:
                    v1 =np.append(array, np.array([1,0,0]))
                    tmp.append(v1)
            elif (selectsurface == 100):
                for array in datain:
                    v1 =np.append(array, np.array([0,1,0]))
                    tmp.append(v1)
            elif (selectsurface == 211):
                for array in datain:
                    v1 =np.append(array, np.array([0,0,1]))
                    tmp.append(v1)
            datain = np.array(tmp)


        #if (self.mode == 3) or (self.mode == 4):

        out =model.predict(datain, batch_size=1)
        #print rp.ECFPname_withsite, out
        outnew = list(out)
        outnew.sort(key=lambda x:x[0])
            
        return outnew[0][0],outnew[0][1]
        
    def RunRpinfoNet_twomodel(self,rp,Str,barriermodel, heatmodel):      
        datain =self.GenECFPwithsite(rp,self.RealIndexToInput_withsite,depth =3)
        if len(datain) == 0:
            return 10,10

        barrierout = barriermodel.predict(datain, batch_size=1)
        heatout =   heatmodel.predict(datain, batch_size=1)
        #print rp.ECFPname_withsite, out


        print barrierout
        _b = [x[0] for x in barrierout]
        _h = [x[0] for x in heatout]

        print _b
        #barrierout = list(barrierout)
        #print barrierout
        #heatout = list(heatout)

        outnew = zip(_b, _h)
        print outnew
        outnew.sort(key=lambda x:x[0])

        return outnew[0][0],outnew[0][1]
 

    def GenProduct_Simple(self):
        for iStr,Str in enumerate(self):
            self[iStr].product = []
            #print Str.predictpatt
            for i,pattinfo in enumerate(Str.predictpatt[:5]):
                #print 'try patt %d'%i
                rp = pattinfo[0]
                if rp.ECFP0name == '626280626280640619640619':
                    continue
                infolist = []
                psl=self.OutProduct(Str,rp,depth =3)
                for p in psl:
                    self[iStr].product.append([p,pattinfo[1]])
        return


    def GenProduct_Singlemodel(self, fixedsurface = 0):
        for iStr,Str in enumerate(self):
            self[iStr].product = []
            #print Str.predictpatt
            for i,pattinfo in enumerate(Str.predictpatt[:15]):
                #print 'try patt %d'%i
                rp = pattinfo[0]
                if rp.ECFP0name == '626280626280640619640619':
                    continue
                infolist = []
                psl=self.OutProduct(Str,rp,depth =3)
                for p in psl:
                    if (fixedsurface):
                        try:
                            add = 0
                            for rp in self.allrpdict0_doublesampled[rp.ECFP0name][:]:
                                if (int(rp.surfaceinfo.split('-')[0]) == fixedsurface):
                                    p.barrier  = rp.barrier
                                    p.heat   = rp.heat
                                    p.rpname   = rp.ECFPname_withsite
                                    p.surfaceinfo   = rp.surfaceinfo
                                    p.Lcheck = 1
                                    add = 1
                                    break
                            if (add == 0):
                                p.barrier  = self.allrpdict0_rangeinfo[rp.ECFP0name].barriermin
                                p.heat   = self.allrpdict0_rangeinfo[rp.ECFP0name].heatmin
                                p.rpname   = self.allrpdict0[rp.ECFP0name][0]
                                p.surfaceinfo   = 0
                                p.Lcheck = 0
                        except:
                            p.barrier  = self.allrpdict0_rangeinfo[rp.ECFP0name].barriermin
                            p.heat   = self.allrpdict0_rangeinfo[rp.ECFP0name].heatmin
                            p.rpname   = self.allrpdict0[rp.ECFP0name][0]
                            p.surfaceinfo   = 0
                            p.Lcheck = 0

                    
                    else:
                        try:
                            p.barrier  = self.allrpdict0_doublesampled[rp.ECFP0name][0].barrier
                            p.heat   = self.allrpdict0_doublesampled[rp.ECFP0name][0].heat
                            p.rpname   = self.allrpdict0_doublesampled[rp.ECFP0name][0].ECFPname_withsite
                            p.surfaceinfo   = self.allrpdict0_doublesampled[rp.ECFP0name][0].surfaceinfo
                            p.Lcheck = 1
                        except:
                            p.barrier  = self.allrpdict0_rangeinfo[rp.ECFP0name].barriermin
                            p.heat   = self.allrpdict0_rangeinfo[rp.ECFP0name].heatmin
                            p.rpname   = self.allrpdict0[rp.ECFP0name][0]
                            p.surfaceinfo   = 0
                            p.Lcheck = 0


                    #try:
                    #    p.barrier = self.purerp0dict[rp.ECFP0name][2]
                    #    p.heat  = self.purerp0dict[rp.ECFP0name][3]
                    #    p.rpname = self.purerp0dict[rp.ECFP0name][0]
                    #    p.surfaceinfo = self.purerp0dict[rp.ECFP0name][1]
                    #except:
                    #    print ('empty rp in dict: %s'%(rp.ECFP0name))
                    #    p.barrier = 0
                    #    p.heat = 0
                    #    p.rpname = 'wrong'
                    #    p.surfaceinfo = 0

                    self[iStr].product.append([p,pattinfo[1]])
        return

    def GenProduct_mode5(self,selectsurface= 0, LaddZPE = 1):
        for iStr,Str in enumerate(self):
            self[iStr].product = []
            for i,pattinfo in enumerate(Str.predictpatt[:15]):
                _rp = pattinfo[0]
                infolist = []
                if _rp.ECFP0name == '626280626280640619640619':
                    continue

                psl=self.OutProduct(Str,_rp,depth =3)
                if len(psl) == 0:
                    continue
                rangeinfo = self.allrpdict0_rangeinfo[_rp.ECFP0name] 
                print(i,_rp.ECFP0name)

                if (self.mode >= 5):
                    barrier,heat = self.RunRpinfoNet(_rp,Str,self.rpmodel, FSstr = psl[0], selectsurface = selectsurface)
                if (barrier < 0) or (barrier < heat):
                    print ('error: wrong barrier/heat')
                    continue
                if (barrier > rangeinfo.barriermax +0.15) or (barrier < rangeinfo.barriermin -0.15):
                    print ('error: barrier out of range')
                    print ('predict barrier: %f , barriermax: %f, barriermin: %f, rp0name: %s' \
                            %(barrier, rangeinfo.barriermax,rangeinfo.barriermin,_rp.ECFP0name))
                    continue
                if (heat > rangeinfo.heatmax +0.15) or (heat < rangeinfo.heatmin -0.15):
                    print ('error: heat out of range')
                    print ('predict heat: %f , heatmax: %f, heatmin: %f, rp0name: %s' %(heat, rangeinfo.heatmax,rangeinfo.heatmin,_rp.ECFP0name))
                    continue
                if (LaddZPE == 1):
                    ZPEinfo = self.allrpdict0_ZPEinfo[_rp.ECFP0name]
                    barrier = barrier + ZPEinfo[0]
                    heat = heat + ZPEinfo[1]


                if (not selectsurface):
                    infolist.append([barrier,heat,_rp.surfaceinfo,_rp.ECFP0name])
                else:
                    infolist.append([barrier,heat,selectsurface,_rp.ECFP0name])
#                try:
#                    add = 0
#                    for rp in self.allrpdict0_doublesampled[_rp.ECFP0name][:15]:
#                        if (self.mode == 5):
#                            barrier,heat = self.RunRpinfoNet(rp,Str,self.rpmodel)
#                        add = 1
#                        if (barrier < 0) or (barrier < heat):
#                            print ('error: wrong barrier/heat')
#                            continue
#                        if (barrier > rangeinfo.barriermax +0.15) or (barrier < rangeinfo.barriermin -0.15):
#                            print ('error: barrier out of range')
#                            continue
#                        if (heat > rangeinfo.heatmax +0.15) or (heat < rangeinfo.heatmin -0.15):
#                            print ('error: heat out of range')
#                            print ('predict heat: %f , heatmax: %f, heatmin: %f, rp0name: %s' %(heat, rangeinfo.heatmax,rangeinfo.heatmin, rp.ECFP0name))
#                            continue
#                        infolist.append([barrier,heat,rp.surfaceinfo,rp.ECFP0name])
#                    if (add == 0):
#                        if (self.mode == 5):
#                            barrier,heat = self.RunRpinfoNet(_rp,Str,self.rpmodel)
#                        if (barrier < 0) or (barrier < heat):
#                            print ('error: wrong barrier/heat')
#                            continue
#                        if (barrier > rangeinfo.barriermax +0.15) or (barrier < rangeinfo.barriermin -0.15):
#                            print ('error: barrier out of range')
#                            continue
#                        if (heat > rangeinfo.heatmax +0.15) or (heat < rangeinfo.heatmin -0.15):
#                            print ('error: heat out of range')
#                            print ('predict heat: %f , heatmax: %f, heatmin: %f, rp0name: %s' %(heat, rangeinfo.heatmax,rangeinfo.heatmin,_rp.ECFP0name))
#                            continue
#                        infolist.append([barrier,heat,_rp.surfaceinfo,_rp.ECFP0name])
#                except:        
#                    if (self.mode == 5):
#                        barrier,heat = self.RunRpinfoNet(_rp,Str,self.rpmodel)
#                    if (barrier < 0) or (barrier < heat):
#                        print ('error: wrong barrier/heat')
#                        continue
#                    if (barrier > rangeinfo.barriermax +0.15) or (barrier < rangeinfo.barriermin -0.15):
#                        print ('error: barrier out of range')
#                        continue
#                    if (heat > rangeinfo.heatmax +0.15) or (heat < rangeinfo.heatmin -0.15):
#                        print ('error: heat out of range')
#                        print ('predict heat: %f , heatmax: %f, heatmin: %f, rp0name: %s' %(heat, rangeinfo.heatmax,rangeinfo.heatmin,_rp.ECFP0name))
#                        continue
#                    infolist.append([barrier,heat,_rp.surfaceinfo,_rp.ECFP0name])


                if len(infolist) != 0:
                    infolist.sort(key = lambda x:x[0])
                    barrier = infolist[0][0]
                    heat    = infolist[0][1]
                    surfaceinfo = infolist[0][2]
                    rpname0 =  infolist[0][3]
                else:
                    continue

                for p in psl:
                    deltaTS= 0
                    deltaTS_ts = 0
                    #deltaTS, deltaTS_ts = self.AddGibbsEnergy(self[iStr],p)

                    p.barrier = barrier + deltaTS_ts
                    p.heat = heat + deltaTS

                    p.surfaceinfo = surfaceinfo
                    p.rpname = rpname0
                    self[iStr].product.append([p,pattinfo[1]])




    def GenProduct(self):
        for iStr,Str in enumerate(self):
            self[iStr].product = []
            #print Str.predictpatt
            for i,pattinfo in enumerate(Str.predictpatt[:15]):
                #print 'try patt %d'%i
                rp = pattinfo[0]
                infolist = []
                if rp.ECFP0name == '626280626280640619640619':
                    continue

                psl=self.OutProduct(Str,rp,depth =3)
                if len(psl) == 0:
                    #print 'wrong rp: fail to match rc'
                    continue
                #print (rp.ECFP0name,len(self.allrpdict0[rp.ECFP0name]))

                #_tmp = allstr()
                for ECFPname_withsite in self.allrpdict0[rp.ECFP0name][:15]:
                    #self.allrpdict_withsite_pairexample[ECFPname_withsite][0][0].label = ECFPname_withsite
                    #_tmp.append(self.allrpdict_withsite_pairexample[ECFPname_withsite][0][0])
                    #_tmp.append(self.allrpdict_withsite_pairexample[ECFPname_withsite][0].TSstr)
                    #_tmp.append(self.allrpdict_withsite_pairexample[ECFPname_withsite][0][1])
                    rp = self.allrpdict_withsite[ECFPname_withsite] 
                    if (self.mode)%2 == 1:
                        barrier,heat = self.RunRpinfoNet(rp,Str,self.rpmodel)
                    if (self.mode)%2 == 0:
                        barrier,heat = self.RunRpinfoNet_twomodel(rp,Str,self.barriermodel,self.heatmodel)
                        

                    #print rp.ECFPname_withsite,rp.rangeinfo.barriermax,rp.rangeinfo.barriermin, rp.rangeinfo.heatmax,rp.rangeinfo.heatmin
                    if (barrier < 0) or (barrier < heat):
                        print ('error: wrong barrier/heat')
                        continue
                    if (barrier > rp.rangeinfo.barriermax +0.15) or (barrier < rp.rangeinfo.barriermin -0.15):
                        print ('error: barrier out of range')
                        continue
                    if (heat > rp.rangeinfo.heatmax +0.15) or (heat < rp.rangeinfo.heatmin -0.15):
                        print ('error: heat out of range')
                        continue
                    infolist.append([barrier,heat,rp.surfaceinfo,rp.ECFPname_withsite])

                #_tmp.printall('%s.arc'%(rp.ECFP0name))
                if len(infolist) != 0:
                    infolist.sort(key = lambda x:x[0])
                    barrier = infolist[0][0]
                    heat    = infolist[0][1]
                    surfaceinfo = infolist[0][2]
                    rpname_withsite =  infolist[0][3]
                else:
                    continue

                for p in psl:
                    deltaTS, deltaTS_ts = self.AddGibbsEnergy(self[iStr],p)
                    #if deltaTS >= 0:
                    #    p.barrier = barrier + deltaTS
                    #    p.heat = heat + deltaTS
                    #if deltaTS < 0:
                    #    p.heat = heat+ deltaTS

                    
                    p.barrier = barrier + deltaTS_ts
                    p.heat = heat + deltaTS

                    #p.barrier = barrier
                    #p.heat = heat

                    p.surfaceinfo = surfaceinfo
                    p.rpname = rpname_withsite
                    #p.barrierinfo = patt.barrierinfo
                    #p.heatinfo = patt.heatinfo
                    #p.barrier = np.mean(p.barrierinfo)
                    #p.heat = np.mean(p.heatinfo)
                    self[iStr].product.append([p,pattinfo[1]])

    def AddGibbsEnergy(self, Str1, Str2):
        
        keyTSdict = {}
        TS_CO2 = 298.15*213.795*0.001*0.01036
        #TS_CO  = 298.15*197.653*0.001*0.01036
        TS_H2  = 298.15*130.680*0.001*0.01036 - 0.26
        TS_H2O = 298.15*188.834*0.001*0.01036 - 0.04

        keyTSdict['H2'] = TS_H2
        keyTSdict['H2O'] = TS_H2O
        keyTSdict['CO2'] = TS_CO2


        keyfpdict = {}
        keyfragdict = {}

        keyfpdict['H2'] = 9193356
        keyfpdict['CO2'] = 7093542
        keyfpdict['H2O'] = 5207851   

        fp1 = Str1.ECFP0 
        fp2 =ECFP0(Str2,3)
        fp2.GenECFP()
        #print fp1.allindex
        #print fp2.allindex
        #print fp1.IndextoFp
        #print fp2.IndextoFp


        for key,fpindex in keyfpdict.items():
            try:
                nkey1 = len(fp1.IndextoFp[fpindex])
                #print len(fp1.IndextoFp[fpindex])
            except:
                #print ("not found %s in fp1" %(fpindex))
                nkey1 = 0

            try:
                nkey2 = len(fp2.IndextoFp[fpindex])
                #print len(fp2.IndextoFp[fpindex])
            except:
                #print ("not found %s in fp2" %(fpindex))
                nkey2 = 0

            keyfragdict[key] = nkey1 -nkey2

        keyfragdict['H2'] = keyfragdict['H2']/2

        #print keyfragdict
        deltaTS = (keyfragdict['H2']*TS_H2) + (keyfragdict['CO2']*TS_CO2) + (keyfragdict['H2O']*TS_H2O)

        deltaTS_ts = 0
        for key,value in keyfragdict.items():
            if value > 0:
                deltaTS_ts = deltaTS_ts + keyTSdict[key]*value   

        #print deltaTS,deltaTS_ts

        return deltaTS,deltaTS_ts

    def GenProduct_fixedsurface(self,selectsurface ):
        for iStr,Str in enumerate(self):
            self[iStr].product = []
            #print Str.predictpatt
            for i,pattinfo in enumerate(Str.predictpatt[:15]):
                #print 'try patt %d'%i
                rp = pattinfo[0]
                if rp.ECFP0name == '626280626280640619640619':
                    continue
                infolist = []

                psl=self.OutProduct(Str,rp,depth =3)
                if len(psl) == 0:
                    #print 'wrong rp: fail to match rc'
                    continue
                #print (rp.ECFP0name,len(self.allrpdict0[rp.ECFP0name]))

                #_tmp = allstr()
                for ECFPname_withsite in self.allrpdict0[rp.ECFP0name][:15]:
                    #self.allrpdict_withsite_pairexample[ECFPname_withsite][0][0].label = ECFPname_withsite
                    #_tmp.append(self.allrpdict_withsite_pairexample[ECFPname_withsite][0][0])
                    #_tmp.append(self.allrpdict_withsite_pairexample[ECFPname_withsite][0].TSstr)
                    #_tmp.append(self.allrpdict_withsite_pairexample[ECFPname_withsite][0][1])
                    rp = self.allrpdict_withsite[ECFPname_withsite]
                    if int(rp.surfaceinfo.split('-')[0]) != selectsurface:
                        print rp.surfaceinfo.split('-')[0]
                        continue

                    if (self.mode)%2 == 1:
                        barrier,heat = self.RunRpinfoNet(rp,Str,self.rpmodel)
                    if (self.mode)%2 == 0:
                        barrier,heat = self.RunRpinfoNet_twomodel(rp,Str,self.barriermodel,self.heatmodel)
                    #print rp.ECFPname_withsite,rp.rangeinfo.barriermax,rp.rangeinfo.barriermin, rp.rangeinfo.heatmax,rp.rangeinfo.heatmin
                    if (barrier < 0) or (barrier < heat):
                        print ('error: wrong barrier/heat')
                        continue
                    if (barrier > rp.rangeinfo.barriermax +0.15) or (barrier < rp.rangeinfo.barriermin -0.15):
                        print ('error: barrier out of range')
                        continue
                    if (heat > rp.rangeinfo.heatmax +0.15) or (heat < rp.rangeinfo.heatmin -0.15):
                        print ('error: heat out of range')
                        continue
                    infolist.append([barrier,heat,rp.surfaceinfo,rp.ECFPname_withsite])

                #_tmp.printall('%s.arc'%(rp.ECFP0name))
                if len(infolist) != 0:
                    infolist.sort(key = lambda x:x[0])
                    barrier = infolist[0][0]
                    heat    = infolist[0][1]
                    surfaceinfo = infolist[0][2]
                    rpname_withsite =  infolist[0][3]
                else:
                    continue

                for p in psl:
                    deltaTS, deltaTS_ts = self.AddGibbsEnergy(self[iStr],p)
                    p.barrier = barrier + deltaTS_ts
                    p.heat = heat + deltaTS
                    #p.barrier = barrier
                    #p.heat = heat

                    p.surfaceinfo = surfaceinfo
                    p.rpname = rpname_withsite
                    #p.barrierinfo = patt.barrierinfo
                    #p.heatinfo = patt.heatinfo
                    #p.barrier = np.mean(p.barrierinfo)
                    #p.heat = np.mean(p.heatinfo)
                    self[iStr].product.append([p,pattinfo[1]])


    def PredictProduct_Singlemodel(self,depth =3):
        self.rpdict= {}
        self.GenNNInput0(self.RealIndexToInput0,depth)
        self.RunPredictNN(self.model,self.PattIndexDict)
        self.GenProduct_Singlemodel()

        plist = []
        self.name2str= {}

        for Str in self:
            Str.scoredict = {}
            for p in Str.product:
                if p[0].ECFPname not in Str.scoredict.keys():
                    Str.scoredict[p[0].ECFPname]= p[1]
                    self.name2str[p[0].ECFPname] = p[0]
                else:
                    Str.scoredict[p[0].ECFPname]= Str.scoredict[p[0].ECFPname] + p[1]
            for key,value in Str.scoredict.items():
                pair = [Str,self.name2str[key]]
                plist.append([pair,value])
            plist.sort(key= lambda X:X[1], reverse = True)
            self.rpdict[Str.ECFPname]= plist



    def PredictProduct_Simple(self,depth =3):
        self.rpdict= {}
        self.GenNNInput0(self.RealIndexToInput0,depth)
        self.RunPredictNN(self.model,self.PattIndexDict)
        self.GenProduct_Simple()

        plist = []
        self.name2str= {}

        for Str in self:
            Str.scoredict = {}
            for p in Str.product:
                if p[0].ECFPname not in Str.scoredict.keys():
                    Str.scoredict[p[0].ECFPname]= p[1]
                    self.name2str[p[0].ECFPname] = p[0]
                else:
                    Str.scoredict[p[0].ECFPname]= Str.scoredict[p[0].ECFPname] + p[1]
            for key,value in Str.scoredict.items():
                pair = [Str,self.name2str[key]]
                plist.append([pair,value])
            plist.sort(key= lambda X:X[1], reverse = True)
            self.rpdict[Str.ECFPname]= plist

    def PredictProduct(self,depth =3, LaddZPE = 1):
        self.rpdict= {}
        self.GenNNInput0(self.RealIndexToInput0,depth)
        self.RunPredictNN(self.model,self.PattIndexDict)
        if (self.mode>= 5):
            self.GenProduct_mode5(LaddZPE =LaddZPE )
        else:
            self.GenProduct()

        plist = []
        self.name2str= {}

        for Str in self:
            Str.scoredict = {}
            for p in Str.product:
                if p[0].ECFPname not in Str.scoredict.keys():
                    Str.scoredict[p[0].ECFPname]= p[1]
                    self.name2str[p[0].ECFPname] = p[0]
                else:
                    Str.scoredict[p[0].ECFPname]= Str.scoredict[p[0].ECFPname] + p[1]
            for key,value in Str.scoredict.items():
                pair = [Str,self.name2str[key]]
                plist.append([pair,value])
            plist.sort(key= lambda X:X[1], reverse = True)
            self.rpdict[Str.ECFPname]= plist


#        for Str in self:
#            for p in Str.product:
#                pair = [Str,p[0]]
#                plist.append([pair,p[1]])
#            plist.sort(key= lambda X:X[1], reverse = True)
        return plist


    def PredictProduct_fixedsurface(self,selectsurface,depth =3, LaddZPE = 1):
        self.rpdict= {}
        self.GenNNInput0(self.RealIndexToInput0,depth)
        self.RunPredictNN(self.model,self.PattIndexDict)
        if (self.mode== 0):
            self.GenProduct_Singlemodel(selectsurface)
        elif (self.mode >= 5):
            self.GenProduct_mode5(selectsurface, LaddZPE =LaddZPE)
        else:
            self.GenProduct_fixedsurface(selectsurface)

        plist = []
        self.name2str= {}

        for Str in self:
            Str.scoredict = {}
            for p in Str.product:
                if p[0].ECFPname not in Str.scoredict.keys():
                    Str.scoredict[p[0].ECFPname]= p[1]
                    self.name2str[p[0].ECFPname] = p[0]
                else:
                    Str.scoredict[p[0].ECFPname]= Str.scoredict[p[0].ECFPname] + p[1]
            for key,value in Str.scoredict.items():
                pair = [Str,self.name2str[key]]
                plist.append([pair,value])
            plist.sort(key= lambda X:X[1], reverse = True)
            self.rpdict[Str.ECFPname]= plist

        return plist

    def LoadNNInfo(self,model,RealIndexToInput0, PattIndexDict,ECFPdepth):
        self.ECFPdepth = ECFPdepth
        self.RealIndexToInput0 = RealIndexToInput0
        self.model = model
        self.PattIndexDict = PattIndexDict


    def LoadRPinfoNet(self,mode, rpmodel,rpRealIndexToInput, RealIndexToInput_withsite, allrpdict0, allrpdict_withsite,allrpdict_withsite_pairexample = False):
        self.rpRealIndexToInput = rpRealIndexToInput
        self.rpmodel = rpmodel
        self.RealIndexToInput_withsite = RealIndexToInput_withsite
        self.allrpdict0 = allrpdict0
        self.allrpdict_withsite = allrpdict_withsite
        self.allrpdict_withsite_pairexample = allrpdict_withsite_pairexample
        self.mode = mode

    def Loadtwomodel(self,mode, barriermodel,heatmodel,rpRealIndexToInput, RealIndexToInput_withsite, allrpdict0, allrpdict_withsite,allrpdict_withsite_pairexample = False):
        self.rpRealIndexToInput = rpRealIndexToInput
        self.barriermodel = barriermodel
        self.heatmodel = heatmodel
        self.RealIndexToInput_withsite = RealIndexToInput_withsite
        self.allrpdict0 = allrpdict0
        self.allrpdict_withsite = allrpdict_withsite
        self.allrpdict_withsite_pairexample = allrpdict_withsite_pairexample
        self.mode = mode

    def Connect_Path(self, name1,name2,ppath,depth, maxdepth = 10):
        if depth > maxdepth :
            return
        pmid=[]
        pmidpath =[]

        try :
            #rppair = self.rpdict[name1.sminame]
            rppair = self.rpdict[name1.ECFPname]
        except:

#        #  for name problem
#        allrecal = 1
#        if (allrecal):
            _tmp = AllReactStr()
            #_tmp.append(CanStr(name1))
            _tmp.append(name1)
            _tmp.LoadNNInfo(self.model,self.RealIndexToInput0,self.PattIndexDict,self.ECFPdepth)
            _tmp.LoadRPinfoNet(self.rpmodel,self.rpRealIndexToInput, self.RealIndexToInput_withsite, self.allrpdict0, self.allrpdict_withsite)
            #rppair =_tmp.PredictProduct()
            _tmp.PredictProduct()
            self.rpdict[name1.ECFPname] =  _tmp.rpdict[name1.ECFPname]
            rppair = self.rpdict[name1.ECFPname]
#            if name1.ECFPname not in self.ECFP2smi.keys()
#                _tmp.GetAllsminame_fromreactstr()
#                self.ECFP2smi[name1.ECFPname] = _tmp[0].sminame


        for i,pairinfo in enumerate(rppair):
            pair = pairinfo[0]
            #pair.score = pairinfo[1]
            if pair[1].ECFPname == name2.ECFPname:
                _ppath = Path()
                _ppath.copy(ppath)
                _ppath.append(pair)
                _ppath.namechain.append(pair[1].ECFPname)
                _ppath.score = ppath.score* pairinfo[1]
                #_ppath.score =min( ppath.score,pairinfo[1])
                self.allpath.append(_ppath)
            elif i < 15:
                if (pair[1].ECFPname not in ppath.namechain ):
                    pmid.append(pair[1])
                    _ppath = Path()
                    _ppath.copy(ppath)
                    _ppath.append(pair)
                    _ppath.namechain.append(pair[1].ECFPname)
                    _ppath.score = ppath.score* pairinfo[1]
                    #_ppath.score =min( ppath.score,pairinfo[1])
                    pmidpath.append(_ppath)
                    self.searchtrace.append(pair[1].ECFPname)

        depth =depth +1
        #piname2 = [name2]*len(pmid)
        #pdepth = [depth]*len(pmid)
        if len(pmid) > 0:
            pname2 = [name2]*len(pmid)
            pdepth = [depth]*len(pmid)
            pmaxdepth = [maxdepth]*len(pmid)
            return map(self.Connect_Path, pmid,pname2,pmidpath,pdepth,pmaxdepth)
        else :
            return



    def FindPath(self,Str1,Str2,maxdepth = 10, printstring = False):
        Str1.GenECFPname()
        Str2.GenECFPname()


        self.rpdict = {}
        self.allpath =[]
        self.searchtrace =[Str1.ECFPname]
        ppath =Path()
        ppath.namechain=[Str1.ECFPname]
        ppath.score = 1
        depth =0
        self.Connect_Path(Str1,Str2,ppath,depth, maxdepth)
        #self.allpath.sort(key = lambda X:len(X))
        for path in self.allpath:
            path.GetBarrier_list()


        #self.allpath.sort(key = lambda X:X.score, reverse= True)
        self.allpath.sort(key = lambda X:X.score)

        inputfile = open('allpath.pkl','wb')
        pickle.dump(self.allpath,inputfile)
        inputfile.close()


        tmpall= AllStr()
        for i,path in enumerate(self.allpath[:10]):
            for pair in path:
                tmpall.append(pair[0])
                tmpall.append(pair[1])
        tmpall.GetAllsminame_fromreactstr(numproc = 2)

        f= open('result','w')
        iStr =0 
        for i,path in enumerate(self.allpath[:10]):
            f.write('-------path-%d-------\n'%(i+1))
            f.write('score-%f----barrier--%f\n'%(path.score,path.barrier))
            for pair in path:
                f.write('%s   %s   %s   %s   %s   %s\n' %(tmpall[iStr+1].surfaceinfo,tmpall[iStr+1].rpname,tmpall[iStr].sminame,tmpall[iStr+1].sminame,pair[1].barrier,pair[1].heat))
                iStr =iStr+2
