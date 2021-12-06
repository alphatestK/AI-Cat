
import tensorflow as tf
from tensorflow import keras
import numpy as np
from Arc2RPdata import RangeInfo
import pickle
#from reactpair import Pair
from allreactstr import AllReactStr
from allstr_new import allstr
from ECFP0 import AllStr
from simple_func_k import Loadpkl,Dumppkl

#model = keras.models.load_model('testmodel.h5')
model = keras.models.load_model('policymodel.h5')

#model2= keras.models.load_model('net2.h5')
model2= keras.models.load_model('infonet.h5')

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


inputfile = 'methanol.arc'
#inputfile = 'Pttest.arc'
#inputfile = 'lowestpath-1'
#inputfile = 'Formate-Cads.arc'

test= AllReactStr()
test.readfile(inputfile)
test.GetAllsminame()

tmp= allstr()
tmp.readfile(inputfile)
tmp.GetAllsminame()



for i,Str in enumerate(test):
    print Str.sminame
    Str.sminame = tmp[i].sminame
    print Str.sminame

ndim=3

predicttest= AllReactStr()
predicttest.LoadNNInfo(model,RealIndexToInput_0,PattIndexDict,ndim)
predicttest.LoadRPinfoNet(model2,rpRealIndexToInput, RealIndexToInput_withsite,allrpdict0, allrpdict_withsite)

test[0].GenECFPname()
predicttest.append(test[0])
#pl =predicttest.PredictRppair()
#pl =predicttest.PredictProduct()

predicttest.FindPath(test[0],test[1])
#test[0].GenECFPname()
#predicttest.append(test[0])
##pl =predicttest.PredictRppair()
#pl =predicttest.PredictProduct()
#
#testall = AllReactStr()
#for p in pl:
#    print p
#    testall.append(p[0][1])
#testall.GetAllsminame_fromreactstr()
#for Str in testall:
#    print Str.sminame
#
#stdpatt =  stdout[test[0].ECFPname]
#testout = []
#for ikey,e in enumerate(stdpatt):
#    if e > 1e-10:
#        print ikey,e
#        testout.append([PattIndexDict[ikey],e])
#    testout.sort(key = lambda X:X[1],reverse = True)
#
#
#
#product =AllStr()
#for pattinfo in testout:
#    patt = pattinfo[0]
#    psl=testall.OutProduct(test[0],patt,depth =3)
#    for p in psl:
#        print(p.ECFPname,pattinfo[1])
#        product.append(p)
#    
#product.GetAllsminame_fromreactstr()
#
#for Str in product:
#    print Str.sminame
