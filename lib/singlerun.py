import tensorflow as tf
from tensorflow import keras
import numpy as np
import pickle
#from reactpair import Pair
from allreactstr import AllReactStr
from allstr_new import allstr
from ECFP0 import AllStr

from MCTS import Node, Loadinfo

#inputfile = 'methanol.arc'
inputfile = 'test.arc'
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


NNset = Loadinfo()
predicttest = Node(test[0], NNset)



pl =predicttest.expand()

testall = AllReactStr()
for p in pl:
    print p,p[0][1].barrier,p[0][1].heat
    testall.append(p[0][1])
testall.GetAllsminame_fromreactstr()
for Str in testall:
    print Str.sminame

stdpatt =  stdout[test[0].ECFPname]
testout = []
for ikey,e in enumerate(stdpatt):
    if e > 1e-10:
        print ikey,e
        testout.append([PattIndexDict[ikey],e])
    testout.sort(key = lambda X:X[1],reverse = True)



product =AllStr()
for pattinfo in testout:
    patt = pattinfo[0]
    psl=testall.OutProduct(test[0],patt,depth =3)
    for p in psl:
        print(p.ECFPname,pattinfo[1])
        product.append(p)
    
product.GetAllsminame_fromreactstr()

for Str in product:
    print Str.sminame
