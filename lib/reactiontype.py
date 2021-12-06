from Arc2RPdata import AllPair
from dataprint import *
from simple_func_k  import Dumppkl
import numpy as np

test = AllPair()

test.Readfile_Simple('allpair.arc',filemode= 2)

dataLabel = []
datarc = []
for pair in test:
    print pair.ISreacttype
    dataLabel.append(pair.ISreacttype)
    dataLabel.append(pair.FSreacttype)
    datarc.append(pair.rclist0)
    datarc.append(pair.rclist0)
OutArray_String(dataLabel, 'reacttype')
OutArray_Int(datarc,'rclist0')


test2 = AllPair()
test2.ReadLabel('testdata')

for ipair,pair in enumerate(test2):
    test2[ipair].ZPE_TS_IS2FS   = test[ipair].ZPE_TS_IS2FS   
    test2[ipair].ZPE_heat_IS2FS = test[ipair].ZPE_heat_IS2FS 
    test2[ipair].ZPE_TS_FS2IS   = test[ipair].ZPE_TS_FS2IS   
    test2[ipair].ZPE_heat_FS2IS = test[ipair].ZPE_heat_FS2IS 
    #test2[ipair].ZPE_TS_IS2FS   = 0
    #test2[ipair].ZPE_heat_IS2FS = 0
    #test2[ipair].ZPE_TS_FS2IS   = 0
    #test2[ipair].ZPE_heat_FS2IS = 0

test2.RPdict0_Sampled()

