import numpy as np
import pickle

def Readarray_Int(inputfile):
    datain = []
    f= open(inputfile,'r')
    iArray = -1
    for line in f:
        if ('Array' in line):
            iArray = iArray +1 
            datain.append([])
        if ('datain' in line):
            for x in line.split()[1:]:
                datain[iArray].append(int(x))
    
    f.close()
    return datain

def Readarray_Float(inputfile):
    datain = []
    f= open(inputfile,'r')
    iArray = -1
    for line in f:
        if ('Array' in line):
            iArray = iArray +1
            datain.append([])
        if ('datain' in line):
            for x in line.split()[1:]:
                datain[iArray].append(float(x))

    f.close()
    return datain



def OutArray_Int(data,outfile):
    f = open(outfile,'w')

    for iarray,array in enumerate(data):
        #print array
        #print len(array)
        f.write('Array   %8d\n'%iarray)
        for i in xrange(0,len(array),16):
            #f.write('datain %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\n' %(array[i:(i+16)]))
            f.write('datain %s\n' %reduce(lambda x,y : x+y, ['%4d' %num for num in array[i:(i+16)]]))
        #if len(array) > (i+1):
        #    f.write('datain %s\n' %reduce(lambda x,y : x+y, ['%4d' %num for num in array[i:len(array)]]))


    f.close()
    return


def OutArray_Float(data,outfile):
    f = open(outfile,'w')

    for iarray,array in enumerate(data):
        #print array
        #print len(array)
        f.write('Array   %8d\n'%iarray)
        for i in xrange(0,len(array),4):
            #f.write('datain %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\n' %(array[i:(i+16)]))
            f.write('datain %s\n' %reduce(lambda x,y : x+y, ['%16.10f' %num for num in array[i:(i+4)]]))
        #if len(array) > (i+1):
        #    f.write('datain %s\n' %reduce(lambda x,y : x+y, ['%4d' %num for num in array[i:len(array)]]))


    f.close()
    return

def OutArray_String(data,outfile):
    f = open(outfile,'w')

    for iarray,array in enumerate(data):
        f.write('Array   %8d\n'%iarray)
        f.write('datain %s\n' %array)

    f.close()
    return

def Readarray_String(inputfile):
    datain = []
    f= open(inputfile,'r')
    iArray = -1
    for line in f:
        if ('Array' in line):
            iArray = iArray +1
            datain.append('')
        if ('datain' in line):
            for x in line.split()[1:]:
                datain[iArray]= datain[iArray]+' '+x

    f.close()
    return datain



if __name__ == "__main__":
    datain=open('datainput_info.pkl','rb')
    data=pickle.load(datain)
    datain.close()
    
    dataout=open('dataout_info.pkl','rb')
    labels=pickle.load(dataout)
    dataout.close()
    
    
    OutArray_Int(data,'datainput_info')
    OutArray_Float(labels,'dataout_info')

    datain=open('datain.pkl','rb')
    data=pickle.load(datain)
    datain.close()

    dataout=open('dataresult.pkl','rb')
    labels=pickle.load(dataout)
    dataout.close()


    OutArray_Int(data,'datainput')
    OutArray_Float(labels,'dataout')

    
    #datain=Readarray_Float('test1')
    #dataout=Readarray_Float('test2')
