import tensorflow as tf
from tensorflow import keras
import numpy as np
import pickle
import time
from dataprint import *
from metrics import MaxError_barrier,MaxError_heat


data1=Readarray_Int('datain_info')
data2=Readarray_Float('dataresult_info')

data = np.array(data1)
labels = np.array(data2)

data_add=Readarray_Float('datain_addinfo')
dataadd = np.array(data_add)
dataset = []
for i in range(data.shape[0]):
    _array = np.append(data[i],dataadd[i])
    dataset.append(_array)
dataset = np.array(dataset)

sampleweight =Readarray_Float('datalabel')
sampleweight = np.array(sampleweight)


#dataset = data

#print dataset
#print labels

permutation = np.random.permutation(dataset.shape[0])

shuffled_data = dataset[permutation,:]
shuffled_label = labels[permutation,:]
shuffled_weight = sampleweight[permutation,0]

print shuffled_data
print shuffled_label

data = shuffled_data
labels = shuffled_label

#data = dataset
#shuffled_weight = sampleweight[:,0]

print data.shape[0]
print data.shape[1]
print labels.shape[1]


#labels = labels[:,0]

#print labels


model = keras.Sequential([
keras.layers.Dense(512, activation='sigmoid',input_shape=(data.shape[1],)),
keras.layers.Dropout(0.2),
#keras.layers.Dense(128, activation='elu',input_shape=(data.shape[1],)),
keras.layers.Dense(256, activation='sigmoid'),
#keras.layers.Dense(128, activation='sigmoid'),
keras.layers.Dropout(0.2),
#keras.layers.Dense(128, activation='sigmoid'),
#keras.layers.Dropout(0.2),
keras.layers.Dense(labels.shape[1], activation='linear')
#keras.layers.Dense(1, activation='linear')
])



model.compile(optimizer=tf.train.AdamOptimizer(0.001),
                loss='mean_squared_error',
                 metrics=['mae',MaxError_barrier(), MaxError_heat()])

#traindata = data
#trainlabels= labels

traindata = data[:]
trainlabels= labels[:]

val_data =  data[:]
val_labels  = labels[:]

t1 =time.clock()
_t1 = time.time()

for i in range(3):

    hist= model.fit(traindata, trainlabels,sample_weight = shuffled_weight[:], epochs=200, batch_size=32,
              validation_data=(val_data, val_labels))
    #print (hist.history)
    t2 =time.clock()
    _t2 = time.time()

    loss_and_metrics = model.evaluate(data,labels, batch_size= 32)

    print ('wall time: %f'%(_t2-_t1))
    print ('cpu  time: %f'%(t2-t1))
    print ('epoch %d end'%(200*(i+1)))
    model.save('modelsave_%d.h5'%(200*(i+1)))


model.save('infonet.h5')
