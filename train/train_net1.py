import tensorflow as tf
from tensorflow import keras
import numpy as np
import pickle
import time
from dataprint import *


data1=Readarray_Int('datain')
data2=Readarray_Float('dataresult')

data = np.array(data1)
labels = np.array(data2)



print data.shape[0]
print data.shape[1]
print labels.shape[1]



model = keras.Sequential([
keras.layers.Dense(512, activation='elu',input_shape=(data.shape[1],)),
keras.layers.Dropout(0.2),
#keras.layers.Dense(128, activation='elu',input_shape=(data.shape[1],)),
keras.layers.Dense(512, activation='elu'),
keras.layers.Dropout(0.2),
#keras.layers.Dense(128, activation='elu'),
#keras.layers.Dropout(0.2),
keras.layers.Dense(labels.shape[1], activation='softmax')
])


model.compile(optimizer=tf.train.AdamOptimizer(0.001),
              loss='categorical_crossentropy',
              metrics=['accuracy'])


#traindata = data[:700]
#trainlabels= labels[:700]

traindata = data[:]
trainlabels= labels[:]
print traindata
print trainlabels

val_data =  data[:]
val_labels  = labels[:]


t1 = time.clock()
_t1 = time.time()

for i in range(5):

    model.fit(traindata, trainlabels, epochs=500, batch_size=100,
              validation_data=(val_data, val_labels))
    
    t2 =time.clock()
    _t2 = time.time()
    
    print ('wall time: %f'%(_t2-_t1))
    print ('cpu  time: %f'%(t2-t1))
    
    model.save('modelsave_%d.h5'%(500*i))

model.save('policymodel.h5')

