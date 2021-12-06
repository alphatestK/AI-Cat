import tensorflow as tf
from tensorflow import keras
import numpy as np


#def MaxError_metric_fn(y_true, y_pred):
#    absE = tf.abs(y_true - y_pred)
#    #return tf.reduce_max(absE, axis = -1)
#    return tf.reduce_max(absE)

#def MaxError_barrier(y_true, y_pred):
#    absE = tf.abs(y_true[:,0] - y_pred[:,0])
#    return tf.reduce_max(absE)
#
#def MaxError_heat(y_true, y_pred):
#    absE = tf.abs(y_true[:,1] - y_pred[:,1])
#    return tf.reduce_max(absE)


class MaxError_barrier(tf.keras.metrics.Metric):
    def __init__(self, name ='maxerror_barrier', **kwargs):
        super(MaxError_barrier, self).__init__(name=name, **kwargs)
        self.maxerror = self.add_weight(name='tttt', initializer='zeros')

    def update_state(self,  y_true, y_pred, sample_weight=None):
        t1 = tf.subtract(0.0,self.maxerror)
        t2 = tf.add(0.0,self.maxerror)
        #t = self.maxerror
        #print (t)
        #print (np.array(t))
        #print (t.value)
        #self.maxerror.assign_add(-t)
        absE = tf.abs(y_true[:,0] - y_pred[:,0])
        _tmpmax = tf.reduce_max(absE)
        nowmax = tf.reduce_max([_tmpmax,t2])
        add = tf.add(nowmax,t1)
        #self.maxerror.assign_add(nowmax)
        self.maxerror.assign_add(add)

    def result(self):
        return self.maxerror

# if tf version is 2.0, this part will be needed
#    def reset_states(self):
#        self.maxerror.assign(0)


class MaxError_heat(tf.keras.metrics.Metric):
    def __init__(self, name ='maxerror_heat', **kwargs):
        super(MaxError_heat, self).__init__(name=name, **kwargs)
        self.maxerror = self.add_weight(name='tttt', initializer='zeros')

    def update_state(self,  y_true, y_pred, sample_weight=None):
        t1 = tf.subtract(0.0,self.maxerror)
        t2 = tf.add(0.0,self.maxerror)
        #t = self.maxerror
        #print (t)
        #print (np.array(t))
        #print (t.value)
        #self.maxerror.assign_add(-t)
        absE = tf.abs(y_true[:,1] - y_pred[:,1])
        _tmpmax = tf.reduce_max(absE)
        nowmax = tf.reduce_max([_tmpmax,t2])
        add = tf.add(nowmax,t1)
        #self.maxerror.assign_add(nowmax)
        self.maxerror.assign_add(add)

    def result(self):
        return self.maxerror




