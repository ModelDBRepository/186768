# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 06:23:41 2022

@author: maria
"""
from neuron import h
import pickle    
with open('netresults.pickle', 'rb') as f:
    allvars = pickle.load(f)
data = allvars[0]
lfp = allvars[1]
dt = allvars[2]
import savedata as sd
sd.graphrasterplot(data)   

sd.graphlfp(lfp,dt) 


