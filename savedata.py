# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 06:00:21 2022

@author: maria
"""
from neuron import h

def data4rasterplot(net,sz=2):

    listocells = []
    for po in net.cells:
        id = h.Vector()
        tv = h.Vector()
        for i in range(po.n):
            id.append(po.lidvec[i])
            tv.append(po.ltimevec[i])
        data = {"id":id.to_python(), "tv":tv.to_python()}
        listocells.append(data)
    return listocells
      
  
def graphrasterplot(rasterdata,sz=2):
    import matplotlib.pyplot as plt

    fig = plt.figure()
    pon  = 0        
    # if h.g[0] == None:
    col = ['r','g','b','o','k']
    for po in rasterdata:
        id = h.Vector()
        tv = h.Vector()
        id.from_python(po['id'])
        tv.from_python(po['tv'])
        if len(tv) > 0:
            plt.plot(tv,id, linestyle="None",marker='o',markersize=sz, color=col[pon])
        pon += 1

def data4pravgrates(self,skipms=200):
    try:
        self.fnq.tog("DB")
    except:
        self.setfnq(skipms)
    ty = 0
    tf = float( h.tstop - skipms )
    for po in self.cells:
        self.fnq.select("ty",ty)
        vf = self.fnq.getcol("freq")
        if vf.size() > 1:
            print("ty: "), ty, " avg rate = ", vf.mean(), "+/-", vf.stderr(), " Hz"
        else:
            print("ty: "), ty, " avg rate = ", vf.mean(), "+/-", 0.0 , " Hz"
        ty += 1

def graphlfp(lfp,dt):
    import matplotlib.pyplot as plt
    import numpy as np
    tvec = np.arange(0,dt*len(lfp),dt)
    figure = plt.figure()
    plt.plot(tvec,lfp)