# $Id: params.py,v 1.17 2010/12/03 16:21:57 samn Exp $ 

from pyinit import *
from geom import *
from networkmsj import *

# run params
h.tstop=3e3
h.dt = 0.1
h.steps_per_ms = 1/h.dt
#h.cvode_local(1)
h.v_init = -65


# network NMDA params
olmSomaNMDA=1
bassomaNMDA=1
pyrBdendNMDA=1
pyrAdend3NMDA=1

def gSetNMDA(net):
    net.olm.set_r("somaNMDA",olmSomaNMDA)
    net.bas.set_r("somaNMDA",bassomaNMDA)
    net.pyr.set_r("BdendNMDA",pyrBdendNMDA)
    net.pyr.set_r("Adend3NMDA",pyrAdend3NMDA)

gSetNMDA(net)

