# $Id: pyinit.py,v 1.6 2010/11/27 18:20:40 samn Exp $ 

from neuron import *

#h.xwindows=1
#h.show_panel=0
#h.pwd()
#h("load_file(\"xgetargs.hoc\")")
#h("load_file(\"grvec.hoc\")")
#h("load_file(\"labels.hoc\")")
#h("load_file(\"syncode.hoc\")")
#
#h("load_file(\"decnqs.hoc\")")
#h("if (! VECST_INSTALLED) install_vecst()")
#h("if (! INSTALLED_stats) install_stats()")
#
#h("install_PLACE()")
#
#h("load_file(\"nqs_utils.hoc\")")
#h("load_file(\"drline.hoc\")")
#h("load_file(\"stats.hoc\")")
#h("load_file(\"infot.hoc\")")
#h("load_file(\"decmat.hoc\")")

my_pylab_loaded = False

from math import sqrt, pi

import numpy

try:
    import pylab
    from pylab import plot, arange, figure
    my_pylab_loaded = True
except ImportError:
    my_pylab_loaded = False

import os
import datetime
import shutil

import pickle

