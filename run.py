# $Id: run.py,v 1.53 2010/12/15 22:20:27 samn Exp $ 

from pyinit import *
from geom import *
from networkmsj import *
from params import *
import sys

# sets up external inputs
if net.noise:
	net.set_noise_inputs(h.tstop) #h.tstop sets duration of inpus for make noise case

# handler for printing out time during simulation run
def fi():
	for i in range(0,int(h.tstop),100):
		h.cvode.event(i, "print " + str(i))

fih = h.FInitializeHandler(1, fi)

# initialize random # generators of NetStims - forces it at beginning of each sim
def myInitNetStims():
	#for i in range(19):
	#	print i,net.pyr.cell[i].soma.v,net.pyr.cell[i].Adend3.v,net.pyr.cell[i].Bdend.v
	#for i in range(19):
	#	print i,net.olm.cell[i].soma.v
	#for i in range(19):
	#	print i,net.bas.cell[i].soma.v
	net.init_NetStims()
		
fihns = h.FInitializeHandler(0, myInitNetStims)
# fihns = h.FInitializeHandler(1, myInitNetStims)

# handler for washin/washout
fiwash = None

olmWash =  [0, 0] # olm NMDA value for washin/washout
basWash =  [0, 0] # basket NMDA value for washin/washout
pyrWashA = [0, 0] # ...
pyrWashB = [0, 0] # ...
washinT  = 0      # washin time
washoutT = 0      # washout time

def dowashin():
	print "washIN at ", washinT, " = ", h.t , " ", olmWash[0], basWash[0], pyrWashB[0], pyrWashA[0]
	net.olm.set_r("somaNMDA",olmWash[0])
	net.bas.set_r("somaNMDA",basWash[0])
	net.pyr.set_r("BdendNMDA",pyrWashB[0])
	net.pyr.set_r("Adend3NMDA",pyrWashA[0])

def dowashout():
	print "washOUT at ", washoutT, " = " , h.t, " ", olmWash[1], basWash[1], pyrWashB[1], pyrWashA[1]	
	net.olm.set_r("somaNMDA",olmWash[1])
	net.bas.set_r("somaNMDA",basWash[1])
	net.pyr.set_r("BdendNMDA",pyrWashB[1])
	net.pyr.set_r("Adend3NMDA",pyrWashA[1])

def setwash():
	print "washinT ", washinT, " washoutT ", washoutT
	h.cvode.event(washinT,"nrnpython(\"dowashin()\")")
	h.cvode.event(washoutT,"nrnpython(\"dowashout()\")")

# example to do washin/washout, after loading sim:
# import run
# h.tstop=100
# run.olmWash =  [0, 1]
# run.basWash =  [1, 1]
# run.pyrWashA = [1, 1]
# run.pyrWashB = [1, 1]
# run.washinT  = 30
# run.washoutT = 60
# fiwash = h.FInitializeHandler(1,setwash)
# h.run()

class Power():
	pass

class Batch:
	def __init__(self,net):
		self.net = net    #the network, cells, synapses, etc.
		self.pow = Power() #the data
		
	def copydata(self,obj):
		self.pow.n     = obj.n
		self.pow.x     = obj.x
		self.pow.timer = obj.timer
		self.pow.tp    = obj.tp
		self.pow.gp    = obj.gp
		self.pow.tf    = obj.tf
		self.pow.gf    = obj.gf
		self.pow.arch  = obj.arch
		
	def save(self):
		file = open('filen.obj', 'w')
		pickle.dump(self.pow,file)
		
	def load(self):
		file = open('filen.obj', 'r')
		self.pow = pickle.load(file)

	#this function is based on loop in r function, to get a string for the sim params
	def getsimstr(self,r1,r2,r3,r4):
		simstr = "olm_somaNMDA_" + str(r1) + "_"
		simstr = simstr + "bas_somaNMDA_" + str(r2) + "_"
		simstr = simstr + "pyr_BdendNMDA_" + str(r3) + "_"
		simstr = simstr + "pyr_Adend3NMDA_" + str(r4) + "_"
		return simstr		
		
	def r(self, n):
		self.pow.n     = n
		self.pow.arch  = Archive()
		self.pow.x     = numpy.linspace(0,1,self.pow.n)
		self.pow.timer = h.startsw()
		self.pow.tp    = numpy.zeros((self.pow.n,self.pow.n,self.pow.n,self.pow.n))
		self.pow.gp    = numpy.zeros((self.pow.n,self.pow.n,self.pow.n,self.pow.n))
		self.pow.tf    = numpy.zeros((self.pow.n,self.pow.n,self.pow.n,self.pow.n))
		self.pow.gf    = numpy.zeros((self.pow.n,self.pow.n,self.pow.n,self.pow.n))
		
		for i1,r1 in enumerate(self.pow.x):
			self.net.olm.set_r("somaNMDA",r1)
			for i2, r2 in enumerate(self.pow.x):
				self.net.bas.set_r("somaNMDA",r2)
				for i3, r3 in enumerate(self.pow.x):
					self.net.pyr.set_r("BdendNMDA",r3)
					for i4, r4 in enumerate(self.pow.x):
						self.net.pyr.set_r("Adend3NMDA",r4)
						simstr = self.getsimstr(r1,r2,r3,r4)
						print "NMDA/AMPA: " + simstr
						h.run()
						print "Time: ", h.startsw() - self.pow.timer
					
						self.pow.arch.reset_time_stamp()
						
						self.net.calc_psd() #calculate lfp,psd and draw it, then save
						self.pow.arch.save_fig(3,simstr+"fft")
						
						self.net.rasterplot()#draw raster for all cells, then save
						self.pow.arch.save_fig(1,simstr+"rasterogram")

						self.pow.arch.save_vec(simstr+"lfp",self.net.vlfp) #save LFP Vector to file

						self.net.setrastervecs() #setup raster Vectors for ALL cells & save
						self.pow.arch.save_vec(simstr+"idvec",self.net.myidvec)
						self.pow.arch.save_vec(simstr+"timevec",self.net.mytimevec)
			
						self.pow.tp[i1,i2,i3,i4] = self.net.tp
						self.pow.gp[i1,i2,i3,i4] = self.net.gp
						self.pow.tf[i1,i2,i3,i4] = self.net.tf
						self.pow.gf[i1,i2,i3,i4] = self.net.gf

						
	def plot_r(self):		
		self.plot_fun(5, "tp","Theta_Power")
		self.plot_fun(6, "gp","Gamma_Power")
		self.plot_fun(7, "tf","Theta_Frequency")
		self.plot_fun(8, "gf","Gamma_Frequency")		
		self.plot_fun(9, "tp","Theta_Power_Mean",    1)
		self.plot_fun(10,"gp","Gamma_Power_Mean",    1)
		self.plot_fun(11,"tf","Theta_Frequency_Mean",1)
		self.plot_fun(12,"gf","Gamma_Frequency_Mean",1)
		
	def plot_fun(self, fig, var, ylabel, mode=0):
		cond = ["olm","bas","pyrB","pyrA3"]
		f = pylab.figure(fig)
		f.clf()
		f.canvas.mpl_connect('pick_event', onpick)
		
		pylab.subplot(2,2,1)
		pylab.xlabel("NMDA/AMPA for " + cond[0])
		pylab.ylabel(ylabel)
		if mode==0:
			for i1, r1 in enumerate(self.pow.x):
				for i2, r2 in enumerate(self.pow.x):
					for i3, r3 in enumerate(self.pow.x):
						print "[:,"+str(i1)+","+str(i2)+","+str(i3)+"]"
						pylab.plot(self.pow.x, self.pow.__dict__[var][:,i1,i2,i3],label="[:,"+str(i1)+","+str(i2)+","+str(i3)+"]", picker=1)
		#pylab.label()
		else:
			pylab.plot(self.pow.x, self.pow.__dict__[var].mean(axis=1).mean(axis=1).mean(axis=1))

		
		pylab.subplot(2,2,2)
		pylab.xlabel("NMDA/AMPA for " + cond[1])
		pylab.ylabel(ylabel)
		if mode==0:	
			for i1, r1 in enumerate(self.pow.x):
				for i2, r2 in enumerate(self.pow.x):
					for i3, r3 in enumerate(self.pow.x):
						pylab.plot(self.pow.x, self.pow.__dict__[var][i1,:,i2,i3],label="["+str(i1)+",:,"+str(i2)+","+str(i3)+"]", picker=1)
						
		#pylab.label()
		else:
			pylab.plot(self.pow.x, self.pow.__dict__[var].mean(axis=0).mean(axis=1).mean(axis=1))
			
		pylab.subplot(2,2,3)
		pylab.xlabel("NMDA/AMPA for " + cond[2])
		pylab.ylabel(ylabel)
		if mode==0:	
			for i1, r1 in enumerate(self.pow.x):
				for i2, r2 in enumerate(self.pow.x):
					for i3, r3 in enumerate(self.pow.x):
						pylab.plot(self.pow.x, self.pow.__dict__[var][i1,i2,:,i3],label="["+str(i1)+","+str(i2)+",:,"+str(i3)+"]", picker=1)
						
		#pylab.label()		
		else:
			pylab.plot(self.pow.x, self.pow.__dict__[var].mean(axis=0).mean(axis=0).mean(axis=1))
		
		
		pylab.subplot(2,2,4)
		pylab.xlabel("NMDA/AMPA for " + cond[3])
		pylab.ylabel(ylabel)
		if mode==0:	
			for i1, r1 in enumerate(self.pow.x):
				for i2, r2 in enumerate(self.pow.x):
					for i3, r3 in enumerate(self.pow.x):
						pylab.plot(self.pow.x, self.pow.__dict__[var][i1,i2,i3,:],label="["+str(i1)+","+str(i2)+","+str(i3)+",:]", picker=1)
						
		#pylab.label()			
		else:
			pylab.plot(self.pow.x, self.pow.__dict__[var].mean(axis=0).mean(axis=0).mean(axis=0))
		
		self.pow.arch.save_fig(fig,ylabel)

def onpick(event):
	print "REWR"
	print str(event.artist.get_label())+" ("+str(event.mouseevent.xdata)+","+str(event.mouseevent.ydata)+")"
	return True

#save vec to fn (fn is path)
def mysvvec(fn,vec):
	fp = h.File()
	fp.wopen(fn)
	if fp.isopen():
		vec.vwrite(fp)
		fp.close()
	else:
		print "savevec ERR: couldn't open " + fn

#this class is for saving output, i.e. figures and py files to backup	
class Archive:
	def __init__(self):
		self.figprefix = "./gif" #prefix for saving figures
		self.datprefix = "./data"
		self.pyprefix = "./backup"
		self.reset_time_stamp()
		self.save_pyfile("par_sim.py")
		self.save_pyfile("Cells.py")
		
	def save_fig(self, fig, name):
		fn = os.path.join(self.figprefix, self.time_stamp+name+".svg")
		pylab.figure(fig)
		pylab.savefig(fn, orientation='landscape', format='svg', dpi=72)
		
	def reset_time_stamp(self):
		ts = datetime.datetime.now().timetuple()
		self.time_stamp = "_"+str(ts.tm_year)+"_"+str(ts.tm_mon)+"_"+str(ts.tm_mday)+"_"+str(ts.tm_hour)+"_"+str(ts.tm_min)+"_"+str(ts.tm_sec)
		
	def save_pyfile(self, fn):
		nfn = os.path.join(self.pyprefix,fn+self.time_stamp+".py")
		shutil.copy(fn, nfn)

	def save_vec(self, fn, vec):
		nfn = os.path.join(self.datprefix,fn+".vec")
		mysvvec(nfn,vec)

#run a sim and save data
def minrunsv(simstr,tstop=1200,dt=0.1):
	h.tstop=tstop
	h.dt=dt
	h.run()
	print "saving output data"
	net.calc_lfp()
	fn = "./data/"+simstr+"_lfp.vec"
	mysvvec(fn,net.vlfp)
	net.setsnq() # make NQS with spike times
	fn = "./data/"+simstr+"_snq.nqs"
	net.snq.sv(fn)
	print "making and saving output figures"

#read a Vector from file, fn is file-path, vec is a Vector
def myrdvec(fn,vec):
	fp=h.File()
	fp.ropen(fn)
	if not fp.isopen():
		print "myrdvec ERRA: Couldn't open " + fn
		return False
	vec.vread(fp)
	fp.close()
	return True

#load data from minrunsv into net.vlfp,net.snq
def loadminrundat(simstr):
	fs = "./data/"+simstr+"_lfp.vec"
	try:
		net.vlfp.resize(0)
	except:
		net.vlfp = h.Vector()
	myrdvec(fs,net.vlfp)
	fs = "./data/"+simstr+"_snq.nqs"
	try:
		h.nqsdel(net.snq)		
	except:
		pass
	try:
		net.snq=h.NQS(fs)
	except:
		print "loadminrundat ERRB: couldn't read snq from " + fs
	net.snq.verbose=0 # next, copy snq into vectors so can plot with net.rasterplot
	for po in net.cells:
		for i in xrange(len(po.lidvec)):
			id = po.cell[i].id
			po.lidvec[i].resize(0)
			po.ltimevec[i].resize(0)
			if net.snq.select("id",id):
				po.lidvec[i].copy(net.snq.getcol("id"))
				po.ltimevec[i].copy(net.snq.getcol("t"))
	net.snq.verbose=1

def testrun():
	net.olm.set_r("somaNMDA",0)
	h.run()
	arch = Archive()
	net.rasterplot(1)
	arch.save_fig(1,"tmp_rasterplot")
	net.psr.cell[0].plot_volt("soma",2)
	arch.save_fig(2,"tmp_psr_soma_volt")
	net.calc_psd(3)
	arch.save_fig(3,"tmp_fft")
	print "\a"

def batchrun():
	bat = Batch(net)
	bat.r(3)
	bat.plot_r()

def myrast(spikes,times,sz=12):	
	if h.g[0] == None:
		h.gg()
	spikes.mark(h.g[0],times,"O",sz,1,1)
	h.g[0].exec_menu("View = plot")

# testsame - for debugging two runs to make sure output is the same
def testsame(ts,v1,v2):
	h.tstop = ts
	v1 = h.Vector()
	h.run()
	net.calc_lfp()
	v1.copy(net.vlfp)
	v2 = h.Vector()
	h.run()
	net.calc_lfp()
	v2.copy(net.vlfp)
	print "same = " , v1.eq(v2)

#gethilbnqs - make two NQS objects out of LFP with phase/amplitude/filered signals in theta and gamma bands
def gethilbnqs(vlfp,minth=3,maxth=12,ming=30,maxg=80,usemlab=True):
  sampr = 1e3/h.dt # sampling rate in Hertz
  if usemlab:
    nqtheta=h.mathilbert(vlfp,sampr,minth,maxth)
    nqgamma=h.mathilbert(vlfp,sampr,ming,maxg)
  else:
    nar = vlfp.to_python() # -> python -> numpy format
    nar = numpy.array(nar)
    nqtheta=filt.gethilbnq(nar,sampr,minth,maxth) # get an NQS with 'theta' 
    nqgamma=filt.gethilbnq(nar,sampr,ming,maxg)# get an NQS with 'gamma'
  return [nqtheta,nqgamma]

#getampphnq - get an nqs with gamma amplitude vs theta phase - uses NQS objects created by gethilbnqs
def getampphnq(nqtheta,nqgamma,phbins=100,skipms=200):
  colp = int(nqgamma.fi("phase")) # column index for phase
  cola = int(nqgamma.fi("amp"))   # column index for amp
  phmin=nqgamma.v[colp].min() # minimum phase of gamma
  phmax=nqgamma.v[colp].max() # maximum phase of gamma
  phrng=phmax-phmin # range of gamma phase
  nq = h.NQS("avgamp","phase","n","err","minamp","maxamp") # output nqs - amp is average amplitude for a phase, vn is # of samples @ the phase
  #minamp is avgamp - stderr, maxamp is avgamp + stderr. those columns just for easier display of avg+/-error
  vamp=nq.v[0] # average amplitude for a given phase
  vph=nq.v[1] # theta phase
  vn=nq.v[2] # number of samples at the given phase
  ve=nq.v[3] # stderror
  vmin=nq.v[4] # avg-stderror
  vmax=nq.v[5] # avg+stderror
  vph.indgen(phmin,phmax,phrng/phbins) # init range of phases
  nq.pad()
  vamp.fill(0)
  vn.fill(0) # init counts to 0
  lv = h.List() # list to keep amplitude samples
  for i in xrange(int(vph.size())):
    lv.append(h.Vector())
  sz=int(nqgamma.v[0].size())
  startx=int(skipms/h.dt)
  for i in xrange(startx,sz,1):
    bin=int(phbins*(nqtheta.v[colp][i]-phmin)/phrng)
    if bin<0:
      print "bin < 0!"
    if bin>=phbins+1:
      print "bin >= phbins+1"
    lv.o(bin).append(nqgamma.v[cola][i])
  for i in xrange(0,int(vamp.size()),1):
    sz = lv.o(i).size()
    if sz > 0: # if no samples, skip
      av = lv.o(i).mean()
      if sz > 1: # make sure can call stderr
        er = lv.o(i).stderr()
      else:
        er = 0
      vamp.x[i] = av
      vn.x[i] = sz
      ve.x[i] = er
      vmin.x[i] = av - er
      vmax.x[i] = av + er
  return nq

# checkbase - compares baseline to OLM activity off
# returns results in a python list
def checkbase(endt=3e3,skipms=200,justone=False):
  vlfp = []
  vtmp = []
  nqp = []
  nqa = []
  nqh = []
  snq = []
  fnq = []
  h.tstop=endt
  j = 0
  dt = h.dt
  for i in xrange(1,-1,-1):
    print "set olm NMDA to ", float(i)
    net.olm.set_r("somaNMDA",float(i))
    print "running for " , endt , " ms "
    h.run()
    net.calc_lfp()
    vlfp.append(net.vlfp)
    vtmp.append(h.Vector())
    vtmp[j].copy(vlfp[j],skipms/dt,vlfp[j].size()-1)
    vtmp[j].sub(vtmp[j].mean())
    nqp.append( h.matpmtm(vtmp[j],1e3/dt) )
    vtmp[j].copy(vlfp[j],skipms/dt,vlfp[j].size()-1)
    nqh.append( gethilbnqs(vtmp[j],3,12,30,80) )
    nqa.append( getampphnq(nqh[j][0],nqh[j][1]) )
    net.setsnq()
    snq.append( h.NQS() )
    snq[j].cp(net.snq)
    net.setfnq(skipms)
    fnq.append( h.NQS() )
    fnq[j].cp(net.fnq)
    net.pravgrates()
    j += 1
    if justone:
	    break
  return [vlfp,vtmp,nqp,nqa,nqh,snq,fnq]

############################
#   setup multithreading   #
pc = h.ParallelContext()   #
pc.nthread(32)             #

#h.load_file('parcom.hoc')
#pc = h.ParallelComputeTool()
#pc.nthread(8)
#p.multisplit(True)

############################

if 0:
	testrun()

if 0:
	h.tstop=200
	net.pyr.cell[0].set_spikes([100],"BdendNMDA", 28*0.04e-3)
	h.run()
	net.pyr.cell[0].plot_volt("soma")

if 0:
	batchrun()
####################################################################################################
