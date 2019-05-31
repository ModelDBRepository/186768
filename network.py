# $Id: network.py,v 1.125 2011/06/10 15:10:05 samn Exp $ 

from pyinit import *
from geom import *
import random

gGID = 0 # global ID for cells

class Population:
	"Population of cells"
	# cell_type -- pyr, bas, olm
	# n -- number of cells in the population
	# x, y, z -- initial position for the first Cell
	# dx -- an increment of the x-position for the cell location
	# amp, dur, delay -- parameters for the IClamp in the soma
	# Spikes are stored in ltimevec (times) and lidvec (cell # within the population)
	def __init__(self, cell_type, n, x, y, z, dx, amp, dur, delay):
		global gGID
		self.cell = [] # List of cells in the population
		self.nc   = [] # NetCon list for recording spikes
		self.n    = n  # number of cells
		self.x    = x
		self.y    = y
		self.z    = z
		self.ltimevec = h.List() # list of Vectors for recording spikes, one per cell
		self.lidvec = h.List()
		self.nssidx = {}
		self.nseidx = {}
		self.ncsidx = {}
		self.nceidx = {}
		for i in range(n):
			self.cell.append(cell_type(x+i*dx,y,z,gGID))
			self.cell[-1].somaInj.amp   = amp
			self.cell[-1].somaInj.dur   = dur
			self.cell[-1].somaInj.delay = delay
			self.nc.append(h.NetCon(self.cell[-1].soma(0.5)._ref_v, None, sec=self.cell[-1].soma))
			self.ltimevec.append(h.Vector()) #NB: each NetCon gets own Vectors for recording. needed to avoid multithreading crash
			self.lidvec.append(h.Vector())
			self.nc[-1].record(self.ltimevec[-1],self.lidvec[-1],gGID) # record cell spikes with gGID
			gGID = gGID + 1 # inc global cell ID
			
	def set_r(self, syn, r):
		for c in self.cell:
			c.__dict__[syn].syn.r = r

class MSpec: # this class uses matlab to make a spectrogram

	def __init__(self,vlfp,maxfreq,nsamp,dodraw): #make a spectrogram using matlab
		h("jjj=name_declared(\"nqspec\")")
		h("if(jjj){nqsdel(nqspec) print \"deleted nqspec\"}")
		h("objref nqspec")
		vslfp = h.Vector()
		vslfp.copy(vlfp)
		vslfp.sub(vlfp.mean())
		h.nqspec = h.matspecgram(vslfp,1e3/h.dt,maxfreq,nsamp,dodraw)
		self.nqspec = h.nqspec

	def powinrange(self,minf,maxf): # get scalar power in range of frequencies
		nn = self.nqspec.select(-1,"f","[]",minf,maxf)
		if nn == 0:
			return 0
		h("jnk = 0")
		h("vec.resize(0)")
		for i in self.nqspec.ind:
			mystr = "vec.copy(nqspec.get(\"pow\","
			mystr += str(int(i))
			mystr += ").o)"
			h(mystr)
			h("jnk += vec.sum()")
		jnk = h.jnk
		return jnk / nn

	def vecinrange(self,minf,maxf): # get vector of power in range of frequencies vs time
		nn = self.nqspec.select(-1,"f","[]",minf,maxf)
		if nn == 0:
			return None
		h("objref vjnk")
		h("vjnk=new Vector()")
		h.vec.resize(0)
		for i in self.nqspec.ind:
			mystr = "vec.copy(nqspec.get(\"pow\","
			mystr += str(int(i))
			mystr += ").o)"
			h(mystr)
			if h.vjnk.size()==0:
				h.vjnk.copy(h.vec)
			else:
				h.vjnk.add(h.vec)
		h.vjnk.div(self.nqspec.ind.size())
		return h.vjnk

class Network:

	def __init__(self,noise=True,connections=True,DoMakeNoise=True,iseed=1234,UseNetStim=True,wseed=4321,scale=1.0,MSGain=1.0,SaveConn=False):
		import math
		print "Setting Cells"
		self.pyr = Population(cell_type=PyrAdr,n=int(math.ceil(800*scale)), x= 0, y=0, z=0, dx=50, amp= 50e-3, dur=1e9, delay=2*h.dt)
		self.bas = Population(cell_type=Bwb,   n=int(math.ceil(200*scale)), x=10, y=0, z=0, dx=50, amp=     0, dur=  0, delay=2*h.dt)
		self.olm = Population(cell_type=Ow,   n=int(math.ceil(200*scale)), x=20, y=0, z=0, dx=50, amp=-25e-3, dur=1e9, delay=2*h.dt)
		
		# psr = sensor cell to estimate the E->E connections
		self.psr = Population(cell_type=PyrAdr,n=1,   x= 0, y=0, z=0, dx=50, amp= 50e-3, dur=1e9, delay=2*h.dt) 
		self.cells = [self.pyr, self.bas, self.olm, self.psr]
		self.iseed = iseed # seed for noise inputs
		self.noise = noise
		self.DoMakeNoise = DoMakeNoise
		self.UseNetStim = UseNetStim
		self.wseed = wseed # seed for 'wiring'
		self.MSGain = MSGain # gain for MS weights
		self.RecPyr = False
		self.SaveConn = SaveConn
		
		if connections:
			print "Setting Connections"
			self.set_all_conns()


	def set_noise_inputs(self,simdur): #simdur only used for make_all_noise
		if self.DoMakeNoise:
			if self.UseNetStim:
				self.make_all_NetStims(simdur,self.iseed)
			else:
				self.make_all_noise(simdur,self.iseed)
		else:
			self.load_all_noise()
		print "Done!"

	def load_all_noise(self): #load noise from data files
		print "Loading Noise"		
		print "to PYR"
		self.b_pyr_somaAMPAf=self.load_spikes("spike_noise_pyr_800_soma_AMPA_ISI_1_N_10000_noise_1.npy",self.pyr,"somaAMPAf",0.05e-3)
		self.b_pyr_Adend3AMPAf=self.load_spikes("spike_noise_pyr_800_Adend3_AMPA_ISI_1_N_10000_noise_1.npy",self.pyr,"Adend3AMPAf",0.05e-3)
		self.b_pyr_somaGABAf=self.load_spikes("spike_noise_pyr_800_soma_GABA_ISI_1_N_10000_noise_1.npy",self.pyr,"somaGABAf",0.012e-3)
		self.b_pyr_Adend3GABAf=self.load_spikes("spike_noise_pyr_800_Adend3_GABA_ISI_1_N_10000_noise_1.npy",self.pyr,"Adend3GABAf",0.012e-3)
		self.b_pyr_Adend3NMDA=self.load_spikes("spike_noise_pyr_800_Adend3_NMDA_ISI_100_N_100_noise_1.npy", self.pyr,"Adend3NMDA",6.5e-3)
		print "to BAS"
		self.b_bas_somaAMPAf=self.load_spikes("spike_noise_bas_200_soma_AMPA_ISI_1_N_10000_noise_1.npy",self.bas,"somaAMPAf",w=0.02e-3)
		self.b_bas_somaGABA=self.load_spikes("spike_noise_bas_200_soma_GABA_ISI_1_N_10000_noise_1.npy",self.bas,"somaGABAf",w=0.2e-3)
		self.b_bas_somaGABAf=self.load_spikes("spike_noise_bas_200_soma_GABAf_ISI_150_N_65_noise_0.npy",self.bas,"somaGABAss",w=1.6e-3)
		print "to OLM"
		self.b_olm_somaAMPAf=self.load_spikes("spike_noise_olm_200_soma_AMPA_ISI_1_N_10000_noise_1.npy",self.olm,"somaAMPAf",w=0.02e-3)
		self.b_olm_somaGABAf=self.load_spikes("spike_noise_olm_200_soma_GABA_ISI_1_N_10000_noise_1.npy",self.olm,"somaGABAf",w=0.2e-3)
		self.b_olm_somaGABAss=self.load_spikes("spike_noise_olm_200_soma_GABAss_ISI_150_N_65_noise_0.npy",self.olm,"somaGABAss",w=1.6e-3)

	#this should be called @ beginning of each sim - done in an FInitializeHandler in run.py
	def init_NetStims(self):
		# h.mcell_ran4_init(self.iseed)
		for i in xrange(len(self.nrl)):
			rds = self.nrl[i]
			sead = self.nrlsead[i]
			rds.MCellRan4(sead,sead)
			rds.negexp(1)			
			# print i,rds,sead

	#creates NetStims (and associated NetCon,Random) - provide 'noise' inputs
	#returns next useable value of sead
	def make_NetStims(self,po,syn,w,ISI,time_limit,sead):
		po.nssidx[syn] = len(self.nsl) #index into net.nsl
		po.ncsidx[syn] = len(self.ncl) #index into net.ncl
		for i in range(po.n):
			cel = po.cell[i]

			ns = h.NetStim()
			ns.interval = ISI
			ns.noise = 1			
			ns.number = (1e3 / ISI) * time_limit
			ns.start = 0

			nc = h.NetCon(ns,cel.__dict__[syn].syn)
			nc.delay = h.dt * 2 # 0
			nc.weight[0] = w

			rds = h.Random()
			rds.negexp(1)            # set random # generator using negexp(1) - avg interval in NetStim
			rds.MCellRan4(sead,sead) # seeds are in order, shouldn't matter			
			ns.noiseFromRandom(rds)  # use random # generator for this NetStim
			
			#ns.start = rds.discunif(0,1e3) # start inputs random time btwn 0-1e3 ms to avoid artificial sync
			#rds.MCellRan4(sead,sead) # reinit rand # generator

			self.nsl.append(ns)
			self.ncl.append(nc)
			self.nrl.append(rds)
			self.nrlsead.append(sead)
			sead = sead + 1

		po.nseidx[syn] = len(self.nsl)-1
		po.nceidx[syn] = len(self.ncl)-1
		
		return sead

	# setup recording of pyramidal cell inputs, assumes using NetCon,NetStims
	def RecPYRInputs(self):
		self.RecPyr = True
		self.NCV = {}
		self.sys = ['somaAMPAf', 'Adend3AMPAf', 'somaGABAf', 'Adend3GABAf']
		sys=self.sys
		for s in sys:
			self.NCV[s] = []
			sidx = self.pyr.ncsidx[s]
			eidx = self.pyr.nceidx[s]
			for i in xrange(sidx,eidx+1):
				self.NCV[s].append(h.Vector())
				self.ncl[i].record(self.NCV[s][-1])

	# make an NQS with pyramidal cell input times
	def setnqin(self):
		try:
			h.nqsdel(self.nqin)
		except:
			pass
		self.nqin = h.NQS("id","sy","vt")
		nqin=self.nqin
		nqin.odec("vt")
		jdx = 0
		for s in self.sys:
			sidx = self.pyr.ncsidx[s]
			eidx = self.pyr.nceidx[s]
			idx = 0
			for i in xrange(0,len(self.NCV[s])):
				nqin.append(idx,jdx,self.NCV[s][i])
				idx = idx + 1
			jdx = jdx + 1
				
	# make a histogram of pyramidal cell spike outputs
	def mkspkh(self,binsz):
		snq=self.snq
		snq.verbose = 0
		self.spkh = h.List()
		for i in xrange(0,800):
			if snq.select("id",i) > 0:
				vt = snq.getcol("t")
				self.spkh.append(vt.histogram(0,h.tstop,binsz))
			else:
				self.spkh.append(h.Vector())
		snq.verbose=1

	def make_all_NetStims(self,simdur,rdmseed):
		print "Making NetStims"
		# h.mcell_ran4_init(self.iseed)
		self.nsl = [] #NetStim List
		self.ncl = [] #NetCon List
		self.nrl = [] #Random List for NetStims
		self.nrlsead = [] #List of seeds for NetStim randoms
		# numpy.random.seed(rdmseed) # initialize random # generator
		print "Making Noise"
		print "to PYR"
		rdtmp = rdmseed # starting sead value - incremented in make_NetStims
		rdtmp=self.make_NetStims(po=self.pyr, syn="somaAMPAf",   w=0.05e-3,  ISI=1,  time_limit=simdur, sead=rdtmp) 
		rdtmp=self.make_NetStims(po=self.pyr, syn="Adend3AMPAf", w=0.05e-3,  ISI=1,  time_limit=simdur, sead=rdtmp)
		rdtmp=self.make_NetStims(po=self.pyr, syn="somaGABAf",   w=0.012e-3, ISI=1,  time_limit=simdur, sead=rdtmp)
		rdtmp=self.make_NetStims(po=self.pyr, syn="Adend3GABAf", w=0.012e-3, ISI=1,  time_limit=simdur, sead=rdtmp)
		rdtmp=self.make_NetStims(po=self.pyr, syn="Adend3NMDA",  w=6.5e-3,   ISI=100,time_limit=simdur, sead=rdtmp)
		print "to BAS"			
		rdtmp=self.make_NetStims(po=self.bas, syn="somaAMPAf",   w=0.02e-3,  ISI=1,  time_limit=simdur, sead=rdtmp)
		rdtmp=self.make_NetStims(po=self.bas, syn="somaGABAf",   w=0.2e-3,   ISI=1,  time_limit=simdur, sead=rdtmp)
		print "to OLM"
		#rdtmp=self.make_NetStims(po=self.olm, syn="somaAMPAf",   w=0.02e-3,  ISI=1,  time_limit=simdur, sead=rdtmp)
		rdtmp=self.make_NetStims(po=self.olm, syn="somaAMPAf",   w=0.0625e-3,  ISI=1,  time_limit=simdur, sead=rdtmp)
		rdtmp=self.make_NetStims(po=self.olm, syn="somaGABAf",   w=0.2e-3,   ISI=1,  time_limit=simdur, sead=rdtmp)
		#setup medial septal inputs to OLM and BASKET cells, note that MSGain can be 0 == no effect
		ns = h.NetStim()
		ns.interval = 150
		ns.noise = 0 # NO randomness for the MS inputs
		ns.number = (1e3 / 150.0) * simdur
		self.nsl.append(ns)
		for i in range(self.bas.n): # MS inputs to BASKET cells
			nc = h.NetCon(ns,self.bas.cell[i].__dict__["somaGABAss"].syn)
			nc.delay = 2*h.dt
			nc.weight[0] = 1.6e-3 * self.MSGain
			self.ncl.append(nc)
		for i in range(self.olm.n): # MS inputs to OLM cells
			nc = h.NetCon(ns,self.olm.cell[i].__dict__["somaGABAss"].syn)
			nc.delay = 2*h.dt
			nc.weight[0] = 1.6e-3 * self.MSGain
			self.ncl.append(nc)
		
	def make_all_noise(self,simdur,rdmseed): # create noise for simdur milliseconds
		numpy.random.seed(rdmseed) # initialize random # generator
		import math
		print "Making Noise"
		fctr = (simdur+simdur/2) / 10000.0		
		print "to PYR"
		self.b_pyr_somaAMPAf=self.make_spikes(self.pyr,"somaAMPAf",0.05e-3,self.pyr.n,"soma",1,math.ceil(10000*fctr),1,simdur)
		self.b_pyr_Adend3AMPAf=self.make_spikes(self.pyr,"Adend3AMPAf",0.05e-3,self.pyr.n,"Adend3",1,math.ceil(10000*fctr),1,simdur)
		self.b_pyr_somaGABAf=self.make_spikes(self.pyr,"somaGABAf",0.012e-3,self.pyr.n,"soma",1,math.ceil(10000*fctr),1,simdur)
		self.b_pyr_Adend3GABAf=self.make_spikes(self.pyr,"Adend3GABAf",0.012e-3,self.pyr.n,"Adend3",1,math.ceil(10000*fctr),1,simdur)
		self.b_pyr_Adend3NMDA=self.make_spikes(self.pyr,"Adend3NMDA",6.5e-3,self.pyr.n,"Adend3",100,math.ceil(100*fctr),1,simdur)
		print "to BAS"			
		self.b_bas_somaAMPAf=self.make_spikes(self.bas,"somaAMPAf",0.02e-3,self.bas.n,"soma",1,math.ceil(10000*fctr),1,simdur)
		self.b_bas_somaGABA=self.make_spikes(self.bas,"somaGABAf",0.2e-3,self.bas.n,"soma",1,math.ceil(10000*fctr),1,simdur)
		self.b_bas_somaGABAf=self.make_spikes(self.bas,"somaGABAss",1.6e-3,self.bas.n,"soma",150,math.ceil(65*fctr),0,simdur)
		print "to OLM"
		self.b_olm_somaAMPAf=self.make_spikes(self.olm,"somaAMPAf",0.02e-3,self.olm.n,"soma",1,math.ceil(10000*fctr),1,simdur)
		self.b_olm_somaGABAf=self.make_spikes(self.olm,"somaGABAf",0.2e-3,self.olm.n,"soma",1,math.ceil(10000*fctr),1,simdur)	
		self.b_olm_somaGABAss=self.make_spikes(self.olm,"somaGABAss",1.6e-3,self.olm.n,"soma",150,math.ceil(65*fctr),0,simdur)		

	def make_conn(self, preN, postN, conv):
		conn = numpy.zeros((postN,conv),dtype=numpy.int16)
		for i in range(postN):
			conn[i,:]=random.sample(range(preN),conv)
		return conn

	def set_all_conns(self):
		random.seed(self.wseed) # initialize random # generator for wiring
		print "PYR -> X , NMDA"   # src, trg, syn, delay, weight, conv
		self.pyr_bas_NM=self.set_connections(self.pyr,self.bas, "somaNMDA", 2, 1.15*1.2e-3, 100)
		self.pyr_olm_NM=self.set_connections(self.pyr,self.olm, "somaNMDA", 2, 1.0*0.7e-3, 10)
		self.pyr_pyr_NM=self.set_connections(self.pyr,self.pyr, "BdendNMDA",2, 1*0.004e-3,  25)

		print "PYR -> X , AMPA"
		self.pyr_bas_AM=self.set_connections(self.pyr,self.bas, "somaAMPAf",2, 0.3*1.2e-3,  100)
		self.pyr_olm_AM=self.set_connections(self.pyr,self.olm, "somaAMPAf",2, 0.3*1.2e-3,  10)
		self.pyr_pyr_AM=self.set_connections(self.pyr,self.pyr, "BdendAMPA",2, 0.5*0.04e-3, 25)
			
		print "BAS -> X , GABA"
		#self.bas_bas_GA=self.set_connections(self.bas,self.bas, "somaGABAf",2, 1.0e-3, 60)#orig 1
		#self.bas_bas_GA=self.set_connections(self.bas,self.bas, "somaGABAf",2, 2  *  1.5*1.0e-3, 60)#new 1
		self.bas_bas_GA=self.set_connections(self.bas,self.bas, "somaGABAf",2, 3  *  1.5*1.0e-3, 60)#new 2
		self.bas_pyr_GA=self.set_connections(self.bas,self.pyr, "somaGABAf",2, 2  *  2*0.18e-3, 50)#new 1

		print "OLM -> PYR , GABA"
		#self.olm_pyr_GA=self.set_connections(self.olm,self.pyr, "Adend2GABAs",2, 3*6.0e-3, 20)#original weight value
		self.olm_pyr_GA=self.set_connections(self.olm,self.pyr, "Adend2GABAs",2, 4.0  *  3*6.0e-3, 20)#new weight value

	        #pyramidal to PSR cell -- for testing only
		print "PYR -> PSR, AMPA/NMDA"
		self.pyr_psr_NM=self.set_connections(self.pyr,self.psr, "BdendNMDA",2, 1*0.004e-3,  25)
		self.pyr_psr_AM=self.set_connections(self.pyr,self.psr, "BdendAMPA",2, 0.5*0.04e-3, 25)


	def set_conn_weight(self, conn, weight):
		for nc in conn:
			nc.weight[0] = weight
			
	def set_connections(self,src,trg,syn,delay,w,conv):
		conn = self.make_conn(src.n,trg.n,conv)
		nc = []
		for post_id, all_pre in enumerate(conn):
			for j, pre_id in enumerate(all_pre):
				nc.append(h.NetCon(src.cell[pre_id].soma(0.5)._ref_v, trg.cell[post_id].__dict__[syn].syn, 0, delay, w, sec=src.cell[pre_id].soma))	
		if self.SaveConn:
			try:
				print self.nqcon.size()
			except:
				self.nqcon = h.NQS("id1","id2","w","syn")
				self.nqcon.strdec("syn")
			for post_id, all_pre in enumerate(conn):
				for j, pre_id in enumerate(all_pre):
					self.nqcon.append(src.cell[pre_id].id,trg.cell[post_id].id,w,syn)	
			
		return nc

	def load_spikes(self,fn,po,syn,w,time_limit=10000):
		fn = os.path.join("data",fn)
		events = numpy.load(fn)
		print "Begin setting events...", po
		print events.shape
		for i,ii in enumerate(events):
			ii=ii[ii<=time_limit]
			po.cell[i].__dict__[syn].append(ii)
			po.cell[i].__dict__[syn].syn.Vwt = w
		print "End setting events"
		return events

	def make_spikes(self,po,syn,w,cellN,comp,ISI,eventN,noise,time_limit):
		events = numpy.random.exponential(ISI, (cellN,eventN))*noise+numpy.repeat(ISI,cellN*eventN).reshape((cellN,eventN))*(1-noise)
		events = numpy.cumsum(events,axis=1)
		print "Begin setting events...", po
		print events.shape
		for i,ii in enumerate(events):
			ii=ii[ii<=time_limit]
			po.cell[i].__dict__[syn].append(ii)
			po.cell[i].__dict__[syn].syn.Vwt = w
		print "End setting events"
		return events		
			
	def rasterplot(self,sz=2):
		pon  = 0		
		if h.g[0] == None:
			h.gg()
		col = [2, 4, 3, 1]
		for po in self.cells:
			id = h.Vector()
			tv = h.Vector()
			for i in xrange(po.n):
				id.append(po.lidvec[i])
				tv.append(po.ltimevec[i])
			id.mark(h.g[0],tv,"O",sz,col[pon],1)
			pon += 1
		h.g[0].exec_menu("View = plot")

	def setrastervecs(self):
		self.myidvec = h.Vector() #IDs and firing times for ALL cells
		self.mytimevec = h.Vector()
		for po in self.cells:
			for i in xrange(po.n):
				self.myidvec.append(po.lidvec[i])
				self.mytimevec.append(po.ltimevec[i])

	# setsnq - make an NQS with ids, spike times, types
	def setsnq(self):
		try:
			h.nqsdel(self.snq)
		except:
			pass
		self.snq = h.NQS("id","t","ty")
		ty = 0
		vec = h.Vector()
		for po in self.cells:
			for i in xrange(po.n):
				self.snq.v[0].append(po.lidvec[i])
				self.snq.v[1].append(po.ltimevec[i])
				vec.resize(po.lidvec[i].size())
				vec.fill(ty)
				self.snq.v[2].append(vec)
			ty += 1


	# setfnq - make an NQS with ids, firing rates, types
	def setfnq(self,skipms=200):
		try: 
			self.snq.tog("DB")
		except:
			self.setsnq()
		try:
			h.nqsdel(self.fnq)
		except:
			pass
		self.fnq = h.NQS("id","freq","ty")
		tf = h.tstop - skipms
		ty = 0
		for po in self.cells:
			for i in xrange(po.n):
				id = po.cell[i].id
				n = float( self.snq.select("t",">",skipms,"id",id) )
				self.fnq.append(id, n*1e3/tf, ty)
			ty += 1

	# pravgrates - print average firing rates using self.fnq
	def pravgrates(self,skipms=200):
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
				print "ty: ", ty, " avg rate = ", vf.mean(), "+/-", vf.stderr(), " Hz"
			else:
				print "ty: ", ty, " avg rate = ", vf.mean(), "+/-", 0.0 , " Hz"
			ty += 1
			
	def calc_lfp(self): # lfp is modeled as a difference between voltages in distal apical and basal compartemnts 
		self.vlfp = h.Vector(self.pyr.cell[0].Adend3_volt.size()) #lfp in neuron Vector
		for cell in self.pyr.cell: 
			self.vlfp.add(cell.Adend3_volt)
			self.vlfp.sub(cell.Bdend_volt)
		self.vlfp.div(len(self.pyr.cell)) # normalize lfp by amount of pyr cells
		self.lfp=numpy.array(self.vlfp.to_python()) # convert to python array (so can do PSD)

	def calc_specgram(self,maxfreq,nsamp,dodraw,skipms=0):
		self.calc_lfp()
		if skipms > 0:
			vtmp = h.Vector()
			vtmp.copy(self.vlfp,skipms/h.dt,self.vlfp.size()-1)
			self.MSpec = MSpec(vtmp,maxfreq,nsamp,dodraw)
		else:
			self.MSpec = MSpec(self.vlfp,maxfreq,nsamp,dodraw)
		
	def calc_psd(self,fig=3):
		self.calc_lfp()
		t0   = 200 # reject first ms of the signal
		fmax = 200 # upper limit for a periodogram frequency
		div  = int(1000/h.dt/(2*fmax)) # downsample the signal
		tr = [3,  12] # Theta frequency range
		gr = [30, 80] # Gamma frequency range
		t0i = int(t0/h.dt)
		if t0i > len(self.lfp):
			print "LFP is too short! (<200 ms)"
			return 0,0,0,0,0,0
		
		pylab.figure(fig)
		pylab.clf()
		
		pylab.subplot(2,1,1) # plot LFP
		pylab.plot(numpy.array(range(len(self.lfp)))*h.dt, self.lfp)
		
		pylab.subplot(2,1,2) # plot periodogram
		data = self.lfp[t0i::div] # downsample data
		Pxx, freqs = pylab.psd(data-data.mean(), Fs=1000/h.dt/div) # calculate FFT
		tind = numpy.where((freqs>=tr[0]) & (freqs<=tr[1]))[0] # index where for theta frequences  
		gind = numpy.where((freqs>=gr[0]) & (freqs<=gr[1]))[0] # index where for gamma frequences
		self.tp = Pxx[tind].mean() * numpy.diff(tr) # integral over theta power
		self.gp = Pxx[gind].mean() * numpy.diff(gr) # integral over gamma power
		self.ti = self.get_lim_max(Pxx, tind) # index of the frequency with a maximal power in theta range  
		self.gi = self.get_lim_max(Pxx, gind) # index of the frequency with a maximal power in gamma range
		self.tf = freqs[self.ti]
		self.gf = freqs[self.gi]
		pylab.scatter(self.tf, 10*numpy.log10(Pxx[self.ti]), 100, 'b','o')
		pylab.scatter(self.gf, 10*numpy.log10(Pxx[self.gi]), 100, 'r','o')
		pylab.xlim(0,fmax)
			
	def get_lim_max(self, data, ind): 	# return the position of the maximal element in data located in the postion indexed by ind
		return  ind[data[ind].argmax()]


#make the Network - use params in rseed.txt if the file exists -- makes it easier to run a batch
#if rseed.txt doesn't exist, the Network is created with default params
try:
	fp = open("./rseed.txt","r")
	ls = fp.readlines()
	ISEED = int(ls[0])
	WSEED = int(ls[1])
	MSG = 1.0
	if len(ls) > 2:
		MSG = float(ls[2])
	fp.close()
        #create the network
	net = Network(noise=True,connections=True,DoMakeNoise=True,iseed=ISEED,UseNetStim=True,wseed=WSEED,scale=1.0,MSGain=MSG) 
	print "set network from rseed.txt : iseed=",ISEED,", WSEED=",WSEED,", MSG = ",MSG
except:
	net = Network()
	print "set network from default constructor"

#setup some variables in hoc
def sethocix():
	h("PYRt=0")
	h("BASKETt=1")
	h("OLMt=2")
	h("PSRt=3")
	h("CTYP.o(PYRt).s=\"PYRt\"")
	h("CTYP.o(BASKETt).s=\"BASKETt\"")
	h("CTYP.o(OLMt).s=\"OLMt\"")
	h("CTYP.o(PSRt).s=\"PSRt\"")
	h("ix[PYRt]=0")
	h("ixe[PYRt]=799")
	h("ix[BASKETt]=800")
	h("ixe[BASKETt]=999")
	h("ix[OLMt]=1000")
	h("ixe[OLMt]=1199")
	h("ix[PSRt]=1200")
	h("ixe[PSRt]=1200")
	h("numc[PYRt]=800")
	h("numc[BASKETt]=200")
	h("numc[OLMt]=200")
	h("numc[PSRt]=1")

sethocix()
