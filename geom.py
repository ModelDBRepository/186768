# $Id: geom.py,v 1.30 2011/01/02 19:56:14 samn Exp $ 

from pyinit import *

class Synapse:
	def __init__(self, sect, loc, tau1, tau2, e):
		self.syn		= h.MyExp2SynBB(loc, sec=sect)
		self.syn.tau1	= tau1
		self.syn.tau2	= tau2
		self.syn.e		= e 
		
class SynapseNMDA:
	def __init__(self, sect, loc, tau1, tau2, tau1NMDA, tau2NMDA, r, e):
		self.syn			= h.MyExp2SynNMDABB(loc, sec=sect)
		self.syn.tau1		= tau1
		self.syn.tau2		= tau2
		self.syn.tau1NMDA	= tau1NMDA
		self.syn.tau2NMDA	= tau2NMDA 
		self.syn.r			= r
		self.syn.e			= e 
		
###############################################################################
#
# General Cell
#
###############################################################################
class Cell:
	"General cell"
	
	def __init__(self,x,y,z,id):
		self.x=x
		self.y=y
		self.z=z
		self.id=id
		self.all_sec = []
		self.add_comp('soma',True)
		self.set_morphology()
		self.set_conductances()
		self.set_synapses()
		self.set_inj()
		self.calc_area()
		
	def set_morphology(self):
		pass
			
	def set_conductances(self):
		pass
		
	def set_synapses(self):
		pass
		
	def set_inj(self):
		self.somaInj = h.IClamp(0.5, sec=self.soma)	
		
	def add_comp(self, name, rec):
		self.__dict__[name] = h.Section()
		self.all_sec.append(self.__dict__[name])
		# Record voltage
		if rec:
			self.__dict__[name+"_volt"] = h.Vector(int(h.tstop/h.dt)+1)
			self.__dict__[name+"_volt"].record(self.__dict__[name](0.5)._ref_v)
	
	def plot_volt(self, name, fig=1):
		figure(fig)
		volt = self.__dict__[name+"_volt"].to_python()
		plot(arange(len(volt))*h.dt, volt)
		
	def calc_area(self):
		self.total_area = 0
		self.n = 0
		for sect in self.all_sec:
			self.total_area += h.area(0.5,sec=sect)
			self.n+=1
			
###############################################################################
#
# Basket Cell -- Bwb
#
###############################################################################

class Bwb(Cell):
	"Basket cell"
	
	def set_morphology(self):
		total_area = 10000 # um2
		self.soma.nseg  = 1
		self.soma.cm    = 1      # uF/cm2
		diam = sqrt(total_area) # um
		L    = diam/pi  # um
			
		h.pt3dclear(sec=self.soma)
		h.pt3dadd(self.x, self.y, self.z,   diam, sec=self.soma)
		h.pt3dadd(self.x, self.y, self.z+L, diam, sec=self.soma)
			
	def set_conductances(self):
		self.soma.insert('pas')
		self.soma.e_pas = -65     # mV
		self.soma.g_pas = 0.1e-3  # S/cm2 
	  
		self.soma.insert('Nafbwb')
		self.soma.insert('Kdrbwb')
	   
	def set_synapses(self):
		self.somaAMPAf 	= Synapse(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, e=0)
		self.somaGABAf 	= Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80)
		self.somaGABAss	= Synapse(sect=self.soma, loc=0.5, tau1=20,   tau2=40, e=-80)#only for septal input
		self.somaNMDA 	= SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15, tau2NMDA=150, r=1, e=0)
		
###############################################################################
#
# OLM Cell -- Ow
#
###############################################################################
class Ow(Cell):
	"OLM cell"
   
	def set_morphology(self):
		total_area = 10000 # um2
		self.soma.nseg  = 1
		self.soma.cm    = 1      # uF/cm2
		diam = sqrt(total_area) # um
		L    = diam/pi  # um

		h.pt3dclear(sec=self.soma)
		h.pt3dadd(self.x, self.y, self.z,   diam, sec=self.soma)
		h.pt3dadd(self.x, self.y, self.z+L, diam, sec=self.soma)
	
	def set_conductances(self):
		self.soma.insert('pas')
		self.soma.e_pas = -65     # mV
		self.soma.g_pas = 0.1e-3  # S/cm2 

		self.soma.insert('Nafbwb')
		self.soma.insert('Kdrbwb')

		self.soma.insert('Iholmw')
		self.soma.insert('Caolmw')
		self.soma.insert('ICaolmw')
		self.soma.insert('KCaolmw')

	def set_synapses(self):
		self.somaGABAf 	= Synapse(sect=self.soma, loc=0.5, tau1=0.07, tau2=9.1, e=-80)
		self.somaAMPAf 	= Synapse(    sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, e=0)
		self.somaGABAss	= Synapse(    sect=self.soma, loc=0.5, tau1=20,	  tau2=40, e=-80)#only for septal input
		self.somaNMDA 	= SynapseNMDA(sect=self.soma, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15, tau2NMDA=150, r=1, e=0)
		
###############################################################################
#
# Pyramidal Cell -- KopAdr
#
###############################################################################
class PyrAdr(Cell):
	"Pyramidal cell"

	def set_morphology(self):
		self.add_comp('Bdend',True)
		self.add_comp('Adend1',False)
		self.add_comp('Adend2',False)
		self.add_comp('Adend3',True)

		h.pt3dclear(sec=self.soma)
		h.pt3dadd(self.x, self.y, self.z,          20, sec=self.soma)
		h.pt3dadd(self.x, self.y, self.z+20,       20, sec=self.soma)

		h.pt3dclear(sec=self.Bdend)
		h.pt3dadd(self.x, self.y, self.z,          2, sec=self.Bdend)
		h.pt3dadd(self.x, self.y, self.z-200,      2, sec=self.Bdend)

		h.pt3dclear(sec=self.Adend1)
		h.pt3dadd(self.x, self.y, self.z+20,       2, sec=self.Adend1)
		h.pt3dadd(self.x, self.y, self.z+20+150,   2, sec=self.Adend1)

		h.pt3dclear(sec=self.Adend2)
		h.pt3dadd(self.x, self.y, self.z+20+150,   2, sec=self.Adend2)
		h.pt3dadd(self.x, self.y, self.z+20+150*2, 2, sec=self.Adend2)

		h.pt3dclear(sec=self.Adend3)
		h.pt3dadd(self.x, self.y, self.z+20+150*2, 2, sec=self.Adend3)
		h.pt3dadd(self.x, self.y, self.z+20+150*3, 2, sec=self.Adend3)

		self.Bdend.connect(self.soma,      0, 0)
		self.Adend1.connect(self.soma,   0.5, 0)
		self.Adend2.connect(self.Adend1,   1, 0)
		self.Adend3.connect(self.Adend2,   1, 0)

	def set_conductances(self):
		for sect in self.all_sec:
			sect.insert('pas')
			sect(0.5).pas.g = 0.0000357
			sect.insert('nacurrent')
			sect.insert('kacurrent')
			sect.insert('kdrcurrent')
			sect.insert('hcurrent')
			sect(0.5).pas.e = -70     # mV
			sect.cm = 1
			sect.Ra = 150

		self.Adend1(0.5).nacurrent.ki = 0.5
		self.Adend1(0.5).kacurrent.g  = 0.072
		self.Adend1(0.5).hcurrent.v50 = -82
		self.Adend1(0.5).hcurrent.g   = 0.0002
		
		self.Adend2(0.5).nacurrent.ki = 0.5
		self.Adend2(0.5).kacurrent.g  = 0
		self.Adend2(0.5).kacurrent.gd = 0.120
		self.Adend2(0.5).hcurrent.v50 = -90
		self.Adend2(0.5).hcurrent.g   = 0.0004
		
		self.Adend3(0.5).cm           = 2
		self.Adend3(0.5).pas.g        = 0.0000714
		self.Adend3(0.5).nacurrent.ki = 0.5
		self.Adend3(0.5).kacurrent.g  = 0
		self.Adend3(0.5).kacurrent.gd = 0.200		
		self.Adend3(0.5).hcurrent.v50 = -90
		self.Adend3(0.5).hcurrent.g   = 0.0007
		
		self.Bdend(0.5).nacurrent.ki  = 1

	def set_synapses(self):
		self.somaGABAf 	 = Synapse(    sect=self.soma,   loc=0.5, tau1=0.07, tau2=9.1, 	  e=-80)
		self.somaAMPAf 	 = Synapse(    sect=self.soma,   loc=0.5, tau1=0.05, tau2=5.3, 	   e=0)
		self.BdendAMPA   = Synapse(    sect=self.Bdend,  loc=1.0, tau1=0.05, tau2=5.3,     e=0)
		self.BdendNMDA   = SynapseNMDA(sect=self.Bdend,  loc=1.0, tau1=0.05, tau2=5.3, tau1NMDA=15, tau2NMDA=150, r=1, e=0)
		self.Adend2GABAs = Synapse(	   sect=self.Adend2, loc=0.5, tau1=0.2,  tau2=20,   e=-80)
		self.Adend3GABAf = Synapse(	   sect=self.Adend3, loc=0.5, tau1=0.07, tau2=9.1,   e=-80)
		self.Adend3AMPAf = Synapse(	   sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3,    e=0)
		self.Adend3NMDA  = SynapseNMDA(sect=self.Adend3, loc=0.5, tau1=0.05, tau2=5.3, tau1NMDA=15, tau2NMDA=150, r=1, e=0)

