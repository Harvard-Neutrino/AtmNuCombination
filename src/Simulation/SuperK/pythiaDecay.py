import pythia8
import numpy as np
from math import sqrt

class pythiaDecay:
	def __init__(self):
		self.UnstableParents = {15:"tau-",-15:"tau+",
			  321:"k+",-321:"k-",211:"pi+",111:"pi0",
			  221:"eta",411:"d+",-411:"d-",421:"d0"}
		self.pythia = pythia8.Pythia()
		self.pythia.readString("ProcessLevel:all = off")
		self.pythia.readString("Random:setSeed = on")
		self.pythia.readString("Random:seed = 0")
		self.pythia.init()

		for part in self.UnstableParents:
			self.pythia.particleData.mayDecay(part,True)

	def decay(self, pdgid, E, p):
		pid = int(pdgid)
		mass = self.pythia.particleData.m0(pid)
		mom = np.linalg.norm(p)
		E = sqrt(mom**2 + mass**2) # re-compute energy as some mass values differ a little
		p4vec = pythia8.Vec4(p[0],p[1],p[2],E)

		self.pythia.event.reset()
		self.pythia.event.append(pid,91,0,0,p4vec,mass)
		self.pythia.forceHadronLevel()

		final_p = np.array([])
		final_pv = np.array([])
		final_E = np.array([])
		final_pid = np.array([])

		i = -1
		for k in range(self.pythia.event.size()):
			if self.pythia.event[k].isFinal() and self.pythia.event[k].tau()/300<250: #decay time cut
				i += 1
				dummyp = np.array([self.pythia.event[k].p()[0], self.pythia.event[k].p()[1], self.pythia.event[k].p()[2]])
				dummymom = np.linalg.norm(dummyp)
				final_E = np.append(final_E,self.pythia.event[k].e())
				final_p = np.append(final_p,dummymom)
				final_pid = np.append(final_pid, self.pythia.event[k].id())
				if i==0:
					final_pv = dummyp
				else:
					final_pv = np.vstack((final_pv,dummyp))

		return final_pid, final_E, final_p, final_pv


	def canDecay(self, pdgid):
		if pdgid not in self.UnstableParents:
			return False
		else:
			return True
