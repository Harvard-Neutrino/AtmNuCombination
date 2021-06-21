import numpy as np

class Event:
	def __init__(self,pdg_lep,e_lep,p_lep_v,p_lep_norm,n_had_fp,pdg_fp,e_fp,dir_fp,p_fp_v,p_fp_norm):
		self.NumberOfHadrons = n_had_fp
		self.HadronPDG = pdg_fp
		self.HadronEnergy = e_fp
		self.HadronDirection = dir_fp
		self.HadronMomentum = p_fp_norm
		self.HadronMomenta = p_fp_v

		self.LeptonPDG = pdg_lep
		self.LeptonEnergy = e_lep
		self.LeptonDirection = p_lep_v / p_lep_norm
		self.LeptonMomentum = p_lep_norm
		self.LeptonMomenta = p_lep_v

		self.NumberOfParticles = 1+n_had_fp
		self.ParticlePDG = np.insert(pdg_fp, 0, pdg_lep)
		self.ParticleEnergy = np.insert(e_fp, 0, e_lep)
		self.ParticleMomentum = np.insert(p_fp_norm, 0, p_lep_norm)
		self.ParticleDirection = np.vstack((self.LeptonDirection, self.HadronDirection))
		self.ParticleMomenta = np.vstack((self.LeptonMomenta, self.HadronMomenta))