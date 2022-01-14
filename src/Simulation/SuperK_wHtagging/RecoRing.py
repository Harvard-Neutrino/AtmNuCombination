import numpy as np
from particle import Particle
import ParticleProperties as pp
import math
from math import log10, pi
import random
import applications as ap

class RecoRing:
	AngResoThr = 0.975 #math.cos(25.*pi/180.)
	MeVtoGeV = 0.001
	def __init__(self, TrueNRing, TrueRingPDG, TrueRingIP, TrueRingE, TrueRingP, TrueRingDir, distros, mode):
		self.distros = distros
		self.NRing = TrueNRing
		self.CC = True if abs(mode)<30 else False
		self.RecoLabels = np.ones(TrueNRing)
		self.HeavyChargedWA(TrueRingPDG, TrueRingP)
		# Ring IP reconstruction
		self.ReconSRIP(TrueRingIP, TrueRingP)
		# Momentum reconstruction of each ring
		self.ReconMomentum(TrueRingPDG, TrueRingP)
		# Direction reconstruction of each ring
		self.ReconDirection(TrueRingDir)
		if self.NRing >1: # Oh, the multi-ringers! Special care needed
			dmer = np.where(self.Momentum == np.amax(self.Momentum))
			mer = dmer[0][0]
			self.ReconMRIP(mer, TrueRingIP[mer], TrueRingPDG[mer])
			mTrueRingPDG, mTrueRingP, mTrueRingDir = self.TooClose(TrueRingPDG, TrueRingP, TrueRingDir) # Put labels on rings which won't be reconstructed due to poor angular resolution
			excessRing = self.NRing - 5
			if excessRing > 0:
				for index in np.argsort(mTrueRingP)[:excessRing]:
					self.RecoLabels[index] = 0
			if np.any(self.RecoLabels==0):
				self.ReconMomentum(mTrueRingPDG, mTrueRingP)
				self.ReconDirection(mTrueRingDir)

		self.TotalVariables()
		self.MuEdk = 0
		self.Neutrons = 0
		self.Imass = 0

	def HeavyChargedWA(self, pdg, p):
		for i,part in enumerate(pdg):
			if abs(part)==211 and p[i]<0.3:
				self.RecoLabels[i] = 0
			elif abs(part)==2212 and p[i]<2.3:
				self.RecoLabels[i] = 0

	def SKType(self, ipnu, cc, ntag):
		
		if self.NRing>1 and self.MERIP==2:
			if cc:
				ip_r = self.distros.Random('other_cc')
				if ipnu>0:
					ip_n = self.distros.Random('mre_nunubar_ccnue')
				else:
					ip_n = self.distros.Random('mre_nunubar_ccnuebar')
			else:
				ip_r = self.distros.Random('other_nc')
				ip_n = self.distros.Random('mre_nunubar_nc')

		if self.NRing==1 and self.MERIP==2:
			pi0m_r = 0.
			pi0d_r = 0.
			pi0s_r = 0.
			if self.Evis<0.25:
				pi0m_cc, wm_cc = self.distros.Random_and_Weight('sgpi0_cc_0')
				pi0m_nc, wm_nc = self.distros.Random_and_Weight('sgpi0_nc_0')
			elif self.Evis<0.4:
				pi0m_cc, wm_cc = self.distros.Random_and_Weight('sgpi0_cc_1')
				pi0s_cc, ws_cc = self.distros.Random_and_Weight('pi0_signal_0_cc')
				pi0d_cc, wd_cc = self.distros.Random_and_Weight('pi0_delta_0_cc')
				pi0m_nc, wm_nc = self.distros.Random_and_Weight('sgpi0_nc_1')
				pi0s_nc, ws_nc = self.distros.Random_and_Weight('pi0_signal_0_nc')
				pi0d_nc, wd_nc = self.distros.Random_and_Weight('pi0_delta_0_nc')
			elif self.Evis<0.63:
				pi0m_cc, wm_cc = self.distros.Random_and_Weight('sgpi0_cc_2')
				pi0s_cc, ws_cc = self.distros.Random_and_Weight('pi0_signal_1_cc')
				pi0d_cc, wd_cc = self.distros.Random_and_Weight('pi0_delta_1_cc')
				pi0m_nc, wm_nc = self.distros.Random_and_Weight('sgpi0_nc_2')
				pi0s_nc, ws_nc = self.distros.Random_and_Weight('pi0_signal_1_nc')
				pi0d_nc, wd_nc = self.distros.Random_and_Weight('pi0_delta_1_nc')
			elif self.Evis<1:
				pi0m_cc, wm_cc = self.distros.Random_and_Weight('sgpi0_cc_3')
				pi0s_cc, ws_cc = self.distros.Random_and_Weight('pi0_signal_2_cc')
				pi0d_cc, wd_cc = self.distros.Random_and_Weight('pi0_delta_2_cc')
				pi0m_nc, wm_nc = self.distros.Random_and_Weight('sgpi0_nc_3')
				pi0s_nc, ws_nc = self.distros.Random_and_Weight('pi0_signal_2_nc')
				pi0d_nc, wd_nc = self.distros.Random_and_Weight('pi0_delta_2_nc')
			else:
				pi0m_cc, wm_cc = self.distros.Random_and_Weight('sgpi0_cc_4')
				pi0s_cc, wd_cc = self.distros.Random_and_Weight('pi0_signal_3_cc')
				pi0d_cc, ws_cc = self.distros.Random_and_Weight('pi0_delta_3_cc')
				pi0m_nc, wm_nc = self.distros.Random_and_Weight('sgpi0_nc_4')
				pi0s_nc, ws_nc = self.distros.Random_and_Weight('pi0_signal_3_nc')
				pi0d_nc, wd_nc = self.distros.Random_and_Weight('pi0_delta_3_nc')

			if cc: 
				if self.Evis<0.25:
					pi0m_r = pi0m_cc
				else:
					pi0m_r = pi0m_cc
					pi0d_r = pi0d_cc
					pi0s_r = pi0s_cc
			else: 
				if self.Evis<0.25:
					pi0m_r = pi0m_nc
				else:
					pi0m_r = pi0m_nc
					pi0d_r = pi0d_nc
					pi0s_r = pi0s_nc

			if self.Evis<0.25:
				if pi0m_r < 100:
					pi0 = 0
				else:
					pi0 = 1
			elif self.Evis>=0.25:
				if pi0m_r < 100.:
					pi0 = 0
				else:
					L_pi0m = math.log(abs(wm_cc)) - math.log(abs(wm_nc))
					L_pi0d = math.log(abs(wd_cc)) - math.log(abs(wd_nc))
					L_pi0s = math.log(abs(ws_cc)) - math.log(abs(ws_nc))
					ll = L_pi0m + L_pi0d + L_pi0s
					if ll>0.1:
						pi0=1
					else:
						pi0=0

		itype=-1
		if self.NRing==1:
			if self.Evis>0.1 and self.Evis<1.330 and self.MERIP==2:
				if self.MuEdk==0 and pi0==0:
					itype=0
				elif self.MuEdk>0 and pi0==0:
					itype=1
				elif self.MuEdk==0 and pi0==1:
					itype=2
			elif self.Evis>0.2 and self.Evis<1.330 and self.MERIP==3:
				if self.MuEdk==0:
					itype=3
				elif self.MuEdk==1:
					itype=4
				elif self.MuEdk>1:
					itype=5
			elif self.Evis>=1.330 and self.MERIP==2:
				if self.MuEdk>0:
					itype=7
				elif self.MuEdk==0:
					itype=8
			elif self.Evis>=1.330 and self.MERIP==3:
				itype=9
		elif self.NRing>1:
			if self.MERIP==2 and self.Evis>1.330:
				if ip_r>-0.25:
					if ip_n>0:
						itype=10
					elif ip_n<0:
						itype=11
				else:
					itype=13
			elif self.MERIP==3 and self.MERMomentum>0.6:
				itype=12
		if self.NRing==2:
			self.Imass = pp.InvariantMass(self.MERMomentum, self.sMERMomentum, self.MERDirection, self.sMERDirection)
			if self.MERIP==2 and self.sMERIP==2 and self.Imass>=0.085 and self.Imass<=0.215:
				itype=6

		self.originalType = itype


		itype=-1
		if ntag:
			if self.NRing==1:
				if self.Evis>0.1 and self.Evis<1.330 and self.MERIP==2:
					if self.MuEdk==0 and pi0==0:
						if self.Neutrons==0:
							itype=1
						else:
							itype=2
					elif self.MuEdk>0 and pi0==0:
						itype=0
					elif self.MuEdk==0 and pi0==1:
						itype=3
				elif self.Evis>0.2 and self.Evis<1.330 and self.MERIP==3:
					if self.MuEdk==1 and self.Neutrons>0:
						itype=5
					else:
						itype=4
				elif self.Evis>=1.330 and self.MERIP==2:
					if self.MuEdk>0:
						itype=7
					elif self.MuEdk==0:
						if self.Neutrons==0:
							itype=8
						else:
							itype=9
				elif self.Evis>=1.330 and self.MERIP==3:
					if self.MuEdk==1 and self.Neutrons>0:
						itype=10
					else:
						itype=11

			elif self.NRing>1:
				if self.MERIP==2 and self.Evis>1.330:
					if ip_r>-0.25:
						if ip_n>0:
							itype=12
						elif ip_n<0:
							itype=13
					else:
						itype=14
				elif self.MERIP==3 and self.MERMomentum>0.6:
					itype=15

			if self.NRing==2:
				self.Imass = pp.InvariantMass(self.MERMomentum, self.sMERMomentum, self.MERDirection, self.sMERDirection)
				if self.MERIP==2 and self.sMERIP==2 and self.Imass>=0.085 and self.Imass<=0.215:
					itype=6

			self.Type = itype
		else:
			self.Type = self.originalType


	def DecayE(self, ipnu, cc, mode):
		recmuedk = 0
		meson_flag = pp.MesonProduction(self.MERIP,abs(mode))

		if self.MERIP==2 and self.Evis<1.33:
			if cc and ipnu>0 and meson_flag:
				recmuedk = self.distros.Random('sge_muedk_ccnu_meson')
			elif cc and ipnu<0 and meson_flag:
				recmuedk = self.distros.Random('sge_muedk_ccnub_meson')
			elif (not cc) and meson_flag:
				recmuedk = self.distros.Random('sge_muedk_nc_meson')
			else:
				recmuedk = 0
				dummy = np.random.rand()
				if ipnu==14 and cc and dummy<0.83:
					recmuedk = 1
				elif ipnu==-14 and cc and dummy<0.96:
					recmuedk = 1
		elif self.MERIP==2 and self.Evis>=1.33:
			if cc and ipnu==12:
				recmuedk = self.distros.Random('mge_muedk_ccnue')
			elif cc and ipnu==-12:
				recmuedk = self.distros.Random('mge_muedk_ccnuebar')
			elif cc and abs(ipnu)>=14:
				recmuedk = self.distros.Random('mge_muedk_ccnumus')
			else:
				recmuedk = self.distros.Random('mge_muedk_nc')
		elif self.MERIP==3 and self.Evis<1.33:
			if cc and ipnu>0 and meson_flag:
				recmuedk = self.distros.Random('sgm_muedk_ccnu_meson')
			elif cc and ipnu<0 and meson_flag:
				recmuedk = self.distros.Random('sgm_muedk_ccnub_meson')
			elif (not cc) and meson_flag:
				recmuedk = self.distros.Random('sgm_muedk_nc_meson')
			else:
				recmuedk = 0
				dummy = np.random.rand()
				if ipnu==14 and cc and dummy<0.83:
					recmuedk = 1
				elif ipnu==-14 and cc and dummy<0.96:
					recmuedk = 1
		elif self.MERIP==3 and self.Evis>=1.33:
			if cc and ipnu>0:
				recmuedk = self.distros.Random('mgm_muedk_ccnu')
			elif cc and ipnu<0:
				recmuedk = self.distros.Random('mgm_muedk_ccnub')
			else:
				recmuedk = self.distros.Random('mgm_muedk_nc')		
		self.MuEdk = recmuedk


	def TotalVariables(self):
		if self.NRing == 1:
			self.MER = 0
			self.sMER = 0
			self.Evis = self.Momentum
			self.MERMomentum = self.Momentum
			self.sMERMomentum = 0
			self.TotDir = self.Direction
			self.MERDirection = self.Direction
			self.sMERDirection = self.Direction
			self.MERIP = self.IP
			self.sMERIP = 0
		elif self.NRing > 1:
			nring = 0
			maxp = 0
			smaxp = 0
			mer = 0
			smer = 0
			direc = np.zeros(3)
			evis = 0
			for i,label in enumerate(self.RecoLabels):
				if label == 1:
					nring += 1
					evis += self.Momentum[i]
					direc += self.Direction[i,:] * self.Momentum[i]
					if self.Momentum[i] > maxp:
						smaxp = maxp
						smer = mer
						maxp = self.Momentum[i]
						mer = i
			self.MER = mer
			self.sMER = smer
			self.MERMomentum = self.Momentum[mer]
			self.sMERMomentum = self.Momentum[smer]
			self.MERDirection = self.Direction[mer,:]
			self.sMERDirection = self.Direction[smer,:]
			self.Evis = evis
			if np.linalg.norm(direc)==0:
				self.NRing = 0
				self.TotDir = np.zeros(3)
			else:
				self.NRing = nring
				self.TotDir = direc / np.linalg.norm(direc)
			self.MERIP = self.IP[mer]
			self.sMERIP = self.IP[smer]


	def TooClose(self, pdg, P, dirv):
		mergeP = P
		mergeDir = dirv
		mergepdg = pdg
		for i in range(self.NRing-1):
			for j in range(i+1,self.NRing):
				cos = np.dot(self.Direction[i],self.Direction[j])
				if cos > self.AngResoThr:
					if self.Momentum[j]<self.Momentum[i]: 
						self.IP[j] = self.IP[i]
						pdg[j] = pdg[i]
					else: self.IP[j] = 2
					
					mergeP[j] = P[i] + P[j]
					mergeDir[j,:] =  dirv[i,:] + dirv[i,:]
					newNorm = np.linalg.norm(mergeDir[j,:])
					mergeDir[j,:] = mergeDir[j,:] / newNorm
					self.RecoLabels[i] = 0

		return mergepdg, mergeP, mergeDir
				

	def ReconDirection(self, dirv):
		ang=2
		self.Direction = dirv
		for i,ip in enumerate(self.IP):
			if  self.RecoLabels[i]==1:
				if ip==2 and self.Momentum[i]<1.3:
					ang = self.distros.Random('ang_sge')
				elif ip==3 and self.Momentum[i]<1.3:
					ang = self.distros.Random('ang_sgm')
				elif ip==2 and self.Momentum[i]>=1.3:
					ang = self.distros.Random('ang_mge')
				elif ip==3 and self.Momentum[i]>=1.3:
					ang = self.distros.Random('ang_mgm')
				ang=ang*pi/180.
			# Rodrigues way
				u = ap.RndVector()
				if self.NRing == 1:
					recodir = ap.RodRot(dirv,u,ang)
					self.Direction = recodir
				elif self.NRing > 1:
					recodir = ap.RodRot(dirv[i,:],u,ang)
					self.Direction[i,:] = recodir
				else:
					self.Direction = 9999*np.ones(3)


	def ReconMomentum(self, pdg, P):
		self.Momentum = P
		for i, p in enumerate(P):
			if  self.RecoLabels[i]==1 and p>0:
				if self.IP[i] == 2:
					reso = 0.01*(0.6+2.6/math.sqrt(p))
					bias = 0.0048 * p - 0.00072
					mHypothesis = Particle.from_pdgid(11).mass * self.MeVtoGeV
				else:
					reso = 0.01*(1.7+0.7/math.sqrt(p))
					bias = 0.0025 * p + 0.0015
					mHypothesis = Particle.from_pdgid(13).mass * self.MeVtoGeV
				
				p = (1-random.gauss(bias,reso)) * p
				pCorrSq = p**2+(Particle.from_pdgid(pdg[i]).mass * self.MeVtoGeV)**2-mHypothesis**2
				if pCorrSq < 0:
					self.RecoLabels[i] = 0
					self.Momentum[i] = 0
				else:
					self.Momentum[i] = math.sqrt(pCorrSq)


	def ReconMRIP(self, mer, ip, pdg):
		if self.CC and abs(pdg)==13:
			pid = self.distros.Random('mr_mer_pid_ccmu')
		else:
			pid = self.distros.Random('mr_mer_pid_rest')
		if pid>0:
			self.IP[mer]=3
		else:
			self.IP[mer]=2


	def ReconSRIP(self, ip, P):
		mis = np.random.rand()
		self.IP = np.array(ip)

		for i,p in enumerate(P):
			if p<0.2:
				if ip[i]==2:
					skmis = 0.0045
				else:#ip[i]==3
					skmis = 0.019
			elif p<0.3:
				if ip[i]==2:
					skmis = 0.015
				else:#ip[i]==3
					skmis = 0.011
			elif p<0.4:
				if ip[i]==2:
					skmis = 0.018
				else:#ip[i]==3
					skmis = 0.007
			elif p<0.5:
				if ip[i]==2:
					skmis = 0.015
				else:#ip[i]==3
					skmis = 0.006
			elif p<0.6:
				if ip[i]==2:
					skmis = 0.014
				else:#ip[i]==3
					skmis = 0.0045
			elif p<0.7:
				if ip[i]==2:
					skmis = 0.011
				else:#ip[i]==3
					skmis = 0.006
			elif p<0.8:
				if ip[i]==2:
					skmis = 0.0065
				else:#ip[i]==3
					skmis = 0.0075
			elif p<0.9:
				if ip[i]==2:
					skmis = 0.0065
				else:#ip[i]==3
					skmis = 0.0075
			elif p<1.0:
				if ip[i]==2:
					skmis = 0.005
				else:#ip[i]==3
					skmis = 0.009
			elif p<1.1:
				if ip[i]==2:
					skmis = 0.005
				else:#ip[i]==3
					skmis = 0.009
			elif p<1.2:
				if ip[i]==2:
					skmis = 0.006
				else:#ip[i]==3
					skmis = 0.0085
			elif p<1.3:
				if ip[i]==2:
					skmis = 0.0045
				else:#ip[i]==3
					skmis = 0.0075
			elif p<1.4:
				if ip[i]==2:
					skmis = 0.003
				else:#ip[i]==3
					skmis = 0.006
			if p>=1.4:
				if ip[i]==2:
					skmis = 0.004
				else:#ip[i]==3
					skmis = 0.008

			if mis<skmis:
				if ip[i]==2:
					self.IP[i]=3
				else:
					self.IP[i]=2

	def HNeutrons(self, ipnu, cc):
		nn = 0
		if self.MERIP==2 and self.Evis<1.33:
			if cc and ipnu>0:
				nn = self.distros.Random('sge_nh_nucc')
			elif cc and ipnu<0:
				nn = self.distros.Random('sge_nh_nubarcc')
			elif not cc:
				nn = self.distros.Random('sge_nh_nc')
		elif self.MERIP==2 and self.Evis>=1.33:
			if cc and ipnu>0:
				nn = self.distros.Random('mge_nh_nucc')
			elif cc and ipnu<0:
				nn = self.distros.Random('mge_nh_nubarcc')
			elif not cc:
				nn = self.distros.Random('mge_nh_nc')
		elif self.MERIP==3 and self.Evis<1.33:
			if cc and ipnu>0:
				nn = self.distros.Random('sgm_nh_nucc')
			elif cc and ipnu<0:
				nn = self.distros.Random('sgm_nh_nubarcc')
			elif not cc:
				nn = self.distros.Random('sgm_nh_nc')			
		elif self.MERIP==3 and self.Evis>1.33:
			if cc and ipnu>0:
				nn = self.distros.Random('mgm_nh_nucc')
			elif cc and ipnu<0:
				nn = self.distros.Random('mgm_nh_nubarcc')
			elif not cc:
				nn = self.distros.Random('mgm_nh_nc')
		self.Neutrons = nn

	def GdNeutrons(self, ipnu, cc):
		nn = 0
		if self.MERIP==2 and self.Evis<1.33:
			if cc and ipnu>0:
				nn = self.distros.Random('sge_ngd_nucc')
			elif cc and ipnu<0:
				nn = self.distros.Random('sge_ngd_nubarcc')
			elif not cc:
				nn = self.distros.Random('sge_ngd_nc')
		elif self.MERIP==2 and self.Evis>=1.33:
			if cc and ipnu>0:
				nn = self.distros.Random('mge_ngd_nucc')
			elif cc and ipnu<0:
				nn = self.distros.Random('mge_ngd_nubarcc')
			elif not cc:
				nn = self.distros.Random('mge_ngd_nc')
		elif self.MERIP==3 and self.Evis<1.33:
			if cc and ipnu>0:
				nn = self.distros.Random('sgm_ngd_nucc')
			elif cc and ipnu<0:
				nn = self.distros.Random('sgm_ngd_nubarcc')
			elif not cc:
				nn = self.distros.Random('sgm_ngd_nc')			
		elif self.MERIP==3 and self.Evis>1.33:
			if cc and ipnu>0:
				nn = self.distros.Random('mgm_ngd_nucc')
			elif cc and ipnu<0:
				nn = self.distros.Random('mgm_ngd_nubarcc')
			elif not cc:
				nn = self.distros.Random('mgm_ngd_nc')
		self.Neutrons = nn
