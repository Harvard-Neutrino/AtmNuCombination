import numpy as np
from particle import Particle
import ParticleProperties as pp
import math
from math import log10, pi
# import ROOT as root
# import mulself.TrueIProcessing
import random
import applications as ap
# from Apps import TrueRingConstructor

class RecoRing:
	AngResoThr = 0.9 #math.cos(25.*pi/180.)
	MeVtoGeV = 0.001
	def __init__(self, TrueNRing, TrueRingPDG, TrueRingIP, TrueRingE, TrueRingP, TrueRingDir, distros, mode):
		
		self.distros = distros

		self.NRing = TrueNRing
		# self.IP = TrueRingIP
		CC = True if abs(mode)<30 else False

		self.RecoLabels = np.ones(TrueNRing)
		self.PiPlusWA(TrueRingPDG, TrueRingP)

		# Ring IP reconstruction
		self.ReconSRIP(TrueRingIP, TrueRingP)

		# Momentum reconstruction of each ring
		self.ReconMomentum(TrueRingPDG, TrueRingP)

		# Direction reconstruction of each ring
		self.ReconDirection(TrueRingDir)


		if self.NRing >1: # Oh, the multi-ringers! Special care needed
			self.ReconMRIP(TrueRingIP, CC, TrueRingPDG)
			mTrueRingPDG, mTrueRingP, mTrueRingDir = self.TooClose(TrueRingPDG, TrueRingP, TrueRingDir) # Put labels on rings which won't be reconstructed due to poor angular resolution
			# print(mTrueRingP)
			# print(mTrueRingPDG)
			# print(self.RecoLabels)
			if np.any(self.RecoLabels==0):
				self.ReconMomentum(mTrueRingPDG, mTrueRingP)
				self.ReconDirection(mTrueRingDir)
			# print(self.RecoDirection)
			# print(self.TrueDirection)

		self.TotalVariables()

		self.MuEdk = 0



	def PiPlusWA(self, pdg, p):
		for i,part in enumerate(pdg):
			if abs(part)>=211 and p[i]<0.3:
				self.RecoLabels[i] = 0

		# self.MER = np.argmax(self.RecoMomentum)
		# self.Evis = np.sum(self.RecoMomentum)

		# if self.RecoRing == 1: self.Direction = self.RecoDirection

	def SKType(self, ipnu, cc):
		
		if self.NRing>1 and self.MERIP==2:
			if cc:
				ip_r = self.distros.Random('other_cc')
				if ipnu>0:
					ip_n = self.distros.Random('mre_nunubar_ccnue')
				else:
					ip_n = self.distros.Random('mre_nunubar_ccnuebar')
				# ip_r = other.GetRandom()
				# ip_n = nunub.GetRandom()
			else:
				ip_r = self.distros.Random('other_nc')
				ip_n = self.distros.Random('mre_nunubar_nc')
				# ip_r = other.GetRandom()
				# ip_n = nunub.GetRandom()

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
					if ll>0:
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
			elif self.MERIP==3 and self.Evis>0.6 and self.MERMomentum>0.6:
				itype=12

		if self.NRing==2:
			invMass = pp.InvariantMass(self.MERMomentum, self.sMERMomentum, self.MERDirection, self.sMERDirection)
			if self.MERIP==2 and self.sMERIP==2 and invMass>=0.085 and invMass<=0.215:
				itype=6

		self.Type = itype




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
			self.TotDir = direc / np.linalg.norm(direc)
			self.Evis = evis
			self.NRing = nring
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
					# if abs(pdg[i]) < abs(pdg[j]): self.PDG[j] = self.PDG[i]


				

	def ReconDirection(self, dirv):
		ang=4
		self.Direction = dirv
		for i,ip in enumerate(self.IP):
			# if ip==2 and self.RecoMomentum[i]<1.3 and self.NRing==1: # SGE
			if  self.RecoLabels[i]==1:
				if ip==2 and self.Momentum[i]<1.3:
					ang = self.distros.Random('ang_sge')
				# elif ip==3 and self.RecoMomentum[i]<1.3 and self.NRing==1: # SGM
				elif ip==3 and self.Momentum[i]<1.3:
					ang = self.distros.Random('ang_sgm')
				# elif ip==2 and self.RecoMomentum[i]>=1.3 and self.NRing==1: # MGE
				elif ip==2 and self.Momentum[i]>=1.3:
					ang = self.distros.Random('ang_mge')
				# elif ip==3 and self.RecoMomentum[i]>=1.3 and self.NRing==1: # MGM
				elif ip==3 and self.Momentum[i]>=1.3:
					ang = self.distros.Random('ang_mgm')

				ang=ang*pi/180.
				# Rodrigues way
				u = ap.RndVector()
				if self.NRing == 1:
					self.Direction = ap.RodRot(dirv,u,ang)
				elif self.NRing > 1:
					# print(u, self.RecoDirection[i])
					self.Direction[i,:] = ap.RodRot(dirv[i,:],u,ang)
					# print(self.RecoDirection[i] )


	def ReconMomentum(self, pdg, P):
		self.Momentum = P
		for i, p in enumerate(P):
			if  self.RecoLabels[i]==1:
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
				# print(mHypothesis)

				if pCorrSq < 0:
					self.RecoLabels[i] = 0
					self.Momentum[i] = 0
				else:
					self.Momentum[i] = math.sqrt(pCorrSq)



	def ReconMRIP(self, ip, cc, pdg):
		self.IP = ip
		for i,part in enumerate(ip):
			# if ((part==2 and P[i]>1.33) or (part==3 and P[i]>0.6)):
			if cc and abs(pdg[i])==14:
				pid = self.distros.Random('mr_mer_pid_ccmu')
			else:
				pid = self.distros.Random('mr_mer_pid_rest')
			if pid>0:
				self.IP[i]=3
			else:
				self.IP[i]=2


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

