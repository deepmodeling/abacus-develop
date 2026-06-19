from util import ND_list
import util
import functools
import itertools
import torch

class Opt_Orbital:
		
	def cal_Q(self,QI,C,info_stru,info_element):
		"""
		  <\psi|\phi> = <\psi|jY> * <jY|\phi>
		  Q[it][il][ib,ia*im*iu]
		  	= sum_{q} QI[it][il][ib*ia*im,ie] * C[it][il][ie,iu]
		"""
		Q = dict()
		for it in info_stru.Na.keys():
			Q[it] = ND_list(info_element[it].Nl)

		for it in info_stru.Na.keys():
			for il in range(info_element[it].Nl):
				Q[it][il] = torch.mm( QI[it][il], C[it][il].to(torch.complex128) ).view(info_stru.Nb_true,-1)
		return Q

		
		
	def cal_S(self,SI,C,info_stru,info_element):
		"""
		  <\phi|\phi> = <\phi|jY> * <jY|jY> * <jY|\phi>
		  S[it1,it2][il1][il2][ia1*im1*iu1,ia2*im2*iu2]
		  	= sum_{ie1 ie2} C^*[it1][il1][ie1,iu1] * SI[it1,it2][il1][il2][ia1,im1,ie1,ia2,im2,ie2] * C[it2][[il2][ie2,iu2]
		"""
		S = dict()
		for it1,it2 in itertools.product( info_stru.Na.keys(), info_stru.Na.keys() ):
			S[it1,it2] = ND_list(info_element[it1].Nl, info_element[it2].Nl)

		for it1,it2 in itertools.product( info_stru.Na.keys(), info_stru.Na.keys() ):
			for il1,il2 in itertools.product( range(info_element[it1].Nl), range(info_element[it2].Nl) ):
				# SI_C[ia1*im1*ie1*ia2*im2,iu2]
				SI_C = torch.mm( 
					SI[it1,it2][il1][il2].view(-1,info_element[it2].Ne), 
					C[it2][il2].to(torch.complex128) )
				# SI_C[ia1*im1,ie1,ia2*im2*iu2]
				SI_C = SI_C.view( info_stru.Na[it1]*util.Nm(il1), info_element[it1].Ne, -1 )
				# Ct[iu1,ie1]
				Ct = C[it1][il1].t().to(torch.complex128)
				C_mm = functools.partial(torch.mm,Ct)
				# C_SI_C[ia1*im1][iu1,ia2*im2*iu2]
				C_SI_C = list(map( C_mm, SI_C ))
				# C_SI_C[ia1*im1*iu1,ia2*im2*iu2]
				C_SI_C = torch.cat( C_SI_C, dim=0 )
#???				C_SI_C = C_SI_C.view(info_stru.Na[it1]*util.Nm(il1)*info_element[it1].Nu[il1],-1)
				S[it1,it2][il1][il2] = C_SI_C
		return S
		
		
		
	def change_index_S(self,S,info_stru,info_element):							# S[it1,it2][il1][il2][ia1*im1*iu1,ia2*im2*iu2]
		"""
		  <\phi|\phi>
		  S_cat[it1*il1*iat*im1*iu1,iat2*il2*ia2*im2*iu2]
		"""
		# S_[it1][il1*ia1*im1*iu1,it2*il2*ia2*im2*iu2]
		S_ = dict()
		for it1 in info_stru.Na.keys():
			# S_t[it2][il1*ia1*im1*iu1,il2*ia2*im2*iu2]
			S_t = dict()
			for it2 in info_stru.Na.keys():
				# S_tt[il1][ia1*im1*iu1,il2*ia2*im2*iu2]
				S_tt = ND_list(info_element[it1].Nl)
				for il1 in range(info_element[it1].Nl):
					S_tt[il1] = torch.cat( S[it1,it2][il1], dim=1 )
				S_t[it2] = torch.cat( S_tt, dim=0 )
			S_[it1] = torch.cat( list(S_t.values()), dim=1 )
		# S_cat[it1*il1*ia1*im1*iu1,it2*il2*ia2*im2*iu2]
		S_cat = torch.cat( list(S_.values()), dim=0 )	
		return S_cat
		
		
		
	def change_index_Q(self,Q,info_stru):					# Q[it][il][ib,ia*im*iu]
		"""
		  <\psi|\phi>
		  Q_cat[ib,it*il*ia*im*iu]
		"""
		# Q_b[ib][0,it*il*ia*im*iu]
		Q_b = ND_list(info_stru.Nb_true)
		for ib in range(info_stru.Nb_true):
			# Q_[it][il*ia*im*iu]
			Q_ = dict()
			for it in info_stru.Na.keys():
				# Q_ts[il][ia*im*iu]
				Q_ts = [ Q_tl[ib] for Q_tl in Q[it] ]
				Q_[it] = torch.cat(Q_ts)
			Q_b[ib] = torch.cat(list(Q_.values())).view(1,-1)
		# Q_cat[ib,it*il*ia*im*iu]
		Q_cat = torch.cat( Q_b, dim=0 )
		return Q_cat
			
			
	
	def cal_coef(self,Q,S):
		# Q[ib,it*il*ia*im*iu]
		# S[it1*il1*ia1*im1*iu1,it2*il2*ia2*im2*iu2]
		"""
		  <\psi|\phi> * <\phi|\phi>^{-1}
		  coef[ib,it*il*ia*im*iu]
			= Q[ib,it1*il1*ia1*im1*iu1] * S{[it1*il1*ia1*im1*iu1,it2*il2*ia2*im2*iu2]}^{-1}
		"""
		S_I = torch.inverse(S)
		coef = torch.mm(Q, S_I)
		return coef
		
		
		
	def cal_V(self,coef,Q):
		# coef[ib,it*il*ia*im*iu]
		# Q[ib,it*il*ia*im*iu]
		"""
		  <\psi|\psi> = <\psi|\phi> * <\phi|\phi>^{-1} * <\phi|psi>
		  V[ib1,ib2]
		  	= sum_{it1,ia1,il1,im1,iu1} sum_{it2,ia2,il2,im2,iu2}
		  	Q[ib1,it1*il1*ia1*im1*iu1] * S{[it1*il1*ia1*im1*iu1,it2*il2*ia2*im2*iu2]}^{-1} * Q[ib2,it2*il2*ia2*im2*iu2]
		"""
		V = torch.mm( coef, Q.t().conj() ).real
		return V


	def cal_V_origin(self,V,V_info):
		# V[ib1,ib2]
		"""
		  <\psi|\psi> = <\psi|\phi> * <\phi|\phi>^{-1} * <\phi|psi>
		  V_origin[ib]	
		  V_origin[ib1,ib2]
		"""			
		if V_info["same_band"]:		V_origin = V.diag().sqrt()
		else:						V_origin = V.sqrt()
		return V_origin		
		
		
	def cal_V_linear(self,coef,Q_linear,S_linear,V,V_info):
		# coef[ib,it*il*ia*im*iu]
		# Q_linear[ib,it*il*ia*im*iu]
		# S_linear[it1*il1*ia1*im1*iu1,it2*il2*ia2*im2*iu2]
		# V[ib1,ib2]
		"""
		  V_linear[ib]
		  V_linear[ib1,ib2]
		"""		
		V_linear_1 = coef.mm(S_linear).mm(coef.t().conj()).real
		V_linear_2 = Q_linear.mm(coef.t().conj()).real
		V_linear_3 = coef.mm(Q_linear.t().conj()).real
		if V_info["same_band"]:
			V_linear_1 = V_linear_1.diag()
			V_linear_2 = V_linear_2.diag()
			V_linear_3 = V_linear_3.diag()
		if V_info["same_band"]:		Z = V.diag().sqrt()
		else:						Z = V.sqrt()
		Z = util.update0(Z)
		V_linear = (-V_linear_1/Z + V_linear_2 + V_linear_3) / Z
		return V_linear
			
	
	def cal_T(self,C,E):
		""" T = 0.5* sum_{it,il,iu} sum_{ie} ( E[it][il,ie] * C[it][il][ie,iu] )**2 """
		T = torch.zeros(1)
		num = 0
		for it,C_t in C.items():
			for il,C_tl in enumerate(C_t):
				for iu in range(C_tl.size()[1]):
					T_tlu = torch.zeros(1)
					Z_tlu = 0
					for ie in range(C_tl.size()[0]):
						T_tlu = T_tlu + ( E[it][il,ie] * C_tl[ie,iu] )**2
						Z_tlu = Z_tlu + E[it][il,ie].item()**2
					T = T + T_tlu/Z_tlu
				num += C_tl.size()[1]
		T = 0.5 * T / num
		return T