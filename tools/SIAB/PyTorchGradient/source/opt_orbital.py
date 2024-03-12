from util import ND_list
import inverse
import torch_complex
import functools
import itertools
import torch

class Opt_Orbital:
		
	def cal_Q(self,QI,C,info):
		"""
		Q[ist][it][il][ib,ia*im*iu]
			= sum_{q} QI[ist][it][il][ib*ia*im,ie] * C[it][il][ie,iu]
		"""
		Q = ND_list(info.Nst,element="dict()")
		for ist in range(info.Nst):
			for it in info.Nt[ist]:
				Q[ist][it] = ND_list(info.Nl[it])

		for ist in range(info.Nst):
			for it in info.Nt[ist]:
				for il in range(info.Nl[it]):
					Q[ist][it][il] = torch_complex.mm( QI[ist][it][il], C[it][il] ).view(info.Nb[ist],-1)
		return Q

		
		
	def cal_S(self,SI,C,info):
		"""
		S[ist][it1,it2][il1][il2][ia1*im1*iu1,ia2*im2*iu2]
			= sum_{ie1 ie2} C^*[it1][il1][ie1,iu1] * SI[ist][it1,it2][il1][il2][ia1,im1,ie1,ia2,im2,ie2] * C[it2][[il2][ie2,iu2]
		"""
		S = ND_list(info.Nst,element="dict()")
		for ist in range(info.Nst):
			for it1,it2 in itertools.product( info.Nt[ist], info.Nt[ist] ):
				S[ist][it1,it2] = ND_list(info.Nl[it1],info.Nl[it2])

		for ist in range(info.Nst):
			for it1,it2 in itertools.product( info.Nt[ist], info.Nt[ist] ):
				for il1,il2 in itertools.product( range(info.Nl[it1]), range(info.Nl[it2]) ):
					# SI_C[ia1*im1*ie1*ia2*im2,iu2]
					SI_C = torch_complex.mm( 
						SI[ist][it1,it2][il1][il2].view(-1,info.Ne[it2]), 
						C[it2][il2] )
					# SI_C[ia1*im1,ie1,ia2*im2*iu2]
					SI_C = SI_C.view( info.Na[ist][it1]*info.Nm(il1), info.Ne[it1], -1 )
					# Ct[iu1,ie1]
					Ct = C[it1][il1].t()
					C_mm = functools.partial(torch_complex.mm,Ct)
					# C_SI_C[ia1*im1][iu1,ia2*im2*iu2]
					C_SI_C = list(map( C_mm, SI_C ))
					# C_SI_C[ia1*im1*iu1,ia2*im2*iu2]
					C_SI_C = torch_complex.cat( C_SI_C, dim=0 )
#???				C_SI_C = C_SI_C.view(info.Na[ist][it1]*info.Nm(il1)*info.Nu[it1][il1],-1)
					S[ist][it1,it2][il1][il2] = C_SI_C
		return S

		
		
	def cal_V(self,Q,S,info,V_info):
		"""
		V[ist][ib]
			= sum_{it1,ia1,il1,im1,iu1} sum_{it2,ia2,il2,im2,iu2}
			Q[ist][it1][il1][ib,ia1*im1*iu1] * S[ist]{[it1][it2][il1][il2][ia1*im1*iu1,ia2*im2*iu2]}^{-1} * Q[ist][it2][il2][ib,ia2*im2*iu2]		
		V[ist][ib1,ib2]
			= sum_{it1,ia1,il1,im1,iu1} sum_{it2,ia2,il2,im2,iu2}
			Q[ist][it1][il1][ib1,ia1*im1*iu1] * S[ist]{[it1][it2][il1][il2][ia1*im1*iu1,ia2*im2*iu2]}^{-1} * Q[ist][it2][il2][ib2,ia2*im2*iu2]
		"""
		V = ND_list(info.Nst)
		for ist in range(info.Nst):
		
			# S_s[it1][il1*ia1*im1*iu1,it2*il2*ia2*im2*iu2]
			S_s = dict()
			for it1 in info.Nt[ist]:
				# S_st[it2][il1*ia1*im1*iu1,il2*ia2*im2*iu2]
				S_st = dict()
				for it2 in info.Nt[ist]:
					# S_stt[il1][ia1*im1*iu1,il2*ia2*im2*iu2]
					S_stt = ND_list(info.Nl[it1])
					for il1 in range(info.Nl[it1]):
						S_stt[il1] = torch_complex.cat( S[ist][it1,it2][il1], dim=1 )
					S_st[it2] = torch_complex.cat( S_stt, dim=0 )
				S_s[it1] = torch_complex.cat( list(S_st.values()), dim=1 )
			# S_cat[it1*il1*iat*im1*iu1,iat2*il2*ia2*im2*iu2]
			S_cat = torch_complex.cat( list(S_s.values()), dim=0 )
			S_I = torch_complex.inverse(S_cat)
			
			if V_info["same_band"]:
				# V_s[ib]
				V_s = ND_list(info.Nb[ist])
				for ib in range(info.Nb[ist]):
					# Q_s[it][il*ia*im*iu]
					Q_s = dict()
					for it in info.Nt[ist]:
						# Q_ts[il][ia*im*iu]
						Q_ts = [ Q_stl[ib] for Q_stl in Q[ist][it] ]
						Q_s[it] = torch_complex.cat(Q_ts)
					# Q_cat[it*il*ia*im*iu]
					Q_cat = torch_complex.cat(list(Q_s.values()))
					V_s[ib] = torch_complex.dot( Q_cat.conj(), torch_complex.mv( S_I, Q_cat ) ).real.view(-1)
				# V[ist][ib]
				V[ist] = torch.cat(V_s)
			else:
				# Q_b[ib][0,it*il*ia*im*iu]
				Q_b = ND_list(info.Nb[ist])
				for ib in range(info.Nb[ist]):
					# Q_s[it][il*ia*im*iu]
					Q_s = dict()
					for it in info.Nt[ist]:
						# Q_ts[il][ia*im*iu]
						Q_ts = [ Q_stl[ib] for Q_stl in Q[ist][it] ]
						Q_s[it] = torch_complex.cat(Q_ts)
					Q_b[ib] = torch_complex.cat(list(Q_s.values())).view(1,-1)
				# Q_cat[ib,it*il*ia*im*iu]
				Q_cat = torch_complex.cat( Q_b, dim=0 )
				# V[ist][ib1,ib2]
				V[ist] = torch_complex.mm( Q_cat.conj(), torch_complex.mm( S_I, Q_cat.t() ) ).real

		return V
		
	
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