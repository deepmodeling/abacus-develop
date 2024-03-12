def print_V(V,file_name):
	""" V[ist][ib] """
	with open(file_name,"w") as file:
		for V_s in V:
			for V_sb in V_s:
				print(1-V_sb.item(),end="\t",file=file)
			print(file=file)
			
def print_S(S,file_name):
	""" S[ist][it1,it2][il1][il2][ia1*im1*in1,ia2*im2*in2] """
	with open(file_name,"w") as file:
		for ist,S_s in enumerate(S):
			for (it1,it2),S_tt in S_s.items():
				for il1,S_ttl in enumerate(S_tt):
					for il2,S_ttll in enumerate(S_ttl):
						print(ist,it1,it2,il1,il2,file=file)
						print(S_ttll.real.numpy(),file=file)
						print(S_ttll.imag.numpy(),"\n",file=file)
						
def print_Q(Q,file_name):
	""" Q[ist][it][il][ib,ia*im*iu] """
	with open(file_name,"w") as file:
		for ist,Q_s in enumerate(Q):
			for it,Q_st in Q_s.items():
				for il,Q_stl in enumerate(Q_st):
					print(ist,it,il,file=file)
					print(Q_stl.real.numpy(),file=file)
					print(Q_stl.imag.numpy(),"\n",file=file)
					