from util import *
import torch
import numpy as np

# C means Coefficient

def random_C_init(info):
	""" C[it][il][ie,iu] """
	C = dict()
	for it in info.Nt_all:
		C[it] = ND_list(info.Nl[it])
		for il in range(info.Nl[it]):
			C[it][il] = torch.tensor(np.random.uniform(-1,1, (info.Ne[it], info.Nu[it][il])), dtype=torch.float32, requires_grad=True)
	return C
	
	
	
def read_C_init(file_name,info):
	""" C[it][il][ie,iu] """
	C = random_C_init(info)

	with open(file_name,"r") as file:
	
		for line in file:
			if line.strip() == "<Coefficient>":	
				line=None
				break
		ignore_line(file,1)
	
		while True:
			line = file.readline().strip()
			if line.startswith("Type"):
				it,il,iu = list(map(int,file.readline().split()));	
				it=info.Nt_all[it-1];	iu-=1
				line = file.readline().split()
				for ie in range(info.Ne[it]):
					if not line:	line = file.readline().split()
					C[it][il].data[ie,iu] = float(line.pop(0))
			elif line.startswith("</Coefficient>"):
				break;
			else:
				raise IOError("unknown line in read_C_init "+file_name+"\n"+line)
	return C

	
	
def copy_C(C,info):
	C_copy = dict()
	for it in info.Nt_all:
		C_copy[it] = ND_list(info.Nl[it])
		for il in range(info.Nl[it]):
			C_copy[it][il] = C[it][il].clone()
	return C_copy
	
	
	
def write_C(file_name,info,C,Spillage):
	with open(file_name,"w") as file:
		print("<Coefficient>", file=file)
		nTotal = 0
		for it,C_t in C.items():
			for il,C_tl in enumerate(C_t):
				for iu in range(C_tl.size()[1]):
					nTotal += 1 
			#nTotal = sum(info["Nu"][it])
		#
		print("\t %s Total number of radial orbitals."%nTotal , file=file) 
		for it,C_t in C.items():
			for il,C_tl in enumerate(C_t):
				for iu in range(C_tl.size()[1]):
					print("\tType\tL\tZeta-Orbital", file=file)
					print(f"\t  {info.Nt_all.index(it)+1} \t{il}\t    {iu+1}", file=file)
					for ie in range(C_tl.size()[0]):
						print("\t", '%18.14f'%C_tl[ie,iu].item(), file=file)
		print("</Coefficient>", file=file)
		print("<Mkb>", file=file)
		print("Left spillage = %.10e"%Spillage.item(), file=file)
		print("</Mkb>", file=file)

	
	
#def init_C(info):
#	""" C[it][il][ie,iu] """
#	C = ND_list(max(info.Nt))
#	for it in range(len(C)):
#		C[it] = ND_list(info.Nl[it])
#		for il in range(info.Nl[it]):
#			C[it][il] = torch.autograd.Variable( torch.Tensor( info.Ne, info.Nu[it][il] ), requires_grad = True )
#			
#	with open("C_init.dat","r") as file:
#		line = []
#		for it in range(len(C)):
#			for il in range(info.Nl[it]):
#				for i_n in range(info.Nu[it][il]):
#					for ie in range(info.Ne[it]):
#						if not line:	line=file.readline().split()
#						C[it][il].data[ie,i_n] = float(line.pop(0))
#	return C
