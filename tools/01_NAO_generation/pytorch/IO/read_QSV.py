from util import *
import torch
import torch_complex
import itertools
import numpy as np

import time

def read_file(info,file_list,V_info):
	""" QI[ist][it][il][ib*ia*im,ie] """
	""" SI[ist][it1][it2][il1][il2][ie1,ia1,im1,ia2,im2,ie2] """
	""" VI[ist][ib] """
	import copy
	info_true = copy.deepcopy(info)
	info_true.Nst = len(file_list)
	info_true.Nt = ND_list(info_true.Nst,element="list()")
	info_true.Na = ND_list(info_true.Nst,element="dict()")
	info_true.Nb = ND_list(info_true.Nst)
	info_true.Nk = ND_list(info_true.Nst)
	info_true.Ne = dict()
	info_true.weight = ND_list(info_true.Nst,element="list()")
	QI=[];	SI=[];	VI=[]

	for ist_true,file_name in enumerate(file_list):
		print(file_name)
		with open(file_name,"r") as file:

			ignore_line(file,4)
			Nt_tmp = int(file.readline().split()[0])
			for it in range(Nt_tmp):
				t_tmp = file.readline().split()[0]
				assert t_tmp in info.Nt_all
				info_true.Nt[ist_true].append( t_tmp )
				info_true.Na[ist_true][t_tmp] = int(file.readline().split()[0])
				ignore_line( file, info_true.Na[ist_true][t_tmp] )
			ignore_line(file,6)
			Nl_ist = int(file.readline().split()[0])+1
			for it,Nl_C in info.Nl.items():
				print(it,Nl_ist,Nl_C)
				assert Nl_ist>=Nl_C
				info_true.Nl[it] = Nl_ist
			info_true.Nk[ist_true] = int(file.readline().split()[0])
			info_true.Nb[ist_true] = int(file.readline().split()[0])
			ignore_line(file,1)
			Ne_tmp = list(map(int,file.readline().split()[:Nt_tmp]))
			for it,Ne in zip(info_true.Nt[ist_true],Ne_tmp):
				assert info_true.Ne.setdefault(it,Ne)==Ne
			ignore_line(file,1)
			for ik in range(info_true.Nk[ist_true]):
				info_true.weight[ist_true].append(float(file.readline().split()[-1]))
			ignore_line(file,3)

			line = None
			for ik in range(info_true.Nk[ist_true]):
				print("read QI:",ist_true,ik)
				qi,line = read_QI(info_true,ist_true,file,line)
				QI.append( qi )
			ignore_line(file,3)
			for ik in range(info_true.Nk[ist_true]):
				print("read SI:",ist_true,ik)
				si,line = read_SI(info_true,ist_true,file,line)
				SI.append( si )
			if V_info["init_from_file"]:
				ignore_line(file,3)
			for ik in range(info_true.Nk[ist_true]):
				print("read VI:",ist_true,ik)
				vi,line = read_VI(info_true,V_info,ist_true,file,line)
				VI.append( vi )

	info.Nst = sum(info_true.Nk,0)
	import itertools
	repeat_Nk = lambda x: list( itertools.chain.from_iterable( map( lambda x:itertools.repeat(*x), zip(x,info_true.Nk) ) ) )
	info.Nt = repeat_Nk(info_true.Nt)
	info.Na = repeat_Nk(info_true.Na)
	info.Nb = repeat_Nk(info_true.Nb)
	info.Nb_true = repeat_Nk(info.Nb_true)
	info.Ne = info_true.Ne
	info.weight = list(itertools.chain.from_iterable([np.array(wk)*ws for wk,ws in zip(info_true.weight,info.weight)]))	# info_true.weight[ist][ik] * info.weight[ist]

	return QI,SI,VI




def read_QI(info,ist,file,line=None):
	""" QI[it][il][ib*ia*im,ie] """
	QI = dict()
	for it in info.Nt[ist]:
		QI[it] = ND_list(info.Nl[it])
		for il in range(info.Nl[it]):
			QI[it][il] = torch_complex.ComplexTensor(
				np.empty((info.Nb[ist],info.Na[ist][it],info.Nm(il),info.Ne[it]),dtype=np.float32),
				np.empty((info.Nb[ist],info.Na[ist][it],info.Nm(il),info.Ne[it]),dtype=np.float32) )
	for ib in range(info.Nb[ist]):
		for it in info.Nt[ist]:
			for ia in range(info.Na[ist][it]):
				for il in range(info.Nl[it]):
					for im in range(info.Nm(il)):
						for ie in range(info.Ne[it]):
							if not line:	line = file.readline().split()
							QI[it][il].real[ib,ia,im,ie] = float(line.pop(0))
							if not line:	line = file.readline().split()
							QI[it][il].imag[ib,ia,im,ie] = float(line.pop(0))
	for it in info.Nt[ist]:
		for il in range(info.Nl[it]):
			QI[it][il] = torch_complex.ComplexTensor(
				torch.from_numpy(QI[it][il].real).view(-1,info.Ne[it]),
				torch.from_numpy(QI[it][il].imag).view(-1,info.Ne[it]))				
	return QI,line



def read_SI(info,ist,file,line=None):
	""" SI[it1,it2][il1][il2][ie1,ia1,im1,ia2,im2,ie2] """
	SI = dict()
	for it1,it2 in itertools.product( info.Nt[ist], info.Nt[ist] ):
		SI[it1,it2] = ND_list(info.Nl[it1],info.Nl[it2])
		for il1,il2 in itertools.product( range(info.Nl[it1]), range(info.Nl[it2]) ):
			SI[it1,it2][il1][il2] = torch_complex.ComplexTensor(
				np.empty((info.Na[ist][it1],info.Nm(il1),info.Ne[it1],info.Na[ist][it2],info.Nm(il2),info.Ne[it2]),dtype=np.float32),
				np.empty((info.Na[ist][it1],info.Nm(il1),info.Ne[it1],info.Na[ist][it2],info.Nm(il2),info.Ne[it2]),dtype=np.float32) )
	for it1 in info.Nt[ist]:
		for ia1 in range(info.Na[ist][it1]):
			for il1 in range(info.Nl[it1]):
				for im1 in range(info.Nm(il1)):
					for it2 in info.Nt[ist]:
						for ia2 in range(info.Na[ist][it2]):
							for il2 in range(info.Nl[it2]):
								for im2 in range(info.Nm(il2)):
									for ie1 in range(info.Ne[it1]):
										for ie2 in range(info.Ne[it2]):
											if not line:	line = file.readline().split()
											SI[it1,it2][il1][il2].real[ia1,im1,ie1,ia2,im2,ie2] = float(line.pop(0))
											if not line:	line = file.readline().split()
											SI[it1,it2][il1][il2].imag[ia1,im1,ie1,ia2,im2,ie2] = float(line.pop(0))
	for it1,it2 in itertools.product( info.Nt[ist], info.Nt[ist] ):
		for il1,il2 in itertools.product( range(info.Nl[it1]), range(info.Nl[it2]) ):	
			SI[it1,it2][il1][il2] = torch_complex.ComplexTensor(
				torch.from_numpy(SI[it1,it2][il1][il2].real),
				torch.from_numpy(SI[it1,it2][il1][il2].imag))
	return SI,line



def read_VI(info,V_info,ist,file,line=None):
	if V_info["same_band"]:
		""" VI[ib] """
		if V_info["init_from_file"]:
			VI = np.empty(info.Nb[ist],dtype=np.float32)
			for ib in range(info.Nb[ist]):
				if not line:	line = file.readline().split()
				VI.data[ib] = float(line.pop(0))
		else:
			VI = np.ones(info.Nb[ist],dtype=np.float32)
	else:
		""" VI[ib1,ib2] """
		if V_info["init_from_file"]:
			VI = np.empty((info.Nb[ist],info.Nb[ist]),dtype=np.float32)
			for ib1,ib2 in itertools.product( range(info.Nb[ist]), range(info.Nb[ist]) ):
				if not line:	line = file.readline().split()
				VI[ib1,ib2] = float(line.pop(0))
		else:
			VI = np.eye(info.Nb[ist],info.Nb[ist],dtype=np.float32)
	return torch.from_numpy(VI),line
