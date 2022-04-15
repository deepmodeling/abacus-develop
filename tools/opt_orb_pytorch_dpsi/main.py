#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import IO.read_QSV
import IO.print_QSV
import IO.func_C
import IO.read_json
import IO.print_orbital
import opt_orbital
import orbital
import torch
import numpy as np
import time
import torch_optimizer 
import IO.cal_weight
import util
import IO.change_info
import pprint

def main():
	seed = int(1000*time.time())%(2**32)
	np.random.seed(seed)
	print("seed:",seed)
	time_start = time.time()

	file_list, info_true, weight_info, C_init_info, V_info = IO.read_json.read_json("INPUT")

	weight = IO.cal_weight.cal_weight(weight_info, V_info["same_band"], file_list["origin"])

	info_kst = IO.read_QSV.read_file_head(info_true,file_list["origin"])

	info_stru, info_element, info_opt = IO.change_info.change_info(info_kst,weight)
	info_max = IO.change_info.get_info_max(info_stru, info_element)

	print("info_kst:", info_kst, sep="\n", end="\n"*2, flush=True)
	print("info_stru:", pprint.pformat(info_stru), sep="\n", end="\n"*2, flush=True)
	print("info_element:", pprint.pformat(info_element,width=40), sep="\n", end="\n"*2, flush=True)
	print("info_opt:", pprint.pformat(info_opt,width=40), sep="\n", end="\n"*2, flush=True)
	print("info_max:", pprint.pformat(info_max), sep="\n", end="\n"*2, flush=True)

	QI,SI,VI_origin = IO.read_QSV.read_QSV(info_stru, info_element, file_list["origin"], V_info)
	if "linear" in file_list.keys():
		QI_linear, SI_linear, VI_linear = list(zip(*( IO.read_QSV.read_QSV(info_stru, info_element, file, V_info) for file in file_list["linear"] )))

	if C_init_info["init_from_file"]:
		C, C_read_index = IO.func_C.read_C_init( C_init_info["C_init_file"], info_element )
	else:
		C = IO.func_C.random_C_init(info_element)
	E = orbital.set_E(info_element)
	orbital.normalize(
		orbital.generate_orbital(info_element,C,E),
		{it:info_element[it].dr for it in info_element},
		C, flag_norm_C=True)

	opt_orb = opt_orbital.Opt_Orbital()

	#opt = torch.optim.Adam(sum( ([c.real,c.imag] for c in sum(C,[])), []), lr=info_opt.lr, eps=1e-8)
	#opt = torch.optim.Adam( sum(C.values(),[]), lr=info_opt.lr, eps=1e-20, weight_decay=info_opt.weight_decay)
	#opt = radam.RAdam( sum(C.values(),[]), lr=info_opt.lr, eps=1e-20 )
	opt = torch_optimizer.SWATS( sum(C.values(),[]), lr=info_opt.lr, eps=1e-20 )


	with open("Spillage.dat","w") as S_file:

		print( "\nSee \"Spillage.dat\" for detail status: " , flush=True )
		if info_opt.cal_T:
			print( '%5s'%"istep", "%20s"%"Spillage", "%20s"%"T.item()", "%20s"%"Loss", flush=True )
		else:
			print( '%5s'%"istep", "%20s"%"Spillage", flush=True )

		loss_old = np.inf
		for istep in range(10000):

			Spillage = 0
			for ist in range(len(info_stru)):

				Q = opt_orb.change_index_Q(opt_orb.cal_Q(QI[ist],C,info_stru[ist],info_element),info_stru[ist])
				S = opt_orb.change_index_S(opt_orb.cal_S(SI[ist],C,info_stru[ist],info_element),info_stru[ist],info_element)
				coef = opt_orb.cal_coef(Q,S)
				V = opt_orb.cal_V(coef,Q)
				V_origin = opt_orb.cal_V_origin(V,V_info)

				if "linear" in file_list.keys():
					V_linear = [None] * len(file_list["linear"])
					for i in range(len(file_list["linear"])):
						Q_linear = opt_orb.change_index_Q(opt_orb.cal_Q(QI_linear[i][ist],C,info_stru[ist],info_element),info_stru[ist])
						S_linear = opt_orb.change_index_S(opt_orb.cal_S(SI_linear[i][ist],C,info_stru[ist],info_element),info_stru[ist],info_element)
						V_linear[i] = opt_orb.cal_V_linear(coef,Q_linear,S_linear,V,V_info)

				def cal_Spillage(V_delta):
					Spillage = (V_delta * weight[ist][:info_stru[ist].Nb_true]).sum()
					return Spillage

				def cal_delta(VI, V):
					return ((VI[ist]-V)/util.update0(VI[ist])).abs()		# abs or **2?
				
				Spillage += 2*cal_Spillage(cal_delta(VI_origin,V_origin))
				if "linear" in file_list.keys():
					for i in range(len(file_list["linear"])):
						Spillage += cal_Spillage(cal_delta(VI_linear[i],V_linear[i]))

			if info_opt.cal_T:
				T = opt_orb.cal_T(C,E)
				if not "TSrate" in vars():	TSrate = torch.abs(0.002*Spillage/T).data[0]
				Loss = Spillage + TSrate*T
			else:
				Loss = Spillage

			if info_opt.cal_T:
				print_content = [istep, Spillage.item(), T.item(), Loss.item()]
			else:
				print_content = [istep, Spillage.item()]
			print(*print_content, sep="\t", file=S_file, flush=True)
			if not istep%100:
				print(*print_content, sep="\t", flush=True)

			if Loss.item() < loss_old:
				loss_old = Loss.item()
				C_old = IO.func_C.copy_C(C,info_element)
				flag_finish = 0
			else:
				flag_finish += 1
				if flag_finish > 50:
					break

			opt.zero_grad()
			Loss.backward()		
			if C_init_info["init_from_file"] and not C_init_info["opt_C_read"]:
				for it,il,iu in C_read_index:
					C[it][il].grad[:,iu] = 0
			opt.step()
			#orbital.normalize(
			#	orbital.generate_orbital(info_element,C,E),
			#	{it:info_element[it].dr for it in info_element},
			#	C, flag_norm_C=True)

	orb = orbital.generate_orbital(info_element,C_old,E)
	if info_opt.cal_smooth:
		orbital.smooth_orbital(
			orb,
			{it:info_element[it].Rcut for it in info_element}, {it:info_element[it].dr for it in info_element},
			0.1)
	orbital.orth(
		orb,
		{it:info_element[it].dr for it in info_element})
	IO.print_orbital.print_orbital(orb,info_element)
	IO.print_orbital.plot_orbital(
		orb,
		{it:info_element[it].Rcut for it in info_element},
		{it:info_element[it].dr for it in info_element})

	IO.func_C.write_C("ORBITAL_RESULTS.txt",C_old,Spillage)

	print("Time (PyTorch):     %s\n"%(time.time()-time_start), flush=True )


if __name__=="__main__":
	import sys
	np.set_printoptions(threshold=sys.maxsize, linewidth=10000)
	print( sys.version, flush=True ) 
	main()
