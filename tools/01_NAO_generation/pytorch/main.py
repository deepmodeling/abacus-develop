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

print("seed:",torch.initial_seed())
time_start = time.time()

file_list, info, C_init_info, V_info = IO.read_json.read_json("input.json")

QI,SI,VI = IO.read_QSV.read_file(info,file_list,V_info)
print(info)

if C_init_info["init_from_file"]:
	C = IO.func_C.read_C_init( C_init_info["C_init_file"], info )
else:
	C = IO.func_C.random_C_init(info)
E = orbital.set_E(info,info.Rcut)
orbital.normalize( orbital.generate_orbital(info,C,E,info.Rcut,info.dr), info.dr,C,flag_norm_C=True)

opt_orb = opt_orbital.Opt_Orbital()

#opt = torch.optim.Adam(sum( ([c.real,c.imag] for c in sum(C,[])), []), lr=info.lr, eps=1e-8)
opt = torch.optim.Adam( sum(C.values(),[]), lr=info.lr, eps=1e-20)

with open("Spillage.dat","w") as S_file:

	loss_old = np.inf

	for i in range(20000):
		Q = opt_orb.cal_Q(QI,C,info)
		S = opt_orb.cal_S(SI,C,info)
		V = opt_orb.cal_V(Q,S,info,V_info)				

		V_delta = ( torch.abs(VIi-Vi) for VIi,Vi in zip(VI,V) )		# abs or **2?
		if V_info["same_band"]:
			Spillage = sum( Vi[:info.Nb_true[ist]].sum() * info.weight[ist] for ist,Vi in enumerate(V_delta) ) / sum( Nb_true*weight for Nb_true,weight in zip(info.Nb_true,info.weight) )
		else:
			Spillage = sum( Vi[:info.Nb_true[ist],:info.Nb_true[ist]].sum() * info.weight[ist] for ist,Vi in enumerate(V_delta) ) / sum( Nb_true**2*weight for Nb_true,weight in zip(info.Nb_true,info.weight) )
		T = opt_orb.cal_T(C,E)
		if not "TSrate" in vars():	TSrate = torch.abs(0.002*Spillage/T).data[0]
		Loss = Spillage + TSrate*T

		print(Spillage.item(),T.item(),Loss.item(),file=S_file)
		if Loss.item() < loss_old:
			loss_old = Loss.item()
			C_old = IO.func_C.copy_C(C,info)
			flag_finish = 0
		else:
			flag_finish += 1
			if flag_finish > 50:	break

		opt.zero_grad()
		Loss.backward()
		opt.step()
#		orbital.normalize( orbital.generate_orbital(info,C,E,info.Rcut,info.dr), info.dr,C,flag_norm_C=True)

	orb = orbital.generate_orbital(info,C_old,E,info.Rcut,info.dr)
	orbital.smooth_orbital(orb,info.Rcut,info.dr,0.1)
	orbital.orth(orb,info.dr)
	IO.print_orbital.print_orbital(orb,info)
	IO.print_orbital.plot_orbital(orb,info.Rcut,info.dr)

	IO.func_C.write_C("C.dat",info,C_old)

	print("time:\t",time.time()-time_start)