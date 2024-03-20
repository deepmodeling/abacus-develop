import os
import re
import numpy as np
import json
import pprint
import textwrap
import utils
import collections
import pathlib
from sklearn.cluster import KMeans
import read_stru

def cal_ABACUS(T1,T2,i_dis):
	folder = pathlib.Path(utils.folder_name(T1,T2,i_dis)).resolve()
	folder.mkdir(parents=True,exist_ok=False)
	
	
	with open(folder/"INPUT","w") as file:
		info = utils.read_info()
		input_dict = read_stru.get_input_dict()
		input_dict["ntype"] = 1 if T1==T2 else 2
		input_dict["exx_hybrid_type"] = 'opt_orb'
		input_dict["nbands"] = (read_stru.get_nw()[T1] if abs(i_dis)<1E-10 else read_stru.get_nw()[T1]+read_stru.get_nw()[T2])
		input_dict["nspin"] = 1
		input_dict["gamma_only"] = 1
		input_dict["pseudo_dir"] = os.path.abspath(input_dict.get("pseudo_dir",r"./"))
		input_dict["exx_opt_orb_lmax"] = len(info["Nu"])-1
		read_stru.print_input(file,input_dict)
		
	with open(folder/"STRU","w") as file:
		Ts = (T1,) if T1==T2 else (T1,T2)
		file.write("ATOMIC_SPECIES\n")
		pseudo_path = read_stru.get_pseudo_path()
		for T in Ts:
			file.write(f"{T}	1	{pseudo_path[T]}\n")
		file.write("\nNUMERICAL_ORBITAL\n")
		lcao_path = read_stru.get_lcao_path()
		for T in Ts:
			file.write(f"{lcao_path[T]}\n")
		file.write(textwrap.dedent(f"""
			LATTICE_CONSTANT
			1\n
			LATTICE_VECTORS
			30 0 0
			0 30 0
			0 0 30\n
			ATOMIC_POSITIONS
			Cartesian
			"""))
		if T1==T2:
			if abs(i_dis)<1E-10:
				file.write(textwrap.dedent(f"""
					{T1}
					0
					1
					0 0 0 0 0 0
					"""))
			else:
				file.write(textwrap.dedent(f"""
					{T1}
					0
					2
					0 0 0 0 0 0
					{i_dis} 0 0 0 0 0
					"""))
		else:
			file.write(textwrap.dedent(f"""
				{T1}
				0
				1
				0 0 0 0 0 0\n
				{T2}
				0
				1
				{i_dis} 0 0 0 0 0
				"""))

	with open(folder/"KPT","w") as file:
		file.write(textwrap.dedent(f"""\
			K_POINTS
			0
			Gamma
			1 1 1 0 0 0
			"""))

	info = utils.read_info()
	if utils.sub=="qsub":
		with open(folder/"sub.sh","w") as file:
			file.write(textwrap.dedent(f"""\
				#!/bin/bash
				#PBS -q gold5120
				#PBS -l nodes=1:ppn=1
				#PBS -l walltime=1:00:00
				#PBS -o job.log
				#PBS -e job.err
				ulimit -s unlimited
				cd $PBS_O_WORKDIR
				EXEC={info["ABACUS"]}
				mpirun -n 1 -env OMP_NUM_THREADS=1 $EXEC
				"""))
	elif utils.sub=="bsub":
		with open(folder/"sub.sh","w") as file:
			file.write(textwrap.dedent(f"""\
				#!/bin/sh
				#BSUB -q renxg
				#BSUB -o job.log -e job.err
				#BSUB -n 1
				EXEC={info["ABACUS"]}
				mpirun -n 1 -env OMP_NUM_THREADS=1 $EXEC
				"""))

	os.chdir(folder)
	if utils.sub=="qsub":
		os.system("qsub sub.sh")
	elif utils.sub=="bsub":
		os.system("bsub < sub.sh")		
	elif utils.sub=="tianh2":
		os.system(f'yhrun -n 1 -c 1 {info["ABACUS"]} >Log.txt')
	os.chdir("../")

# dis_opt[T1,T2] = [..., dis, ...]		
def get_dis_opt(dis):
	opt_mode = "kmeans"
	dis_opt = dict()
	info = utils.read_info()
	for T1,T2 in dis:
		dis_tmp = read_stru.delete_zero(dis[T1,T2])
		if len(dis_tmp)<=info["dimer_num"]:
			dis_opt[T1,T2] = list(dis_tmp.keys())
		else:
			if opt_mode=="linspace":
				dis_opt[T1,T2] = np.linspace( min(dis_tmp), max(dis_tmp), info["dimer_num"] )
			elif opt_mode=="kmeans":
				kmeans = KMeans(n_clusters=info["dimer_num"])
				label = kmeans.fit_predict(
					np.array(list(dis_tmp.keys())).reshape(-1,1), 
					sample_weight = [num/i_dis**2 for i_dis,num in dis_tmp.items()])
				dis_opt[T1,T2] = list(kmeans.cluster_centers_.reshape(-1))
				pprint.pprint(dict(zip(dis_tmp.keys(),label)))
		if T1==T2:
			dis_opt[T1,T2].append(0.0)
	return dis_opt
	
# dis_weight[T1,T2] = {..., i_dis_opt:weight, ...}
def cal_dis_weight(dis_opt,dis_all):
	def weight_func(x,D):
		return (2/D**3)*x**3 + (-3/D**2)*x**2 + 1
	dis_weight = collections.defaultdict(dict)
	for T1,T2 in dis_opt:
		dis_opt_TT = sorted(dis_opt[T1,T2])
		for index,i_dis_opt in enumerate(dis_opt_TT):
			i_weight = 0
			i_dis_low = dis_opt_TT[index-1] if index>0 else -np.infty
			i_dis_up = dis_opt_TT[index+1] if index<len(dis_opt_TT)-1 else np.infty
			for i_dis,num in dis_all[T1,T2].items():
				if i_dis_low<i_dis<i_dis_opt:
					i_weight += weight_func( i_dis_opt-i_dis, i_dis_opt-i_dis_low ) * num
				elif i_dis==i_dis_opt:
					i_weight += num
				elif i_dis_opt<i_dis<i_dis_up:
					i_weight += weight_func( i_dis_opt-i_dis, i_dis_opt-i_dis_up ) * num
			dis_weight[T1,T2][i_dis_opt] = i_weight
	return dis_weight
	
def cal():
	dis = read_stru.cut_dis(read_stru.cal_dis(read_stru.change_R(read_stru.get_R())))
	dis_decimal = read_stru.round_dis(dis,1E-6)
	pprint.pprint(dis_decimal)
	dis_opt = get_dis_opt(dis_decimal)
	for T1,T2 in dis_opt:
		for i_dis in dis_opt[T1,T2]:
			cal_ABACUS(T1,T2,i_dis)
	dis_weight = cal_dis_weight(dis_opt,dis)
	read_stru.print_folders(dis_weight)
			
if __name__=="__main__":
	cal()