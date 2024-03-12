import json
import os
import re
import textwrap
import pathlib
import read_stru
import utils

def cp_matrix(folders_matrix):
	pathlib.Path(utils.folder_opt+"/"+utils.folder_opt_matrix).mkdir(parents=True,exist_ok=False)
	folders_same_atom = set()
	for folder_matrix in folders_matrix:
		T1,T2 = folder_matrix.split("_")[0].split("-")
		dis = float(folder_matrix.split("_")[1])
		if dis:
			if T1==T2:
				matrix_file = "matrix_0_0_0_1" 
			else:
				matrix_file = "matrix_0_0_1_0"
		else:
			matrix_file = "matrix_0_0_0_0"
		os.system(f"cp {folder_matrix}/{matrix_file} {utils.folder_opt}/{utils.folder_opt_matrix}/{folder_matrix}")
					
def add_nbands(folders_matrix):
	for folder_matrix in folders_matrix:
		with open(f"{utils.folder_opt}/{utils.folder_opt_matrix}/{folder_matrix}","r") as file:
			nband = int(re.compile(r"(\d+)\s+nbands").search(file.read()).group(1))
			folders_matrix[folder_matrix] = (folders_matrix[folder_matrix], nband)

def set_input(folders_matrix):
	Ts = read_stru.get_T()
	input_dict = read_stru.get_input_dict()
	Ecut = float(input_dict["ecutwfc"])
	info = utils.read_info()
	input = {
		"file_list": [f"{utils.folder_opt_matrix}/"+folder_matrix for folder_matrix in folders_matrix],
		"info":
		{
			"Nt_all":	Ts,
			"Nu":		{T:info["Nu"] for T in Ts},
			"Nb_true":	[nbands for weight,nbands in folders_matrix.values()],
			"weight":	[weight for weight,nbands in folders_matrix.values()],
			"Rcut":		read_stru.get_Rcut(),
			"dr":		{T:utils.dr for T in Ts},
			"Ecut":		{T:Ecut for T in Ts},
			"lr":		utils.lr
		},
		"C_init_info":{ "init_from_file": False },
		"V_info":
		{
			"init_from_file":	True,
			"same_band":		False
		}
	}
	return json.dumps(input,indent=4)
	
def cal(input):
	info = utils.read_info()
	with open(f"{utils.folder_opt}/input.json","w") as file:
		file.write(input)
	if utils.sub=="qsub":
		with open(f"{utils.folder_opt}/sub.sh","w") as file:
			file.write(textwrap.dedent(f"""\
				#!/bin/bash
				#PBS -q gold5120
				#PBS -l nodes=1:ppn=1
				#PBS -l walltime=2:00:00
				#PBS -o job.log
				#PBS -e job.err
				ulimit -s unlimited
				cd $PBS_O_WORKDIR
				export OMP_NUM_THREADS=1
				EXEC={info["opt_orb"]}
				python3 -u $EXEC
				"""))
	elif utils.sub=="tianhe2":
		with open(f"{utils.folder_opt}/sub.sh","w") as file:
			file.write(textwrap.dedent(f"""\
				#!/bin/bash
				EXEC={info["opt_orb"]}
				python3 -u $EXEC >Log.txt
				"""))
	os.chdir(utils.folder_opt)
	if utils.sub=="qsub":
		os.system("qsub sub.sh")
	elif utils.sub=="tianhe2":
		os.system("yhbatch -N 1 sub.sh")
#	os.system(f'python3 -u {info["opt_orb"]}')
	os.chdir("../")
		
def all():
	pathlib.Path(utils.folder_opt).mkdir(parents=True,exist_ok=False)
	with open("folders","r") as file:
		folders_matrix = json.loads(file.read())
	cp_matrix(folders_matrix.keys())
	add_nbands(folders_matrix)
	input = set_input(folders_matrix)
	cal(input)
	
if __name__=="__main__":
	all()