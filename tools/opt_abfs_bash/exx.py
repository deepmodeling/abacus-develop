import os
import re
import textwrap
import pathlib
import utils
import read_stru

def cal():
	pathlib.Path(utils.folder_exx).mkdir(parents=True,exist_ok=False)
	
	os.system(f"cp INPUT {utils.folder_exx}/")
	os.system(f"cp KPT {utils.folder_exx}/")
	
	with open(f"{utils.folder_exx}/INPUT","w") as file:
		input_dict = read_stru.get_input_dict()
		input_dict["pseudo_dir"] = os.path.abspath(input_dict.get("pseudo_dir",r"./"))
		read_stru.print_input(file,input_dict,1)
	
	with open("STRU","r") as file:
		strus = re.compile("LATTICE_CONSTANT").split(file.read())
	with open(f"{utils.folder_exx}/STRU","w") as file:
		Ts = read_stru.get_T()
		file.write("ATOMIC_SPECIES\n")
		pseudo_path = read_stru.get_pseudo_path()
		for T in Ts:
			file.write(f"{T}	12	{pseudo_path[T]}\n")
		file.write("\nNUMERICAL_ORBITAL\n")
		lcao_path = read_stru.get_lcao_path()
		for T in Ts:
			file.write(f"{lcao_path[T]}\n")
		file.write("\nABFS_ORBITAL\n")
		for T in read_stru.get_T():
			file.write(f"../{utils.folder_opt}/orb_{T}.dat\n")
		file.write("\nLATTICE_CONSTANT")
		file.write(strus[1])
	
	info = utils.read_info()
	if utils.sub=="qsub":
		with open(f"{utils.folder_exx}/sub.sh","w") as file:
			file.write(textwrap.dedent(f"""\
				#!/bin/bash
				#PBS -q gold5120
				#PBS -l nodes=2:ppn=28
				#PBS -l walltime=99:99:99
				#PBS -o job.log
				#PBS -e job.err
				ulimit -s unlimited
				cd $PBS_O_WORKDIR
				EXEC={info["ABACUS"]}
				mpirun -n 2 -env OMP_NUM_THREADS=28 $EXEC
				"""))
	elif utils.sub=="bsub":
		with open(f"{utils.folder_exx}/sub.sh","w") as file:
			file.write(textwrap.dedent(f"""\
				#!/bin/sh
				#BSUB -q renxg
				#BSUB -o job.log -e job.err
				#BSUB -n 6
				mpirun -n 2 -env OMP_NUM_THREADS=28 {info['ABACUS']}
				"""))

	os.chdir(utils.folder_exx)
	if utils.sub=="qsub":
		os.system("qsub sub.sh")
	elif utils.sub=="bsub":
		os.system(f"bsub < sub.sh")		
	elif utils.sub=="tianhe2":
		os.system(f'yhrun -N 1 -n 1 -c 24 -t 1440 {info["ABACUS"]} >Log.txt 2>&1 &')
	else:
		raise ValueError("utils.sub")
	os.chdir("../")
	
if __name__=="__main__":
	cal()
