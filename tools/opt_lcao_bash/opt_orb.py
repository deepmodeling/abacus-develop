import os
import sys
import pathlib
from read_info import read_info
from print_file import print_file_pw, print_file_opt
from dis import dis
import utils

def cal_pw():
	info = read_info("info.json")
	for distance in dis[info["input"]["element"]]:
		folder = f'{info["input"]["element"]}-{info["input"]["rcut"]}-{distance}'
		pathlib.Path(folder).mkdir(parents=True,exist_ok=False)
		os.chdir(folder)
		print_file_pw(info,distance)
		if utils.sub in ["qsub", "sbatch"]:
			os.system(f"{utils.sub} sub.sh")
		elif utils.sub=="bsub":
			os.system("bsub < sub.sh")
		else:
			raise KeyError("utils.sub = ",utils.sub)		
		os.chdir("../")
		
def cal_opt():
	info = read_info("info.json")
	pathlib.Path("opt_orb").mkdir(parents=True,exist_ok=False)
	os.chdir("opt_orb")
	print_file_opt(info,dis)
	if utils.sub in ["qsub", "sbatch"]:
		os.system(f"{utils.sub} sub.sh")
	elif utils.sub=="bsub":
		os.system("bsub < sub.sh")
	else:
		raise KeyError("utils.sub = ",utils.sub)
	os.chdir("../")
		
if __name__=="__main__":
	if sys.argv[1]=='1':
		cal_pw()
	elif sys.argv[1]=='2':
		cal_opt()