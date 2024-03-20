import json
import functools

@functools.lru_cache(maxsize=None)
def read_info(flag=True):
	"""
	{
		"ABACUS":"dir/ABACUS/ABACUS.mpi.1.0.0",
		"opt_orb":"dir/opt_orb/main.py",
		"Nu":[7,10,8,3,1],
		"dimer_num":5
	}
	"""
	with open("info","r") as file:
		return json.loads(file.read())

def ignore_lines(file,n):
	for i in range(n):
		file.readline()

def search_sentence(file,sentence):
	sentence = sentence.strip()
	for line in file:
		if line.strip() == sentence:
			break
			
def floatlist_str(l):
	return "-".join(list(map(str,l)))

sub = "qsub"
dr  = 0.01
lr  = 0.01

def folder_name(T1,T2,i_dis):
	return f"{T1}-{T2}_{i_dis}"
try:
	info = read_info(False)
	folder_opt = "opt_orb_"+floatlist_str(info["Nu"])
	folder_exx = "exx_"+floatlist_str(info["Nu"])
	folder_opt_matrix = "matrix"
except FileNotFoundError:
	print("info FileNotFoundError")