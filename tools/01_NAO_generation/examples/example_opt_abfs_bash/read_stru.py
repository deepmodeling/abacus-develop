import re
import numpy as np
import itertools
import numba
import functools
import operator
import os
import collections
import decimal
import json
import utils

def skip_notes(line):
	line = re.compile(r"#.*").sub("",line)
	line = re.compile(r"//.*").sub("",line)
	line = line.strip()
	return line

@functools.lru_cache(maxsize=None)
def get_input_dict():
	input_dict = dict()
	with open("INPUT","r") as file:
		for line in file:
			line = skip_notes(line)
			if not line or line=="INPUT_PARAMETERS":
				continue
			label, value = re.compile(r"\s+").split( line, maxsplit=1 )
			input_dict[label] = value
	return input_dict

def print_input(file,input_dict):
	print("INPUT_PARAMETERS",file=file)
	max_len = max([len(label) for label in input_dict])			
	for label,value in input_dict.items():
		print( label, " "*(max_len-len(label)+4), value, file=file)
	
	
"""
@functools.lru_cache(maxsize=None)
def search_in_input(label):
	with open("INPUT","r") as file:
		regex = label+"\s+(\d+)"
		value = re.compile(regex).search(file.read()).group(1)
	return value
"""

@functools.lru_cache(maxsize=None)
def get_k():
	with open("KPT","r") as file:
		utils.search_sentence(file,"Gamma")
		return list(map(int,skip_notes(file.readline()).split()[:3]))

@functools.lru_cache(maxsize=None)
def get_T():
	ntype = int(get_input_dict()["ntype"])
	T = []
	with open("STRU","r") as file:
		utils.search_sentence(file,"ATOMIC_SPECIES")
		for it in range(ntype):
			line = skip_notes(file.readline())
			T.append(line.split()[0])
	return T

@functools.lru_cache(maxsize=None)
def get_Rcut():
	Ts = get_T()
	lcao_path = get_lcao_path()
	Rcut = dict()
	for T in Ts:
		with open(lcao_path[T],"r") as file:
			for line in file:
				line = skip_notes(line)
				if line.startswith("Radius Cutoff(a.u.)"):
					Rcut[T] = float(line.split()[-1])
					break
	return Rcut

@functools.lru_cache(maxsize=None)
def get_pseudo_path():
	Ts = get_T()
	path = dict()
	with open("STRU","r") as file:
		utils.search_sentence(file,"ATOMIC_SPECIES")
		for T in Ts:
			path[T] = skip_notes(file.readline()).split()[-1]
	return path
	
@functools.lru_cache(maxsize=None)
def get_lcao_path():
	Ts = get_T()
	path = dict()
	with open("STRU","r") as file:
		utils.search_sentence(file,"NUMERICAL_ORBITAL")
		for T in Ts:
			path[T] = os.path.abspath(skip_notes(file.readline()))
	return path
	
@functools.lru_cache(maxsize=None)
def get_lattice():
	with open("STRU","r") as file:
		utils.search_sentence(file,"LATTICE_CONSTANT")
		lat0 = float(skip_notes(file.readline()).split()[0])
	with open("STRU","r") as file:
		utils.search_sentence(file,"LATTICE_VECTORS")
		lat_vec = []
		for i in range(3):
			lat_vec.append(list(map(float,skip_notes(file.readline()).split())))
		lat_vec = np.array(lat_vec)
	with open("STRU","r") as file:
		utils.search_sentence(file,"ATOMIC_POSITIONS")
		position = skip_notes(file.readline())
	return lat0, lat_vec, position
	
@functools.lru_cache(maxsize=None)
def get_lcao():
	lcao = collections.defaultdict(list)
	lcao_path = get_lcao_path()
	for T in lcao_path:
		with open(lcao_path[T],"r") as file:
			for line in file:
				line = skip_notes(line)
				if line.startswith("Number of"):
					lcao[T].append(int(line.split()[-1]))
	return dict(lcao)
	
@functools.lru_cache(maxsize=None)
def get_nw():
	nw = dict()
	lcao = get_lcao()
	for T in lcao:
		nw[T] = functools.reduce(operator.add,((2*l+1)*n for l,n in enumerate(lcao[T])))
	return nw

# R[T] = [..., [xi,yi,zi], ...]
def get_R():
	Ts = get_T()
	R = dict()
	with open("STRU","r") as file:
		for T in Ts:
			utils.search_sentence(file,T)
			utils.ignore_lines(file,1)
			na = int(skip_notes(file.readline()).split()[0])
			R_tmp = []
			for i in range(na):
				R_tmp.append(list(map(float,skip_notes(file.readline()).split()[:3])))
			R[T] = np.array(R_tmp)
	return R

def change_R(R):
	lat0, lat_vec, position = get_lattice()
	if position == "Direct":
		for T in R:
			R[T] = np.dot(R[T],lat_vec)
	elif position == "Cartesian_angstrom":
		for T in R:
			R[T] /= 0.529166
	for T in R:
		R[T] *= lat0

	nx,ny,nz = get_k()
	lat_vec *= lat0
	for T in R:
		R_new = []
		for ix in range(nx):
			for iy in range(ny):
				for iz in range(nz):
					R_new.append( R[T] + np.dot(np.array([ix,iy,iz]), lat_vec) )
		R[T] = np.concatenate(R_new)

	return R

# dis[T1,T2] = {..., i_dis:num, ...}
def cal_dis(R):
	dis = dict()
	for T1,T2 in itertools.combinations_with_replacement(R,2):
		dis_TT = collections.defaultdict(int)
		for ia1,ia2 in itertools.product(R[T1],R[T2]):
			i_dis = np.linalg.norm(ia1-ia2)
			dis_TT[i_dis] += 1
		dis[T1,T2] = dict(dis_TT)
	return dis
	
def round_dis(dis,precision):
	dis_round = dict()
	for T1,T2 in dis:
		dis_TT = collections.defaultdict(int)
		for i_dis,num in dis[T1,T2].items():
			i_dis = float(decimal.Decimal(i_dis).quantize(decimal.Decimal(str(precision)), rounding=decimal.ROUND_HALF_UP))		# test
			dis_TT[i_dis] += num
		dis_round[T1,T2] = dict(dis_TT)
	return dis_round

def cut_dis(dis):
	Rcut = get_Rcut()
	for T1,T2 in dis:
		Rcut_sum = Rcut[T1]+Rcut[T2]
		dis[T1,T2] = { i_dis:num for i_dis,num in dis[T1,T2].items() if i_dis<Rcut_sum }
	return dis
	
"""
def cal_cut_dis(R):
	Rcut = get_Rcut()
	dis = dict()
	for T1,T2 in itertools.combinations_with_replacement(R,2):
		Rcut_sum = Rcut[T1]+Rcut[T2]
		dis[T1,T2] = set()
		for ia1,ia2 in itertools.product(R[T1],R[T2]):
			delta = abs(ia1-ia2)
			if not any(delta>Rcut_sum):
				norm = np.linalg.norm(delta)
				if norm>0 and norm<Rcut_sum:
					dis[T1,T2].add(norm)
	return dis
"""
	
def delete_zero(dis):
	dis_new = dis.copy()
	while 0 in dis_new:
		if type(dis_new) in [list,set]:
			dis_new.remove(0)
		elif type(dis_new) in [dict]:
			dis_new.pop(0)
		else:
			raise TypeError
	return dis_new
	
def print_folders(dis_weight):
	folders = dict()
	for T1,T2 in dis_weight:
		for i_dis,i_weight in dis_weight[T1,T2].items():
			folders[utils.folder_name(T1,T2,i_dis)] = i_weight
	with open("folders","w") as file:
		print(json.dumps(folders,indent=4), file=file)