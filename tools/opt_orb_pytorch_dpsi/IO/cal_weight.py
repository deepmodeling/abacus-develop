import IO.read_istate
import torch
import re
import functools
import operator

def cal_weight(info_weight, flag_same_band, stru_file_list=None):
	""" weight[ist][ib] """

	if "bands_file" in info_weight.keys():
		if "bands_range" in info_weight.keys():
			raise IOError('"bands_file" and "bands_range" only once')

		weight = []														# weight[ist][ib]
		for weight_stru, file_name in zip(info_weight["stru"], info_weight["bands_file"]):
			occ = IO.read_istate.read_istate(file_name)
			weight += [occ_k * weight_stru for occ_k in occ]

	elif "bands_range" in info_weight.keys():
		k_weight = read_k_weight(stru_file_list)						# k_weight[ist][ik]
		nbands = read_nbands(stru_file_list)							# nbands[ist]

		st_weight = []													# st_weight[ist][ib]
		for weight_stru, bands_range, nbands_ist in zip(info_weight["stru"], info_weight["bands_range"], nbands):
			st_weight_tmp = torch.zeros((nbands_ist,))
			st_weight_tmp[:bands_range] = weight_stru
			st_weight.append( st_weight_tmp )

		weight = []														# weight[ist][ib]
		for ist,_ in enumerate(k_weight):
			for ik,_ in enumerate(k_weight[ist]):
				weight.append(st_weight[ist] * k_weight[ist][ik])

	else:
		raise IOError('"bands_file" and "bands_range" must once')


	if not flag_same_band:
		for ist,_ in enumerate(weight):
			weight[ist] = torch.tensordot(weight[ist], weight[ist], dims=0)


	normalization = functools.reduce(operator.add, map(torch.sum, weight), 0)
	weight = list(map(lambda x:x/normalization, weight))

	return weight


def read_k_weight(stru_file_list):
	""" weight[ist][ik] """
	weight = []												# weight[ist][ik]
	for file_name in stru_file_list:
		weight_k = []										# weight_k[ik]
		with open(file_name,"r") as file:
			data = re.compile(r"<WEIGHT_OF_KPOINTS>(.+)</WEIGHT_OF_KPOINTS>", re.S).search(file.read()).group(1).split("\n")
			for line in data:
				line = line.strip()
				if line:
					weight_k.append(float(line.split()[-1]))
		weight.append(weight_k)
	return weight


def read_nbands(stru_file_list):
	""" nbands[ib] """
	nbands = []
	for file_name in stru_file_list:
		with open(file_name,"r") as file:
			nbands.append(int(re.compile(r"(\d+)\s+nbands").search(file.read()).group(1)))
	return nbands
