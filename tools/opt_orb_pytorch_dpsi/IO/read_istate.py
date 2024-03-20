import re
import torch
import itertools

# occ[ik][ib]
def read_istate(file_name):
	nspin0 = get_nspin0(file_name)
	if nspin0==1:	occ = [[]]
	elif nspin0==2:	occ = [[],[]]
	with open(file_name,"r") as file:
		content = file.read().split("BAND")
		for content_k in content[1:]:
			content_k = content_k.split("\n")
			k = get_k(content_k[0])
			for ispin in range(nspin0):
				occ[ispin].append([])
			for line in content_k[1:]:
				line = line.strip()
				if line:
					line = line.split()
					if nspin0==1:
						occ[0][-1].append(float(line[2]))
					elif nspin0==2:
						occ[0][-1].append(float(line[2]))
						occ[1][-1].append(float(line[4]))
			for ispin in range(nspin0):
				occ[ispin][-1] = torch.Tensor(occ[ispin][-1])
	occ = list(itertools.chain(*occ))
	return occ

def get_k(line):
	k = re.compile(r"Kpoint\s*=\s*(\d+)").search(line).group(1)
	return int(k)

def get_nspin0(file_name):
	with open(file_name,"r") as file:
		file.readline()
		line = file.readline()
		lens = len(line.split())
		if lens == 3:	return 1
		elif lens == 5:	return 2
		else:	raise