import json
import os

def read_info(file_name):
	with open(file_name,"r") as file:
		info = file.read()
	info = json.loads(info)
	info["input"]["pseudo_dir"] = os.path.abspath(info["input"].get("pseudo_dir","./"))
	info["exe"]["qsub"] = info["exe"].get("qsub",[1,1])
	return info

""" info.json
{
	"exe":
	{
		"exe_pw":		"/home/linpz/ABACUS/ABACUS.1.0.1_2016-12-19/bin/ABACUS.mpi.1.0.0",
		"exe_orbital":	"/home/linpz/ABACUS/ABACUS_RI_git/tools/opt_orb_pytorch/main.py",
		"qsub":			[1,8]
	},
	"input":
	{
		"element":		"N",
		"nbands":		8,
		"ecut":			50,
		"rcut":			6,
		"pseudo_dir":	"../",
		"pseudo":		"N_ONCV_PBE-1.0.upf",
		"smearing_sigma":		0.01,
		"smooth":		false,
	},
	"orbital":	[2,2,1]
}
"""

