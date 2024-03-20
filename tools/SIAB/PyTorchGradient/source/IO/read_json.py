import json
from util import Info

def read_json(file_name):

	with open(file_name,"r") as file:
		input = file.read()
	input = json.loads(input)
	
	info = Info()
	for info_attr,info_value in input["info"].items():
		info.__dict__[info_attr] = info_value
	info.Nl = { it:len(Nu) for it,Nu in info.Nu.items() }
		
	return input["file_list"], info, input["C_init_info"], input["V_info"]

	""" file_name
	{
		"file_list": 
		[
			"/public/udata/linpz/try/SIA/pytorch/test/many_atoms/S2/S2.10.dat",
			"/public/udata/linpz/try/SIA/pytorch/test/many_atoms/SO2/SO2.10.dat",
			"/public/udata/linpz/try/SIA/pytorch/test/many_atoms/SO/SO.10.dat"
		],
		"info":
		{
			"Nt_all":	["S","O"],
			"Nu":		{"S":[3,3,2],"O":[3,3,2]},
			"Nb_true":	[7,9,7],
			"weight":	[2,3,1.5],
			"Rcut":		{"S":10,"O":10},
			"dr":		{"S":0.01,"O":0.01},
			"Ecut":		{"S":100,"O":100},
			"lr":		0.01
		},
		"C_init_info":
		{
			"init_from_file":	false,
			"C_init_file":		"/public/udata/linpz/try/SIA/pytorch/test/many_atoms/SIA/ORBITAL_RESULTS.txt"
		},
		"V_info":
		{
			"init_from_file":	true,
			"same_band":		false
		}
	}
	"""
	
	""" info
	Nt_all	["S", "O"]
	Nu		{"S":[3,3,2], "O":[3,3,2]}
	Nb_true	[7, 9, 7]
	weight	[2, 3, 1.5]
	Rcut	{"S":10, "O":10}
	dr		{"S":0.01, "O":0.01}
	Ecut	{"S":100, "O":100}
	lr		0.01
	Nl		{"S":2, "O":2}
	Nst		3
	Nt		[["S"], ["S","O"], ["S","O"]]
	Na		[{"S":2}, {"S":1,"O":2}, {"S":1,"O":1}]
	Nb		[7, 9, 7]
	Ne		{"S":22, "O":19}
	"""