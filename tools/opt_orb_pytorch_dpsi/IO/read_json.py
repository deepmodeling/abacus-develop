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
		
	return input["file_list"], info, input["weight"], input["C_init_info"], input["V_info"]

	""" file_name
	{
		"file_list": {
			"origin": [
				"~/C_bulk/orb_matrix/test.0.dat",
				"~/CO2/orb_matrix/test.0.dat"
			],
			"linear": [
				[
					"~/C_bulk/orb_matrix/test.1.dat",
					"~/CO2/orb_matrix/test.1.dat"
				],
				[
					"~/C_bulk/orb_matrix/test.2.dat",
					"~/CO2/orb_matrix/test.2.dat"
				],
			]
		},
		"info": {
			"Nt_all": [ "C", "O" ],
			"Nu":   { "C":[2,2,1], "O":[3,2,1] },
			"Rcut": { "C":6,       "O":6       },
			"dr":   { "C":0.01,    "O":0.01    },
			"Ecut": { "C":200,     "O":200     },
			"lr": 0.01,
			"cal_T": false,
			"cal_smooth": false
		},
		"weight":
		{
			"stru":	[1, 2.3],
			"bands_range":	[10, 15],					# "bands_range" and "bands_file" only once
			"bands_file":
			[
				"~/C_bulk/OUT.ABACUS/istate.info",
				"~/CO2/OUT.ABACUS/istate.info"
			]
		},	
		"C_init_info": {
			"init_from_file": false,
			"C_init_file": "~/CO/ORBITAL_RESULTS.txt",
			"opt_C_read": false
		},
		"V_info": {
			"init_from_file": true,
			"same_band": true
		}
	}
	"""
	
	""" info
		Nt_all		['C', 'O']
		Nu			{'C': [2, 2, 1], 'O': [3, 2, 1]}
		Rcut		{'C': 6, 'O': 6}
		dr			{'C': 0.01, 'O': 0.01}
		Ecut		{'C': 200, 'O': 200}
		lr			0.01
		cal_T		False
		cal_smooth	False
		Nl			{'C': 3, 'O': 3}
		Nst			3
		Nt			[['C'], ['C'], ['C', 'O']]
		Na			[{'C': 1}, {'C': 1}, {'C': 1, 'O': 2}]
		Nb			[6, 6, 10]
		Ne			{'C': 19, 'O': 19}
	"""