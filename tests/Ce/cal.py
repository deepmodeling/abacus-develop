import os
import textwrap

def cal(nspin, xc):
	folder = f"Ce_nspin{nspin}_{xc}"
	print(folder)
	os.mkdir(folder)
	os.chdir(folder)

	with open("INPUT","w") as file:
		file.write(textwrap.dedent(f"""\
			INPUT_PARAMETERS
			#Parameters (1.General)
			symmetry		1
			gamma_only      0
			nspin			{nspin}

			#Parameters (2.Iteration)
			ecutwfc			20
			scf_thr			1E-6
			scf_nmax		100

			calculation		scf
			cal_stress		1

			#Parameters (3.Basis)
			basis_type		pw

			#Parameters (4.Smearing)
			smearing_method		gauss
			smearing_sigma		0.002

			#Parameters (5.Mixing)
			mixing_type		pulay

			#Parameters (7.Hybrid)
			dft_functional				{"PBE" if xc=="xc" else "GGA_X_PBE+GGA_C_PBE"}

			pseudo_dir                     ../element
			orbital_dir                    ../element
			"""))
	
	
	with open("STRU","w") as file:
		file.write(textwrap.dedent(f"""\
			ATOMIC_SPECIES
			Ce 140.115 58_Ce.UPF upf201

			LATTICE_CONSTANT
			8.92

			LATTICE_VECTORS
			0.5 0.5 0.0
			0.5 0.0 0.5
			0.0 0.5 0.5

			ATOMIC_POSITIONS
			Direct

			Ce
			0.0
			1
			0.00 0.00 0.00 1 1 1
			"""))
	
	
	with open("KPT","w") as file:
		file.write(textwrap.dedent(f"""\
			K_POINTS
			0
			Gamma
			2 2 2 0 0 0
			"""))
	

	os.system("mpirun -n 20 ../../../ABACUS/abacus-develop/bin/abacus >job.log 2>job.err")
	os.chdir("../")


for nspin in [1,2,4]:
	for xc in ["xc","libxc"]:
		cal(nspin, xc)