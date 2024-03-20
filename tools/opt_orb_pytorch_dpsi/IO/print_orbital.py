periodtable = {   'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7,
                  'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13,
               'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
               'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,
               'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31,
               'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
               'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
               'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49,
               'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
               'Ba': 56, #'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61,
                    ## 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67,
                    ## 'Er': 68, 'Tm': 69, 'Yb': 70, 
                    ## 'Lu': 71, 
               'Hf': 72, 'Ta': 73,
               'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
               'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 
                    ## 'Po': 84, #'At': 85,
                    ## 'Rn': 86, #'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91,
                    ## 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97,
                    ## 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103,
                    ## 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108,
                    ## 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Uut': 113,
                    ## 'Fl': 114, 'Uup': 115, 'Lv': 116, 'Uus': 117, 'Uuo': 118
               } 

def print_orbital(orb,info_element):
	""" orb[it][il][iu][r] """
	for it,orb_t in orb.items():
		#with open("orb_{0}.dat".format(it),"w") as file:
		with open("ORBITAL_{0}U.dat".format( periodtable[it] ),"w") as file:
			print_orbital_head(file,info_element,it)
			for il,orb_tl in enumerate(orb_t):
				for iu,orb_tlu in enumerate(orb_tl):
					print("""                Type                   L                   N""",file=file)
					print("""                   0                   {0}                   {1}""".format(il,iu),file=file)
					for ir,orb_tlur in enumerate(orb_tlu):
						print( '%.14e'%orb_tlur, end="  ",file=file)
						if ir%4==3:	print(file=file)
					print(file=file)
					
					
def plot_orbital(orb,Rcut,dr):
	for it,orb_t in orb.items():
		#with open("orb_{0}_plot.dat".format(it),"w") as file:
		with open("ORBITAL_PLOTU.dat", "w") as file:
			Nr = int(Rcut[it]/dr[it])+1
			for ir in range(Nr):
				print( '%10.6f'%(ir*dr[it]),end="  ",file=file)
				for il,orb_tl in enumerate(orb_t):
					for orb_tlu in orb_tl:
						print( '%18.14f'%orb_tlu[ir],end="  ",file=file)
				print(file=file)
				
				
def print_orbital_head(file,info_element,it):
	print( "---------------------------------------------------------------------------", file=file )
	print( "Element                     {0}".format(it), file=file )
	print( "Energy Cutoff(Ry)           {0}".format(info_element[it].Ecut), file=file )
	print( "Radius Cutoff(a.u.)         {0}".format(info_element[it].Rcut), file=file )
	print( "Lmax                        {0}".format(info_element[it].Nl-1), file=file )
	l_name = ["S","P","D"]+list(map(chr,range(ord('F'),ord('Z')+1)))
	for il,iu in enumerate(info_element[it].Nu):
		print( "Number of {0}orbital-->       {1}".format(l_name[il],iu), file=file )
	print( "---------------------------------------------------------------------------", file=file )
	print( "SUMMARY  END", file=file )
	print( file=file )
	print( "Mesh                        {0}".format(int(info_element[it].Rcut/info_element[it].dr)+1), file=file )
	print( "dr                          {0}".format(info_element[it].dr), file=file )

	
