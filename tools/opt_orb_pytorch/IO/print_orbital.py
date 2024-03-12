def print_orbital(orb,info):
	""" orb[it][il][iu][r] """
	for it,orb_t in orb.items():
		with open("orb_{0}.dat".format(it),"w") as file:
			print_orbital_head(file,info,it)
			for il,orb_tl in enumerate(orb_t):
				for iu,orb_tlu in enumerate(orb_tl):
					print("""                Type                   L                   N""",file=file)
					print("""                   0                   {0}                   {1}""".format(il,iu),file=file)
					for ir,orb_tlur in enumerate(orb_tlu):
						print(orb_tlur,end="\t",file=file)
						if ir%4==3:	print(file=file)
					print(file=file)
					
					
def plot_orbital(orb,Rcut,dr):
	for it,orb_t in orb.items():
		with open("orb_{0}_plot.dat".format(it),"w") as file:
			Nr = int(Rcut[it]/dr[it])+1
			for ir in range(Nr):
				print(ir*dr[it],end="\t",file=file)
				for il,orb_tl in enumerate(orb_t):
					for orb_tlu in orb_tl:
						print(orb_tlu[ir],end="\t",file=file)
				print(file=file)
				
				
def print_orbital_head(file,info,it):
	print( "---------------------------------------------------------------------------", file=file )
	print( "Element                     {0}".format(it), file=file )
	print( "Energy Cutoff(Ry)           {0}".format(info.Ecut[it]), file=file )
	print( "Radius Cutoff(a.u.)         {0}".format(info.Rcut[it]), file=file )
	print( "Lmax                        {0}".format(info.Nl[it]-1), file=file )
	l_name = ["S","P","D"]+list(map(chr,range(ord('F'),ord('Z')+1)))
	for il,iu in enumerate(info.Nu[it]):
		print( "Number of {0}orbital-->       {1}".format(l_name[il],iu), file=file )
	print( "---------------------------------------------------------------------------", file=file )
	print( "SUMMARY  END", file=file )
	print( file=file )
	print( "Mesh                        {0}".format(int(info.Rcut[it]/info.dr[it])+1), file=file )
	print( "dr                          {0}".format(info.dr[it]), file=file )

	