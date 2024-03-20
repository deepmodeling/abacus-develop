from util import ND_list
import numpy as np
from scipy.special import spherical_jn
from scipy.integrate import simps
from scipy.optimize import fsolve
import functools
import torch

def generate_orbital(info_element,C,E):
	""" C[it][il][ie,iu] """
	""" orb[it][il][iu][r] = \suml_{ie} C[it][il][ie,iu] * jn(il,ie*r) """
	orb = dict()
	for it in info_element:
		Nr = int(info_element[it].Rcut/info_element[it].dr)+1
		orb[it] = ND_list(info_element[it].Nl)
		for il in range(info_element[it].Nl):
			orb[it][il] = ND_list(info_element[it].Nu[il])
			for iu in range(info_element[it].Nu[il]):
				orb[it][il][iu] = np.zeros(Nr)
				for ir in range(Nr):
					r = ir * info_element[it].dr
					for ie in range(info_element[it].Ne):
						orb[it][il][iu][ir] += C[it][il][ie,iu].item() * spherical_jn(il,E[it][il,ie].item()*r)
	return orb
						
						
def smooth_orbital(orb,Rcut,dr,smearing_sigma):
	for it,orb_t in orb.items():
		for orb_tl in orb_t:
			for orb_tlu in orb_tl:
				for ir in range(orb_tlu.shape[0]):
					assert orb_tlu.shape[0] == int(Rcut[it]/dr[it])+1
					r = ir * dr[it]
					orb_tlu[ir] *= 1-np.exp( -(r-Rcut[it])**2/(2*smearing_sigma**2) )
	
	
	
def inner_product( orb1, orb2, dr ):
	assert orb1.shape == orb2.shape
	r = np.array(range(orb1.shape[0]))*dr
	return simps( orb1 * orb2 * r * r, dx=dr )
	
def normalize(orb,dr,C=None,flag_norm_orb=False,flag_norm_C=False):
	""" C[it][il][ie,iu] """
	""" orb[it][il][iu][r] = \suml_{ie} C[it][il][ie,iu] * jn(il,ie*r) """
	for it,orb_t in orb.items():
		for il,orb_tl in enumerate(orb_t):
			for iu,orb_tlu in enumerate(orb_tl):
				norm = np.sqrt(inner_product(orb_tlu,orb_tlu,dr[it]))
				if flag_norm_orb:  orb_tlu[:]           = orb_tlu              / norm
				if flag_norm_C:    C[it][il].data[:,iu] = C[it][il].data[:,iu] / norm
				
def orth(orb,dr):
	""" |n'> = 1/Z ( |n> - \sum_{i=0}^{n-1} |i><i|n> ) """
	""" orb[it][il][iu,r] """
	for it,orb_t in orb.items():
		for il,orb_tl in enumerate(orb_t):	
			for iu1,orb_tlu1 in enumerate(orb_tl):
				for iu2 in range(iu1):
					orb_tlu1[:] -= orb_tl[iu2] * inner_product(orb_tlu1,orb_tl[iu2],dr[it])
				orb_tlu1[:] = orb_tlu1 / np.sqrt(inner_product(orb_tlu1,orb_tlu1,dr[it]))
				
def find_eigenvalue(Nl,Ne):
	""" E[il,ie] """
	E = np.zeros((Nl,Ne+Nl+1))
	for ie in range(1,Ne+Nl+1):
		E[0,ie] = ie*np.pi
	for il in range(1,Nl):
		jl = functools.partial(spherical_jn,il)
		for ie in range(1,Ne+Nl+1-il):
			E[il,ie] = fsolve( jl, (E[il-1,ie]+E[il-1,ie+1])/2 )
	return E[:,1:Ne+1]
	
def set_E(info_element):
	""" E[it][il,ie] """
	eigenvalue = { it:find_eigenvalue(info_element[it].Nl,info_element[it].Ne) for it in info_element }
	E = dict()
	for it in info_element:
		E[it] = torch.from_numpy(( eigenvalue[it]/info_element[it].Rcut ).astype("float64"))
	return E