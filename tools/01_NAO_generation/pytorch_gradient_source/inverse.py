import torch


def inverse_DC(A,B,C,D):
	tmp_A  = inverse(A)						# A^{-1}
	#print("tmp_A",tmp_A)
	tmp_AB = torch.mm(tmp_A,B)				# A^{-1} B
	#print("tmp_AB",tmp_AB)
	tmp_CA = torch.mm(C,tmp_A)				# C A^{-1}
	#print("tmp_CA",tmp_CA)
	tmp_X  = inverse(D-torch.mm(tmp_CA,B))	# ( D - C A^{-1} B )^{-1}
	#print("tmp_X",tmp_X)
	tmp_ABX = torch.mm(tmp_AB,tmp_X)		# A^{-1} B ( D - C A^{-1} B )^{-1}
	#print("tmp_ABX",tmp_ABX)
	tmp_XCA = torch.mm(tmp_X,tmp_CA)		# ( D - C A^{-1} B )^{-1} C A^{-1}
	#print("tmp_XCA",tmp_XCA)
	tmp_ABXCA = torch.mm(tmp_ABX,tmp_CA)	# A^{-1} B ( D - C A^{-1} B )^{-1} C A^{-1}
	#print("tmp_ABXCA",tmp_ABXCA)
	
	tmp_up   = torch.cat( [ tmp_A+tmp_ABXCA, -tmp_ABX ], dim=1 )
	tmp_down = torch.cat( [ -tmp_XCA, tmp_X ], dim=1 )
	I = torch.cat( [tmp_up,tmp_down], dim=0 )
	return I	

	
	
def inverse(M):
	
#	assert len(M.size()) == 2, "inverse must be 2D"
#	assert M.size()[0] == M.size()[1], "inverse row != column"
	
	L = M.size()[0]
	
	if L==1:
		return 1/M
		
	elif L==2:
		det  = torch.cat(list( M[:1,:1]*M[1:,1:] - M[:1,1:]*M[1:,:1] )*4).view(2,2)
		I_up   = torch.cat([M[1:,1:],-M[:1,1:]],dim=1)
		I_down = torch.cat([-M[1:,:1],M[:1,:1]],dim=1)
		I_all  = torch.cat([I_up,I_down],dim=0)
		return I_all/det
	
	elif L==3:
		threshold = 1e-10
		if M[0,0].abs() > threshold:
			return inverse_DC( M[:1,:1], M[:1,1:], M[1:,:1], M[1:,1:] )
		elif M[2,2].abs() > threshold:
			return inverse_DC( M[:2,:2], M[:2,2:], M[2:,:2], M[2:,2:] )
		else:
			raise ZeroDivisionError("matrix inverse")
			
	else:
		L2 = L//2
		return inverse_DC( M[:L2,:L2], M[:L2,L2:], M[L2:,:L2], M[L2:,L2:] )			