import torch

class ComplexTensor:
	def __init__(self,real,imag):
		self.real = real
		self.imag = imag
	
	def view(self,*args,**kwargs):
		return ComplexTensor( self.real.view(*args,**kwargs), self.imag.view(*args,**kwargs) )
	def t(self,*args,**kwargs):
		return ComplexTensor( self.real.t(*args,**kwargs), self.imag.t(*args,**kwargs) )
#	def transpose(self,*args,**kwargs):
#		return ComplexTensor( self.real.transpose(*args,**kwargs), self.imag.transpose(*args,**kwargs) )
	def __getitem__(self,*args,**kwargs):
		return ComplexTensor( self.real.__getitem__(*args,**kwargs), self.imag.__getitem__(*args,**kwargs) )
	def __str__(self):
		return "<{0};{1}>".format(self.real, self.imag)
	__repr__=__str__
#	def size(self,*args,**kwargs):
#		return ComplexTensor( self.real.size(*args,**kwargs), self.imag.size(*args,**kwargs) )

	def conj(self):
		return ComplexTensor( self.real, -self.imag )
		
	def mm( self,x2, *args,**kwargs ):
		return mm( self,x2, *args,**kwargs )
		
		
def dot( x1,x2, *args,**kwargs ):
	if isinstance(x1,ComplexTensor):
		if isinstance(x2,ComplexTensor):
			return ComplexTensor( torch.dot( x1.real,x2.real, *args,**kwargs ) - torch.dot( x1.imag,x2.imag, *args,**kwargs ), torch.dot( x1.real,x2.imag, *args,**kwargs ) + torch.dot( x1.imag,x2.real, *args,**kwargs ) )
		else:
			return ComplexTensor( torch.dot( x1.real,x2, *args,**kwargs ), torch.dot( x1.imag,x2, *args,**kwargs ) )
	else:
		if isinstance(x2,ComplexTensor):
			return ComplexTensor( torch.dot( x1,x2.real, *args,**kwargs ), torch.dot( x1,x2.imag, *args,**kwargs ) )
		else:
			return torch.dot( x1,x2, *args,**kwargs )
def mv( x1,x2, *args,**kwargs ):
	if isinstance(x1,ComplexTensor):
		if isinstance(x2,ComplexTensor):
			return ComplexTensor( torch.mv( x1.real,x2.real, *args,**kwargs ) - torch.mv( x1.imag,x2.imag, *args,**kwargs ), torch.mv( x1.real,x2.imag, *args,**kwargs ) + torch.mv( x1.imag,x2.real, *args,**kwargs ) )
		else:
			return ComplexTensor( torch.mv( x1.real,x2, *args,**kwargs ), torch.mv( x1.imag,x2, *args,**kwargs ) )
	else:
		if isinstance(x2,ComplexTensor):
			return ComplexTensor( torch.mv( x1,x2.real, *args,**kwargs ), torch.mv( x1,x2.imag, *args,**kwargs ) )
		else:
			return torch.mv( x1,x2, *args,**kwargs )
def mm( x1,x2, *args,**kwargs ):
	if isinstance(x1,ComplexTensor):
		if isinstance(x2,ComplexTensor):
			return ComplexTensor( torch.mm( x1.real,x2.real, *args,**kwargs ) - torch.mm( x1.imag,x2.imag, *args,**kwargs ), torch.mm( x1.real,x2.imag, *args,**kwargs ) + torch.mm( x1.imag,x2.real, *args,**kwargs ) )
		else:
			return ComplexTensor( torch.mm( x1.real,x2, *args,**kwargs ), torch.mm( x1.imag,x2, *args,**kwargs ) )
	else:
		if isinstance(x2,ComplexTensor):
			return ComplexTensor( torch.mm( x1,x2.real, *args,**kwargs ), torch.mm( x1,x2.imag, *args,**kwargs ) )
		else:
			return torch.mm( x1,x2, *args,**kwargs )

			
def cat( xs, *args,**kwargs ):
	if isinstance(xs[0],ComplexTensor):
		xs_real = [];	xs_imag = []
		for x in xs:
			xs_real.append(x.real)
			xs_imag.append(x.imag)
		return ComplexTensor( torch.cat(xs_real,*args,**kwargs), torch.cat(xs_imag,*args,**kwargs) ) 
	else:
		return torch.cat(xs,*args,**kwargs)
	
	
	
def inverse(M):
	if isinstance(M,ComplexTensor):
		A=M.real
		B=M.imag
		tmp_AB = torch.mm(A.inverse(),B)						# A^{-1} B
		tmp_X  = (A+torch.mm(B,tmp_AB)).inverse()				# ( A + B A^{-1} B )^{-1}
		return ComplexTensor( tmp_X, -torch.mm(tmp_AB,tmp_X) )
	else:
		return M.inverse()