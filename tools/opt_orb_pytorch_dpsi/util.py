def ND_list(*sizes,element=None):
	size_1,*size_other = sizes
	l = [element] * size_1
	if size_other:
		for i in range(len(l)):
			l[i] = ND_list(*size_other,element=element)
	else:
		if element in ["dict()","list()"]:	
			for i in range(size_1):	
				l[i] = eval(element)
	return l
	

def ignore_line(file,N):
	for _ in range(N):	
		file.readline()
		
		
class Info:
	def Nm(self,il): return 2*il+1
	def __str__(self):
		return "\n".join([name+"\t"+str(value) for name,value in self.__dict__.items()])
	__repr__=__str__
	
def change_to_cuda(s):
	if isinstance(s,list):
		return [change_to_cuda(x) for x in s]
	elif isinstance(s,dict):
		return {i:change_to_cuda(x) for i,x in s.items()}
	elif isinstance(s,torch.Tensor):
		return s.cuda()
	elif isinstance(s,torch_complex.ComplexTensor):
		return torch_complex.ComplexTensor( change_to_cuda(s.real), change_to_cuda(s.imag) ) 
	else:
		print(s)
		raise TypeError("change_to_cuda")

def update0(t):
	return t.masked_fill(mask=(t==0), value=1E-10)

def Nm(il):
	return 2*il+1