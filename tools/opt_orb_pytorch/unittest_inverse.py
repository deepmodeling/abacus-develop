import unittest
import inverse
import torch

class unittest_inverse(unittest.TestCase):

	def inverse_test(self,a,ai_true):
	
		a=torch.Tensor(a)
		a=torch.autograd.Variable(a)
		
		ai_test=inverse.inverse(a)
		
		ai_true = torch.Tensor(ai_true)
		ai_true=torch.autograd.Variable(ai_true)
		
		self.assertFalse((ai_test!=ai_true).data.sum())
		
	def test_inverse_1(self):
		self.inverse_test( 
			[[1,2],[2,3]], 
			[[-3,2],[2,-1]] )
	def test_inverse_2(self):
		self.inverse_test( 
			[[1,2,3],[2,4,5],[3,5,6]], 
			[[1,-3,2],[-3,3,-1],[2,-1,0]] )

if __name__ == '__main__':
    unittest.main()