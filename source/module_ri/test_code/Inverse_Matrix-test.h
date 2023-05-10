//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef INVERSE_MATRIX_TEST_H
#define INVERSE_MATRIX_TEST_H

#include "module_ri/Inverse_Matrix.h"
#include <LibRI/unittests/global/Tensor-test.h>

namespace Inverse_Matrix_Test
{
	template<typename Tdata>
	Tensor<Tdata> init_Tensor(const std::vector<size_t> &shape)
	{
		Tensor<Tdata> t(shape);
		for(size_t i=0; i<t.data->size(); ++i)
			t.ptr()[i] = i;
		return t;
	}
	
	template<typename Tdata>
	Tensor<Tdata> init_Tensor2(const std::vector<size_t> &shape)
	{
		Tensor<Tdata> t(shape);
		for(size_t i0=0; i0<t.shape[0]; ++i0)
			for(size_t i1=0; i1<t.shape[1]; ++i1)
				t(i0,i1) = i0==i1 ? i0+1 : (i0+i1)/10.0;
		return t;
	}

	template<typename Tdata>
	void test_input_output()
	{
		Inverse_Matrix<Tdata> inv;

		const size_t n_all = 5;
		const std::vector<size_t> n0 = {2,3};
		const std::vector<size_t> n1 = {1,2,2};
		
		Tensor<Tdata> m = init_Tensor<Tdata>({n_all,n_all});
		
		std::vector<std::vector<Tensor<Tdata>>> ms(n0.size(), std::vector<Tensor<Tdata>>(n1.size()));
		for(size_t Im0=0; Im0<n0.size(); ++Im0)
			for(size_t Im1=0; Im1<n1.size(); ++Im1)
				ms[Im0][Im1] = init_Tensor<Tdata>({n0[Im0], n1[Im1]});

		inv.input(m);
		std::cout<<inv.output()<<std::endl;
		
		inv.input(m);
		std::cout<<inv.output(n0, n1)<<std::endl;
		
		inv.input(ms);
		std::cout<<inv.output()<<std::endl;

		inv.input(ms);
		std::cout<<inv.output(n0, n1)<<std::endl;

		/*
			0       1       2       3       4
			5       6       7       8       9
			10      11      12      13      14
			15      16      17      18      19
			20      21      22      23      24
		*/
		/*
			0
			5

			1       2
			6       7

			3       4
			8       9

			10
			15
			20

			11      12
			16      17
			21      22

			13      14
			18      19
			23      24
		*/
		/*
			0       0       1       0       1
			1       2       3       2       3
			0       0       1       0       1
			1       2       3       2       3
			2       4       5       4       5
		*/
		/*
			0
			1

			0       1
			2       3

			0       1
			2       3

			0
			1
			2

			0       1
			2       3
			4       5

			0       1
			2       3
			4       5
		*/
	}

	template<typename Tdata>
	void test_inverse()
	{
		Tensor<Tdata> t = init_Tensor2<Tdata>({5,5});
		Inverse_Matrix<Tdata> inv;
		inv.input(t);
		inv.cal_inverse(Inverse_Matrix<Tdata>::Method::potrf);
		//inv.cal_inverse(Inverse_Matrix<Tdata>::Method::syev);
		Tensor<Tdata> tI = inv.output();

		std::cout<<t<<std::endl;
		std::cout<<tI<<std::endl;
		std::cout<<t*tI<<std::endl;
		std::cout<<tI*t<<std::endl;
	}
}

#endif