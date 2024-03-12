//==========================================================
// AUTHOR : mohan
// DATE : 2009-4-2
// Last Modify: 2009-08-28
//==========================================================
#ifndef NUMERICAL_BASIS_H
#define NUMERICAL_BASIS_H
#include "../src_spillage/common.h"
#include "bessel_basis.h"
//==========================================================
// CLASS :
// NAME :  Numerical_Basis 
//==========================================================
class Numerical_Basis
{
	public:
	Numerical_Basis();
	~Numerical_Basis();

	void init_table(void);

	static Bessel_Basis bessel_basis;

	private:


};

#endif
