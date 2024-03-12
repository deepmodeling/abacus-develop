#include "ylm_real.h"
#include "numerical_basis.h"

Bessel_Basis Numerical_Basis::bessel_basis;

Numerical_Basis::Numerical_Basis(){}

Numerical_Basis::~Numerical_Basis(){}

//============================================================
// MEMBER FUNCTION : 
// NAME : init
// Firstly, use check(0) to call bessel_basis,init
// to generate TableOne
// generate overlap, use jlq3d.
//============================================================

// The function is called in run_fp.cpp.
void Numerical_Basis::init_table(void)
{
	TITLE("Numerical_Basis","output_overlap");

	Numerical_Basis::bessel_basis.init();
	
	return;
}

