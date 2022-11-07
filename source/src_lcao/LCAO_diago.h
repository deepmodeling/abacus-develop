#ifndef DIAGO_LCAO_MATRIX_H
#define DIAGO_LCAO_MATRIX_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include "src_lcao/LCAO_matrix.h"
#include "src_lcao/local_orbital_wfc.h"
#include "src_pdiag/pdiag_double.h"

class Diago_LCAO_Matrix : Pdiag_Double
{
	public:
	Diago_LCAO_Matrix(LCAO_Matrix* lm);
	~Diago_LCAO_Matrix();

private:
    
    LCAO_Matrix* LM;

	void using_LAPACK(const int &ik, Local_Orbital_wfc &lowf, double* ekb_ik)const;
	void using_LAPACK_complex(const int &ik, std::complex<double> **wfc_k_grid ,ModuleBase::ComplexMatrix &wfc_k)const;
	void using_CG(const int &ik, double **c)const;


};

#endif
