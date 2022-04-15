#ifndef EVOLVE_LCAO_MATRIX_H
#define EVOLVE_LCAO_MATRIX_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/complexmatrix.h"
#include "src_lcao/local_orbital_wfc.h"
#include "src_lcao/LCAO_matrix.h"

class Evolve_LCAO_Matrix
{
	public:
	Evolve_LCAO_Matrix(LCAO_Matrix *lm);
	~Evolve_LCAO_Matrix();
	
    void evolve_complex_matrix(const int& ik, Local_Orbital_wfc &lowf)const;

private:
    LCAO_Matrix* LM;

	void using_LAPACK_complex(const int &ik, std::complex<double>*** wfc_k_grid, std::complex<double> **c_init)const;
#ifdef __MPI
	int using_ScaLAPACK_complex(const int &ik, ModuleBase::ComplexMatrix &wfc_2d)const;
#endif
};

#endif
