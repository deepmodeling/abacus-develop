#ifndef HS_MATRIX_H
#define HS_MATRIX_H

#include "../src_pw/tools.h"

// mohan add this file 2010-09-10
namespace HS_Matrix
{
	void saving_HS(const double *Hloc, const double* Sloc, bool bit, const int &out_hs);

	void save_HS(const double *H, const double *S, bool bit);

	void save_HS_complex(const complex<double> *H, const complex<double> *S, bool bit);

	void save_HSR_tr(const int current_spin); //LiuXh add 2019-07-15

	// jingan add 2021-6-4
	void save_HR_sparse(const int &current_spin);
	void save_SR_sparse(const int &current_spin);

// mohan comment out 2021-02-10
//	void save_HS_ccf(const int &iter, const int &Hnnz, const int *colptr_H, const int *rowind_H, 
//		const double *nzval_H, const double *nzval_S, bool bit);

	void saving_HS_complex(complex<double> *Hloc, complex<double>* Sloc, bool bit, const int &out_hs); //LiuXh, 2017-03-21

	void save_HS_complex(complex<double> *H, complex<double> *S, bool bit); //LiuXh, 2017-03-21
}

#endif
