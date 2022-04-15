//=====================
// AUTHOR : Peize Lin
// DATE : 2021-11-02
//=====================

#ifndef DIAG_SCALAPACK_GVX
#define DIAG_SCALAPACK_GVX

#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include <complex>
#include <utility>
#include <vector>

class Diag_Scalapack_gvx
{

public:
	void pdsygvx_diag(
		const int*const desc,
		const int ncol,
		const int nrow,
		const double*const h_mat,
		const double*const s_mat,
		double*const ekb,
		ModuleBase::matrix &wfc_2d);
	void pzhegvx_diag(
		const int*const desc,
		const int ncol,
		const int nrow,
		const std::complex<double>*const h_mat,
		const std::complex<double>*const s_mat,
		double*const ekb,
		ModuleBase::ComplexMatrix &wfc_2d);

	std::pair<int,std::vector<int>> pdsygvx_once(
		const int*const desc,
		const int ncol,
		const int nrow,
		const double*const h_mat,
		const double*const s_mat,
		double*const ekb,
		ModuleBase::matrix &wfc_2d) const;
	std::pair<int,std::vector<int>> pzhegvx_once(
		const int*const desc,
		const int ncol,
		const int nrow,
		const std::complex<double>*const h_mat,
		const std::complex<double>*const s_mat,
		double*const ekb,
		ModuleBase::ComplexMatrix &wfc_2d) const;

private:
	int degeneracy_max = 12;			// For reorthogonalized memory. 12 followes siesta.

	void post_processing(const int info, const std::vector<int> &vec);
};

#endif