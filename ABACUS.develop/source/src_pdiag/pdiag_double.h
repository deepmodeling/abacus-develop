#ifndef PDIAG_DOUBLE_H
#define PDIAG_DOUBLE_H

#include "src_pw/tools.h"
#include "pdiag_basic.h"

class Pdiag_Double : public Pdiag_Basic
{
	public:
	Pdiag_Double();
	~Pdiag_Double();

	LocalMatrix MatrixInfo;

	int nrow;
	int ncol;
	long nloc;

#ifdef __MPI
    int blacs_ctxt;
	int mpi_comm_rows, mpi_comm_cols;
    int desc[9];
#endif

	// this subroutine needs reconstruction in near future -- mohan note 2021-03
	void diago_double_begin(const int &ik, double **wfc, matrix &wfc_2d,
		double *h_mat, double *s_mat, double *ekb);			// Peize Lin add wfc_2d 2019-01-17
	
	// this subroutine needs reconstruction in near future -- mohan note 2021-03
	void diago_complex_begin(const int &ik, complex<double> **wfc, ComplexMatrix &wfc_2d,
		complex<double> *ch_mat, complex<double> *cs_mat, double *ekb);			// Peize Lin add wfc_2d 2019-01-17

#ifdef __MPI
	void readin(const string &fa, const string &fb, const int &nlocal, double *eigen, double *eigvr);
#endif

	// be called
	void divide_HS_2d
	(
#ifdef __MPI
		MPI_Comm DIAG_WORLD
#endif
	);

#ifdef __MPI
	MPI_Comm comm_2D;
#endif

	int nb;
	int dim0;
	int dim1;

};

#endif
