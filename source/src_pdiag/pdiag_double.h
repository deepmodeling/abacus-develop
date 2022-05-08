#ifndef PDIAG_DOUBLE_H
#define PDIAG_DOUBLE_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include "diag_scalapack_gvx.h"
#include "module_orbital/parallel_orbitals.h"
#include "src_lcao/local_orbital_wfc.h"

class Pdiag_Double
{

public:
    static int out_wfc_lcao;
    static int out_mat_hs; // mohan add 2010-09-02
    static int out_mat_hsR; // LiuXh add 2019-07-16
    Pdiag_Double();
    ~Pdiag_Double();
    
protected:
	// this subroutine needs reconstruction in near future -- mohan note 2021-03
	//void diago_double_begin(const int &ik, double **wfc, ModuleBase::matrix &wfc_2d,
	//	double *h_mat, double *s_mat, double *ekb);			// Peize Lin add wfc_2d 2019-01-17
    void diago_double_begin(const int& ik, Local_Orbital_wfc &lowf,
        double* h_mat, double* s_mat,
        double* Stmp, double* ekb); //LiuXh add 2021-09-06, clear memory, totwfc not used now

	// this subroutine needs reconstruction in near future -- mohan note 2021-03
    void diago_complex_begin(const int& ik, Local_Orbital_wfc &lowf,
        std::complex<double>* ch_mat, std::complex<double>* cs_mat,
        std::complex<double>* Stmp, double* ekb);			// Peize Lin add wfc_2d 2019-01-17

    Diag_Scalapack_gvx diag_scalapack_gvx;			// Peize Lin add 2021.11.02

    /// output control parameters in diago (static)
    // mohan add 2010-09-10
	// output local wave functions.
	// put it here because if we 
	// use HPSEPS, the wave functions
	// is needed to be collected first.

private:

    const Parallel_Orbitals* ParaV;
    
#ifdef __MPI
	//void gath_eig(MPI_Comm comm,int n,double **c,double *Z);
	void gath_eig(MPI_Comm comm,int n,double *Z); //LiuXh add 2021-09-06, clear memory, totwfc not used now
	void gath_eig_complex(MPI_Comm comm,int n,std::complex<double> **c,std::complex<double> *Z, const int &ik); //mohan add 2012-01-09
	void gath_full_eig(MPI_Comm comm,int n,double **c,double *Z);
	void gath_full_eig_complex(MPI_Comm comm,int n,std::complex<double> **c, std::complex<double> *Z);
#endif

};

#endif
