#ifndef HAMILT_PW_H
#define HAMILT_PW_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/complexmatrix.h"

class Hamilt_PW
{

public:

    Hamilt_PW();
    ~Hamilt_PW();

    static int moved;

	void allocate(
		const int &npwx, // number of plane wave (max)
		const int &npol, // polarization
		const int &nkb,  // number of non-local pseudopotential projectors
		const int &nrxx); // number of grids on this processor

    void cal_err
    (
        const int &npw,
        ModuleBase::ComplexMatrix &psi,
        const int &nband,
        double *em,
        double *err
    );

    void init_k(const int ik);

	friend class Diago_David;
	friend class Diago_CG;
	friend class Exx_Lip;
	friend class Hamilt;
    friend class Stochastic_Iter;

    void diagH_subspace(const int ik,
                  const int nstart,
                  const int nbnd,
                  const ModuleBase::ComplexMatrix &psi,
                  ModuleBase::ComplexMatrix &evc,
                  double *en);

    void h_1psi(
        const int npw,
        const std::complex<double> *psi1d,
        std::complex<double> *hpsi,
        std::complex<double> *spsi);

    void h_psi(
		const std::complex<double> *psi,
		std::complex<double> *hpsi,
		const int m = 1); // qianrui add a default parameter 2021-3-31

    void s_1psi(
        const int npw,
        const std::complex < double> *psi,
        std::complex < double> *spsi);

    int *GR_index;

    std::complex<double> *psi_one;

    // hpsi and spsi
    std::complex<double> *hpsi;
    std::complex<double> *spsi;

	// add contributions of h*psi from
	// non-local pseduopotentials
    void add_nonlocal_pp(
		std::complex<double> *hpsi,
		const std::complex<double> *becp,
		const int m);

    void diag_zheev
    (
        const int& npw,
        ModuleBase::ComplexMatrix& psi,
        const int& nband,
        double *em,
        double *err
    ) ;

};

#endif
