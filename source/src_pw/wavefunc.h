#ifndef WAVEFUNC_H
#define WAVEFUNC_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include "wf_atomic.h"

class wavefunc : public WF_atomic
{
	public:

    wavefunc();
    ~wavefunc();

    // allocate memory
    void allocate(const int nks);
    void allocate_ekb_wg(const int nks);

    int out_wfc_pw; //qianrui modify 2020-10-19
    int out_wfc_r=0; // Peize Lin add 2021.11.21

    // et    : (nks,nbnd),eigenvalues of the hamiltonian
    // wg	 : the weight of each k point and band
	double **ekb;		   // band energy at each k point, each band.
	bool   allocate_ekb;   // flag
    ModuleBase::matrix wg;

    // init_wfc : "random",or "atomic" or "file"
    std::string init_wfc;
	int mem_saver; //1: save evc when doing nscf calculation.
    void wfcinit(void);// from wfcinit.f90
    void wfcinit_k();
    int get_starting_nw(void)const;

	// wanf2: save PAO orbitals,
	//void PAO_in_pw_k(const int &ik, ModuleBase::ComplexMatrix &wvf);

	// wanf2: save given localized orbitals. 
	// mohan_to_qianrui: move the LCAO_in_pw_k and LCAO_in_pw_k_q to wavefunc_in_pw.h 
	void LCAO_in_pw_k(const int &ik, ModuleBase::ComplexMatrix &wvf);
	void LCAO_in_pw_k_q(const int &ik, ModuleBase::ComplexMatrix &wvf, ModuleBase::Vector3<double> q);   // pengfei 2016-11-23

	// evc: get the initial wave functions from diagnalized the PAO
	// orbitals first.
	void diago_PAO_in_pw_k(const int &ik, ModuleBase::ComplexMatrix &wvf);

	// used if k dependent staff is ready.
	void prepare_k(void);

	void diago_PAO_in_pw_k2(const int &ik, ModuleBase::ComplexMatrix &wvf);

    int get_R(int ix, int iy, int iz);     // pengfei 2016-11-23

    int iw2it(int iw);
    int iw2ia(int iw);

    void init_after_vc(const int nks); //LiuXh 20180515

    private: // pengfei 2016-11-23

    ModuleBase::Vector3<int> ***R;
    int ** Rmax;
};

#endif //wavefunc
