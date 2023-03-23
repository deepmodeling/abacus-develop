#ifndef WF_IGK_H
#define WF_IGK_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/intarray.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_basis/module_pw/pw_basis_k.h"

class WF_igk 
{
	public:

    WF_igk();
    ~WF_igk();

	//---------------------------------------------------
	// npwx: max npw
	// npw
	// igk: [nks, npw_max]
	//---------------------------------------------------
    int npwx;
    int npw;
    // ModuleBase::IntArray igk;
#ifdef __CUDA
    double *d_g2kin;
#endif

	public:

    // Calculate kinetic energy
    double* get_qvec_cartesian(const int &ik);

    ModuleBase::Vector3<double> get_1qvec_cartesian(const int ik,const int ig)const;

    std::complex<double>* get_sk(const int ik, const int it, const int ia, ModulePW::PW_Basis_K* wfc_basis)const;

    template<typename FPTYPE, typename Device>
    void get_sk(Device * ctx, const int ik, ModulePW::PW_Basis_K* wfc_basis, std::complex<FPTYPE> * sk)const;

    // pengfei 2016-11-23
    std::complex<double>* get_skq(int ik, int it, int ia, ModuleBase::Vector3<double> q);

};
#endif
