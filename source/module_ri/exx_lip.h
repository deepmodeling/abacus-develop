//==========================================================
// AUTHOR : Peize Lin
// DATE : 2015-03-10
//==========================================================
#ifndef EXX_LIP_H
#define EXX_LIP_H

#include "module_base/complexmatrix.h"
#include "module_base/vector3.h"
#include "module_hamilt_general/module_xc/exx_info.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_elecstate/elecstate.h"
#include "module_cell/module_symmetry/symmetry.h"

class K_Vectors;
class wavefunc;
class UnitCell;

class Exx_Lip
{
public:
	Exx_Lip( const Exx_Info::Exx_Info_Lip &info_in );
	~Exx_Lip();

	const Exx_Info::Exx_Info_Lip &info;

    void init(const ModuleSymmetry::Symmetry& symm,
              K_Vectors* kv_ptr_in,
              wavefunc* wf_ptr_in,
              const ModulePW::PW_Basis_K* wfc_basis_in,
              const ModulePW::PW_Basis* rho_basis_in,
              const Structure_Factor& sf,
              const UnitCell* ucell_ptr_in,
              const elecstate::ElecState* pelec_in);
    // void cal_exx(const int& nks);
    const std::complex<double>* const* const* get_exx_matrix() const
    {
        return exx_matrix;
    }
    double get_exx_energy() const
    {
        return exx_energy;
    }

    void write_q_pack() const;

  private:
    bool init_finish;

    int gzero_rank_in_pool;

    struct k_package
    {
        K_Vectors* kv_ptr;
        wavefunc* wf_ptr;
        ModuleBase::matrix wf_wg;
        ModuleBase::ComplexMatrix* hvec_array;
        const elecstate::ElecState* pelec;
    } *k_pack, *q_pack;

    int iq_vecik;

    std::complex<double>** phi;
    std::complex<double>*** psi;
    double *recip_qkg2;
	double sum2_factor;
	std::complex<double> *b;
	std::complex<double> *b0;
	std::complex<double> *sum1;
	std::complex<double> **sum3;

	std::complex<double> ***exx_matrix;
	double exx_energy;

	void wf_wg_cal();
	void phi_cal(k_package *kq_pack, int ikq);
    // void psi_cal();
    void judge_singularity( int ik);
	void qkg2_exp(int ik, int iq);
	void b_cal(int ik, int iq, int ib);
	void sum3_cal(int iq, int ib);
	void b_sum(int iq, int ib);
	void sum_all(int ik);
	void exx_energy_cal();
    void read_q_pack(const ModuleSymmetry::Symmetry& symm,
                     const ModulePW::PW_Basis_K* wfc_basis,
                     const Structure_Factor& sf);

  public:
    const ModulePW::PW_Basis* rho_basis;
    const ModulePW::PW_Basis_K* wfc_basis;

    const UnitCell* ucell_ptr;
};


#endif
