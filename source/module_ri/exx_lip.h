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
#include "module_hamilt_pw/hamilt_pwdft/wfinit.h"
class K_Vectors;
class UnitCell;
template<typename T, typename Device = base_device::DEVICE_CPU>
class Exx_Lip
{
    using Real = typename GetTypeReal<T>::type;
public:
    Exx_Lip(const Exx_Info::Exx_Info_Lip& info_in);
    ~Exx_Lip();

    const Exx_Info::Exx_Info_Lip& info;

    Exx_Lip(const Exx_Info::Exx_Info_Lip& info_in,
        const ModuleSymmetry::Symmetry& symm,
        K_Vectors* kv_ptr_in,
        //   wavefunc* wf_ptr_in,
        psi::WFInit<T, Device>* wf_ptr_in,
        psi::Psi<T, Device>* kspw_psi_ptr_in,
        const ModulePW::PW_Basis_K* wfc_basis_in,
        const ModulePW::PW_Basis* rho_basis_in,
        const Structure_Factor& sf,
        const UnitCell* ucell_ptr_in,
        const elecstate::ElecState* pelec_in);
    
    // void cal_exx(const int& nks);
    void cal_exx();
    const std::vector<std::vector<std::vector<T>>>& get_exx_matrix() const
    {
        return exx_matrix;
    }
    Real get_exx_energy() const
    {
        return exx_energy;
    }

    void write_q_pack() const;

    void set_hvec(const int ik, const T* const hvec, const int naos, const int nbands)
    {
        memcpy(&(*this->k_pack->hvec_array)(ik, 0, 0), hvec, sizeof(T) * naos * nbands);
    }
    psi::Psi<T,Device> get_hvec() const
    {
        return *k_pack->hvec_array;
    }

private:

    int gzero_rank_in_pool;

    // template<typename T, typename Device = base_device::DEVICE_CPU>
    struct k_package
    {
        K_Vectors* kv_ptr = nullptr;
        // wavefunc* wf_ptr;
        psi::Psi<T, Device>* kspw_psi_ptr = nullptr;  ///< PW  wavefunction
        psi::WFInit<T, Device>* wf_ptr = nullptr;
        ModuleBase::matrix wf_wg;

        /// @brief LCAO wavefunction, the eigenvectors from lapack diagonalization
        psi::Psi<T, Device>* hvec_array = nullptr;
        const elecstate::ElecState* pelec = nullptr;
    } *k_pack = nullptr, * q_pack = nullptr;

    int iq_vecik;

    std::vector<std::vector<T>> phi;
    std::vector<std::vector<std::vector<T>>> psi;
    std::vector<Real> recip_qkg2;
    Real sum2_factor;
    std::vector<T> b;
    std::vector<T> b0;
    std::vector<T> sum1;
    std::vector<std::vector<T>> sum3;

    std::vector<std::vector<std::vector<T>>> exx_matrix;
    Real exx_energy = 0.0;

	void wf_wg_cal();
    void phi_cal(k_package* kq_pack, int ikq);
    void psi_cal();
    void judge_singularity( int ik);
	void qkg2_exp(int ik, int iq);
	void b_cal(int ik, int iq, int ib);
	void sum3_cal(int iq, int ib);
	void b_sum(int iq, int ib);
	void sum_all(int ik);
	void exx_energy_cal();
    // void read_q_pack(const ModuleSymmetry::Symmetry& symm,
    //                  const ModulePW::PW_Basis_K* wfc_basis,
    //                  const Structure_Factor& sf);

    //2*pi*i
    T two_pi_i = Real(ModuleBase::TWO_PI) * T(0.0, 1.0);
public:
    const ModulePW::PW_Basis* rho_basis = nullptr;
    const ModulePW::PW_Basis_K* wfc_basis = nullptr;

    const UnitCell* ucell_ptr = nullptr;
};


#endif