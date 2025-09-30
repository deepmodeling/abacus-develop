#ifndef __WRITE_VXC_R_H_
#define __WRITE_VXC_R_H_
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_operator_lcao/op_dftu_lcao.h"
#include "source_lcao/module_operator_lcao/veff_lcao.h"
#include "source_lcao/spar_hsr.h"
#include "source_io/write_HS_sparse.h"
#ifdef __EXX
#include "source_lcao/module_operator_lcao/op_exx_lcao.h"
#include "source_lcao/module_ri/RI_2D_Comm.h"
#endif

#ifndef TGINT_H
#define TGINT_H
template <typename T>
struct TGint;

template <>
struct TGint<double>
{
    using type = Gint_Gamma;
};

template <>
struct TGint<std::complex<double>>
{
    using type = Gint_k;
};
#endif

namespace ModuleIO
{

#ifndef SET_GINT_POINTER_H
#define SET_GINT_POINTER_H
template <typename T>
void set_gint_pointer(Gint_Gamma& gint_gamma, Gint_k& gint_k, typename TGint<T>::type*& gint);

template <>
void set_gint_pointer<double>(Gint_Gamma& gint_gamma, Gint_k& gint_k, typename TGint<double>::type*& gint)
{
    gint = &gint_gamma;
}

template <>
void set_gint_pointer<std::complex<double>>(Gint_Gamma& gint_gamma,
                                            Gint_k& gint_k,
                                            typename TGint<std::complex<double>>::type*& gint)
{
    gint = &gint_k;
}
#endif

template <typename TR> std::set<Abfs::Vector3_Order<int>> get_R_range(const hamilt::HContainer<TR>& hR)
{
    std::set<Abfs::Vector3_Order<int>> all_R_coor;

    return all_R_coor;
}

template <typename TR>
std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, TR>>>
cal_HR_sparse(const hamilt::HContainer<TR>& hR,
    const int current_spin,
    const double sparse_thr);

template <>
std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>
cal_HR_sparse(const hamilt::HContainer<double>& hR,
    const int current_spin,
    const double sparse_thr)
{
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> target;
    sparse_format::cal_HContainer_d(*hR.get_paraV(), current_spin, sparse_thr, hR, target);
    return target;
}
template <>
std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>
cal_HR_sparse(const hamilt::HContainer<std::complex<double>>& hR,
    const int current_spin,
    const double sparse_thr)
{
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> target;
    sparse_format::cal_HContainer_cd(*hR.get_paraV(), current_spin, sparse_thr, hR, target);
    return target;
}

/// @brief  write the Vxc matrix in KS orbital representation, usefull for GW calculation
/// including terms: local/semi-local XC, EXX, DFTU
template <typename TK, typename TR>
void write_Vxc_R(const int nspin,
    const Parallel_Orbitals* pv,
    const UnitCell& ucell,
    Structure_Factor& sf,
    surchem& solvent,
    const ModulePW::PW_Basis& rho_basis,
    const ModulePW::PW_Basis& rhod_basis,
    const ModuleBase::matrix& vloc,
    const Charge& chg,
    Gint_Gamma& gint_gamma,
    Gint_k& gint_k,
    const K_Vectors& kv,
    const std::vector<double>& orb_cutoff,
    Grid_Driver& gd,
#ifdef __EXX
    const std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* const Hexxd,
    const std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* const Hexxc,
#endif
const double sparse_thr=1e-10)
{
    ModuleBase::TITLE("ModuleIO", "write_Vxc_R");
    // 1. real-space xc potential
    // ModuleBase::matrix vr_xc(nspin, chg.nrxx);
    double etxc = 0.0;
    double vtxc = 0.0;
    // elecstate::PotXC* potxc(&rho_basis, &etxc, vtxc, nullptr);
    // potxc.cal_v_eff(&chg, &ucell, vr_xc);
    elecstate::Potential* potxc
        = new elecstate::Potential(&rhod_basis, &rho_basis, &ucell, &vloc, &sf, &solvent, &etxc, &vtxc);
    std::vector<std::string> compnents_list = {"xc"};
    potxc->pot_register(compnents_list);
    potxc->update_from_charge(&chg, &ucell);

    // 2. allocate H(R)
    // (the number of hR: 1 for nspin=1, 4; 2 for nspin=2)
    int nspin0 = (nspin == 2) ? 2 : 1;
    std::vector<hamilt::HContainer<TR>> vxcs_R_ao(nspin0, hamilt::HContainer<TR>(ucell, pv));   // call move constructor
#ifdef __EXX
    std::array<int, 3> Rs_period = { kv.nmp[0], kv.nmp[1], kv.nmp[2] };
    const auto cell_nearest = hamilt::init_cell_nearest(ucell, Rs_period);
#endif
    for (int is = 0; is < nspin0; ++is)
    {
        if (std::is_same<TK, double>::value) { vxcs_R_ao[is].fix_gamma(); }
#ifdef __EXX
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            GlobalC::exx_info.info_ri.real_number ?
                hamilt::reallocate_hcontainer(*Hexxd, &vxcs_R_ao[is], &cell_nearest) :
                hamilt::reallocate_hcontainer(*Hexxc, &vxcs_R_ao[is], &cell_nearest);
        }
#endif
    }

    // 3. calculate the Vxc(R)
    hamilt::HS_Matrix_K<TK> vxc_k_ao(pv, 1); // only hk is needed, sk is skipped
    typename TGint<TK>::type* gint = nullptr;
    set_gint_pointer<TK>(gint_gamma, gint_k, gint);
    std::vector<hamilt::Veff<hamilt::OperatorLCAO<TK, TR>>*> vxcs_op_ao(nspin0);
    for (int is = 0; is < nspin0; ++is)
    {
        vxcs_op_ao[is] = new hamilt::Veff<hamilt::OperatorLCAO<TK, TR>>(gint,
            &vxc_k_ao, kv.kvec_d, potxc, &vxcs_R_ao[is], &ucell, orb_cutoff, &gd, nspin);
        vxcs_op_ao[is]->contributeHR();
#ifdef __EXX
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            GlobalC::exx_info.info_ri.real_number ?
                RI_2D_Comm::add_HexxR(is, GlobalC::exx_info.info_global.hybrid_alpha, *Hexxd, *pv, ucell.get_npol(), vxcs_R_ao[is], &cell_nearest) :
                RI_2D_Comm::add_HexxR(is, GlobalC::exx_info.info_global.hybrid_alpha, *Hexxc, *pv, ucell.get_npol(), vxcs_R_ao[is], &cell_nearest);
        }
#endif
    }

    // test: fold Vxc(R) and check whether it is equal to Vxc(k)
    // for (int ik = 0; ik < kv.get_nks(); ++ik)
    // {
    //     vxc_k_ao.set_zero_hk();
    //     dynamic_cast<hamilt::OperatorLCAO<TK, TR>*>(vxcs_op_ao[kv.isk[ik]])->contributeHk(ik);

    //     // output Vxc(k) (test)
    //     const TK* const hk = vxc_k_ao.get_hk();
    //     std::cout << "ik=" << ik << ", Vxc(K): " << std::endl;
    //     for (int i = 0; i < pv->get_row_size(); i++)
    //     {
    //         for (int j = 0; j < pv->get_col_size(); j++)
    //         {
    //             std::cout << hk[j * pv->get_row_size() + i] << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    // }

    // 4. write Vxc(R) in csr format
    for (int is = 0; is < nspin0; ++is)
    {
        std::set<Abfs::Vector3_Order<int>> all_R_coor = sparse_format::get_R_range(vxcs_R_ao[is]);
        const std::string filename = "Vxc_R_spin" + std::to_string(is);
        ModuleIO::save_sparse(
            cal_HR_sparse(vxcs_R_ao[is], is, sparse_thr),
            all_R_coor,
            sparse_thr,
            false, //binary
            PARAM.globalv.global_out_dir + filename + ".csr",
            *pv,
            filename,
            -1,
            true);  //all-reduce
    }
}
} // namespace ModuleIO
#endif
