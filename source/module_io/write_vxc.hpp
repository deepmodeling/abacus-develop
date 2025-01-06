#ifndef __WRITE_VXC_H_
#define __WRITE_VXC_H_
#include "module_parameter/parameter.h"
#include "module_base/parallel_reduce.h"
#include "module_base/module_container/base/third_party/blas.h"
#include "module_base/scalapack_connector.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/op_dftu_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/veff_lcao.h"
#include "module_psi/psi.h"
#include "module_io/write_HS.h"

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

inline void gint_vl(Gint_Gamma& gg, Gint_inout& io)
{
    gg.cal_vlocal(&io, false);
};

inline void gint_vl(Gint_k& gk, Gint_inout& io, ModuleBase::matrix& wg)
{
    gk.cal_gint(&io);
};

inline void set_para2d_MO(const Parallel_Orbitals& pv, const int nbands, Parallel_2D& p2d)
{
    std::ofstream ofs;
#ifdef __MPI
    p2d.set(nbands, nbands, pv.nb, pv.blacs_ctxt);
#else
    p2d.set_serial(nbands, nbands);
#endif
}

template <typename T>
inline std::vector<T> cVc(T* V,
    T* c,
    const int nbasis,
    const int nbands,
    const Parallel_Orbitals& pv,
    const Parallel_2D& p2d)
{
    std::vector<T> Vc(pv.nloc_wfc, 0.0);
    char transa = 'N';
    char transb = 'N';
    const T alpha = (T)1.0;
    const T beta = (T)0.0;
#ifdef __MPI
    const int i1 = 1;
    ScalapackConnector::gemm(transa, transb,
        nbasis, nbands, nbasis,
        alpha, V, i1, i1, pv.desc,
        c, i1, i1, pv.desc_wfc,
        beta, Vc.data(), i1, i1, pv.desc_wfc);
#else
    container::BlasConnector::gemm(transa, transb, nbasis, nbands, nbasis, alpha, V, nbasis, c, nbasis, beta, Vc.data(), nbasis);
#endif
    std::vector<T> cVc(p2d.nloc, 0.0);
    transa = (std::is_same<T, double>::value ? 'T' : 'C');
#ifdef __MPI
    ScalapackConnector::gemm(transa, transb,
        nbands, nbands, nbasis,
        alpha, c, i1, i1, pv.desc_wfc,
        Vc.data(), i1, i1, pv.desc_wfc,
        beta, cVc.data(), i1, i1, p2d.desc);
#else
    container::BlasConnector::gemm(transa, transb, nbands, nbands, nbasis, alpha, c, nbasis, Vc.data(), nbasis, beta, cVc.data(), nbasis);
#endif
    return cVc;
}

inline double get_real(const std::complex<double>& c)
{
    return c.real();
}

inline double get_real(const double& d)
{
    return d;
}

template <typename T>
double all_band_energy(const int ik, const std::vector<T>& mat_mo, const Parallel_2D& p2d, const ModuleBase::matrix& wg)
{
    double e = 0.0;
    for (int i = 0; i < p2d.get_row_size(); ++i)
    {
        for (int j = 0; j < p2d.get_col_size(); ++j)
        {
            if (p2d.local2global_row(i) == p2d.local2global_col(j))
            {
                e += get_real(mat_mo[j * p2d.get_row_size() + i]) * wg(ik, p2d.local2global_row(i));
            }
        }
    }
    Parallel_Reduce::reduce_all(e);
    return e;
}

template <typename T>
std::vector<double> orbital_energy(const int ik, const int nbands, const std::vector<T>& mat_mo, const Parallel_2D& p2d)
{
#ifdef __DEBUG
    assert(nbands >= 0);
#endif
    std::vector<double> e(nbands, 0.0);
    for (int i = 0; i < nbands; ++i)
    {
        if (p2d.in_this_processor(i, i))
        {
            const int index = p2d.global2local_col(i) * p2d.get_row_size() + p2d.global2local_row(i);
            e[i] = get_real(mat_mo[index]);
        }
    }
    Parallel_Reduce::reduce_all(e.data(), nbands);
    return e;
}

// mohan update 2024-04-01
template <typename T>
void set_gint_pointer(Gint_Gamma& gint_gamma, Gint_k& gint_k, typename TGint<T>::type*& gint);

// mohan update 2024-04-01
template <>
void set_gint_pointer<double>(Gint_Gamma& gint_gamma, Gint_k& gint_k, typename TGint<double>::type*& gint)
{
    gint = &gint_gamma;
}

// mohan update 2024-04-01
template <>
void set_gint_pointer<std::complex<double>>(Gint_Gamma& gint_gamma,
                                            Gint_k& gint_k,
                                            typename TGint<std::complex<double>>::type*& gint)
{
    gint = &gint_k;
}

inline void write_orb_energy(const K_Vectors& kv,
    const int nspin0, const int nbands,
    const std::vector<std::vector<double>>& e_orb,
    const std::string& term, const std::string& label, const bool app = false)
{
    assert(e_orb.size() == kv.get_nks());
    const int nk = kv.get_nks() / nspin0;
    std::ofstream ofs;
    ofs.open(PARAM.globalv.global_out_dir + term + "_" + (label == "" ? "out.dat" : label + "_out.dat"),
        app ? std::ios::app : std::ios::out);
    ofs << nk << "\n" << nspin0 << "\n" << nbands << "\n";
    ofs << std::scientific << std::setprecision(16);
    for (int ik = 0; ik < nk; ++ik)
    {
        for (int is = 0; is < nspin0; ++is)
        {
            for (auto e : e_orb[is * nk + ik])
            { // Hartree and eV
                ofs << e / 2. << "\t" << e * ModuleBase::Ry_to_eV << "\n";
            }
        }
    }
}

/// @brief  write the Vxc matrix in KS orbital representation, usefull for GW calculation
/// including terms: local/semi-local XC, EXX, DFTU
template <typename TK, typename TR>
void write_Vxc(const int nspin,
               const int nbasis,
               const int drank,
               const Parallel_Orbitals* pv,
               const psi::Psi<TK>& psi,
               const UnitCell& ucell,
               Structure_Factor& sf,
               surchem& solvent,
               const ModulePW::PW_Basis& rho_basis,
               const ModulePW::PW_Basis& rhod_basis,
               const ModuleBase::matrix& vloc,
               const Charge& chg,
               Gint_Gamma& gint_gamma, // mohan add 2024-04-01
               Gint_k& gint_k,         // mohan add 2024-04-01
               const K_Vectors& kv,
               const std::vector<double>& orb_cutoff,
               const ModuleBase::matrix& wg,
               Grid_Driver& gd
#ifdef __EXX
               ,
               std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd = nullptr,
               std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc = nullptr
#endif
)
{
    ModuleBase::TITLE("ModuleIO", "write_Vxc");
    int nbands = wg.nc;
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

    // 2. allocate AO-matrix
    // R (the number of hR: 1 for nspin=1, 4; 2 for nspin=2)
    int nspin0 = (nspin == 2) ? 2 : 1;
    std::vector<hamilt::HContainer<TR>> vxcs_R_ao(nspin0, hamilt::HContainer<TR>(pv));
    for (int is = 0; is < nspin0; ++is) {
        vxcs_R_ao[is].set_zero();
        if (std::is_same<TK, double>::value) { vxcs_R_ao[is].fix_gamma(); }
    }
    // k (size for each k-point)
    hamilt::HS_Matrix_K<TK> vxc_k_ao(pv, 1); // only hk is needed, sk is skipped

    // 3. allocate operators and contribute HR
    // op (corresponding to hR)
    typename TGint<TK>::type* gint = nullptr;

    set_gint_pointer<TK>(gint_gamma, gint_k, gint);

    std::vector<hamilt::Veff<hamilt::OperatorLCAO<TK, TR>>*> vxcs_op_ao(nspin0);
    for (int is = 0; is < nspin0; ++is)
    {
        vxcs_op_ao[is] = new hamilt::Veff<hamilt::OperatorLCAO<TK, TR>>(gint,
            &vxc_k_ao, kv.kvec_d, potxc, &vxcs_R_ao[is], &ucell, orb_cutoff, &gd, nspin);

        vxcs_op_ao[is]->contributeHR();
    }
    std::vector<std::vector<double>> e_orb_locxc; // orbital energy (local XC)
    std::vector<std::vector<double>> e_orb_tot;   // orbital energy (total)
#ifdef __EXX
    hamilt::OperatorEXX<hamilt::OperatorLCAO<TK, TR>> vexx_op_ao(&vxc_k_ao,
        &vxcs_R_ao[0],ucell,/*for paraV*/ kv, Hexxd, Hexxc, hamilt::Add_Hexx_Type::k);
    hamilt::HS_Matrix_K<TK> vexxonly_k_ao(pv, 1); // only hk is needed, sk is skipped
    hamilt::OperatorEXX<hamilt::OperatorLCAO<TK, TR>> vexxonly_op_ao(&vexxonly_k_ao,
        &vxcs_R_ao[0],ucell,/*for paraV*/ kv, Hexxd, Hexxc, hamilt::Add_Hexx_Type::k);
    std::vector<std::vector<double>> e_orb_exx; // orbital energy (EXX)
#endif
    hamilt::OperatorDFTU<hamilt::OperatorLCAO<TK, TR>> vdftu_op_ao(&vxc_k_ao, kv.kvec_d, nullptr, kv.isk);

    // 4. calculate and write the MO-matrix Exc
    Parallel_2D p2d;
    set_para2d_MO(*pv, nbands, p2d);

    // ======test=======
    // double total_energy = 0.0;
    // double exx_energy = 0.0;
    // ======test=======
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        vxc_k_ao.set_zero_hk();
        int is = kv.isk[ik];
        dynamic_cast<hamilt::OperatorLCAO<TK, TR>*>(vxcs_op_ao[is])->contributeHk(ik);
        const std::vector<TK>& vlocxc_k_mo = cVc(vxc_k_ao.get_hk(), &psi(ik, 0, 0), nbasis, nbands, *pv, p2d);

#ifdef __EXX
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            e_orb_locxc.emplace_back(orbital_energy(ik, nbands, vlocxc_k_mo, p2d));
            ModuleBase::GlobalFunc::ZEROS(vexxonly_k_ao.get_hk(), pv->nloc);
            vexx_op_ao.contributeHk(ik);
            vexxonly_op_ao.contributeHk(ik);
            std::vector<TK> vexx_k_mo = cVc(vexxonly_k_ao.get_hk(), &psi(ik, 0, 0), nbasis, nbands, *pv, p2d);
            e_orb_exx.emplace_back(orbital_energy(ik, nbands, vexx_k_mo, p2d));
            // ======test=======
            // exx_energy += all_band_energy(ik, vexx_k_mo, p2d, wg);
            // ======test=======
        }
#endif
        if (PARAM.inp.dft_plus_u)
        {
            vdftu_op_ao.contributeHk(ik);
        }
        const std::vector<TK>& vxc_tot_k_mo = cVc(vxc_k_ao.get_hk(), &psi(ik, 0, 0), nbasis, nbands, *pv, p2d);
        e_orb_tot.emplace_back(orbital_energy(ik, nbands, vxc_tot_k_mo, p2d));
        // write
        ModuleIO::save_mat(-1,
                           vxc_tot_k_mo.data(),
                           nbands,
                           false /*binary*/,
                           PARAM.inp.out_ndigits,
                           true /*triangle*/,
                           false /*append*/,
                           "Vxc",
                           "k-" + std::to_string(ik),
                           p2d,
                           drank);
        // ======test=======
        // total_energy += all_band_energy(ik, vxc_tot_k_mo, p2d, wg);
        // ======test=======
    }
    // ======test=======
    // total_energy -= 0.5 * exx_energy;
    // std::cout << "total energy: " << total_energy << std::endl;
    // std::cout << "etxc: " << etxc << std::endl;
    // std::cout << "vtxc_cal: " << total_energy - 0.5 * exx_energy << std::endl;
    // std::cout << "vtxc_ref: " << vtxc << std::endl;
    // std::cout << "exx_energy: " << 0.5 * exx_energy << std::endl;
    // ======test=======
    delete potxc;
    for (int is = 0; is < nspin0; ++is)
    {
        delete vxcs_op_ao[is];
    }

    if (GlobalV::MY_RANK == 0)
    {
        write_orb_energy(kv, nspin0, nbands, e_orb_tot, "vxc", "");
#ifdef __EXX
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            write_orb_energy(kv, nspin0, nbands, e_orb_locxc, "vxc", "local");
            write_orb_energy(kv, nspin0, nbands, e_orb_exx, "vxc", "exx");
        }
#endif
    }
}
} // namespace ModuleIO
#endif
