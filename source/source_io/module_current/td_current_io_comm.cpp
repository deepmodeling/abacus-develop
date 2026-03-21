#include "td_current_io.h"
#include "source_base/timer.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_rt/td_info.h"
#include "source_base/libm/libm.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_threading.h"
#include "source_base/vector3.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_estate/module_pot/H_TDDFT_pw.h"
#include "source_lcao/LCAO_domain.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_rt/td_folding.h"
#ifdef __EXX
#include "source_lcao/module_operator_lcao/op_exx_lcao.h"
#include "source_lcao/module_ri/Exx_LRI.h"
#endif
#ifdef __LCAO
template <typename TR, typename TA>
void ModuleIO::init_from_hR(const hamilt::HContainer<TR>* hR, hamilt::HContainer<TA>* aimR)
{
    ModuleBase::TITLE("ModuleIO", "init_from_hR");
    ModuleBase::timer::tick("ModuleIO", "init_from_hR");
    for (int i = 0; i < hR->size_atom_pairs(); i++)
    {
        hamilt::AtomPair<TR> atom_ij = hR->get_atom_pair(i);
        const int iat1 = atom_ij.get_atom_i();
        const int iat2 = atom_ij.get_atom_j();
        for (int iR = 0; iR < atom_ij.get_R_size(); iR++)
        {
            const ModuleBase::Vector3<int> r_index = atom_ij.get_R_index(iR);
            hamilt::AtomPair<TA> atom_ij_ta(iat1, iat2, r_index, hR->get_paraV());
            aimR->insert_pair(atom_ij_ta);
        }
    }
    aimR->allocate(nullptr, true);

    ModuleBase::timer::tick("ModuleIO", "init_from_hR");
}
void ModuleIO::init_from_adj(const UnitCell& ucell,
                             const Grid_Driver& GridD,
                             const LCAO_Orbitals& orb,
                             const Parallel_Orbitals* pv,
                             std::vector<AdjacentAtomInfo>& adjs_all,
                             ModuleBase::Vector3<hamilt::HContainer<double>*>& rR)
{
    ModuleBase::TITLE("ModuleIOTD_mixing_pot", "init_from_adj");
    ModuleBase::timer::tick("ModuleIO", "init_from_adj");

    auto orb_cutoff_ = orb.cutoffs();
    adjs_all.clear();
    adjs_all.reserve(ucell.nat);

    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        auto tau1 = ucell.get_tau(iat1);
        int T1, I1;
        ucell.iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        GridD.Find_atom(ucell, tau1, T1, I1, &adjs);
        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T2 = adjs.ntype[ad1];
            const int I2 = adjs.natom[ad1];
            const int iat2 = ucell.itia2iat(T2, I2);
            if (pv->get_row_size(iat1) <= 0 || pv->get_col_size(iat2) <= 0)
            {
                continue;
            }
            const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad1];
            // choose the real adjacent atoms
            // Note: the distance of atoms should less than the cutoff radius,
            // When equal, the theoretical value of matrix element is zero,
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (ucell.cal_dtau(iat1, iat2, R_index2).norm() * ucell.lat0
                < orb_cutoff_[T1] + orb_cutoff_[T2])
            {
                is_adj[ad1] = true;
            }
        }
        filter_adjs(is_adj, adjs);
        adjs_all.push_back(adjs);
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            int iat2 = ucell.itia2iat(T2, I2);
            ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            hamilt::AtomPair<double> tmp(iat1, iat2, R_index, pv);
            for (size_t i_alpha = 0; i_alpha != 3; ++i_alpha)
            {
                rR[i_alpha]->insert_pair(tmp);
            }
        }
    }
    // allocate the memory of BaseMatrix in HR, and set the new values to zero
    for (size_t i_alpha = 0; i_alpha != 3; ++i_alpha)
    {
        rR[i_alpha]->allocate(nullptr, true);
    }
    ModuleBase::timer::tick("ModuleIO", "init_from_adj");
}

void ModuleIO::set_rR_from_hR(const UnitCell& ucell,
                              const Grid_Driver& GridD,
                              const LCAO_Orbitals& orb,
                              const Parallel_Orbitals* pv,
                              cal_r_overlap_R& r_calculator,
                              const hamilt::HContainer<std::complex<double>>* hR,
                              ModuleBase::Vector3<hamilt::HContainer<double>*>& rR)
{
    ModuleBase::TITLE("ModuleIO", "set_rR_from_hR");
    ModuleBase::timer::tick("ModuleIO", "set_rR_from_hR");

    // init
    std::vector<AdjacentAtomInfo> adjs_all;
    init_from_adj(ucell, GridD, orb, pv, adjs_all, rR);

    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        auto tau1 = ucell.get_tau(iat1);
        int T1, I1;
        ucell.iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo& adjs = adjs_all[iat1];
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            const int iat2 = ucell.itia2iat(T2, I2);
            const ModuleBase::Vector3<int>& r_index = adjs.box[ad];
            ModuleBase::Vector3<double> dtau = ucell.cal_dtau(iat1, iat2, r_index);

            Atom& atom1 = ucell.atoms[T1];
            Atom& atom2 = ucell.atoms[T2];
            const int npol = ucell.get_npol();

            const int* iw2l1 = atom1.iw2l.data();
            const int* iw2n1 = atom1.iw2n.data();
            const int* iw2m1 = atom1.iw2m.data();
            const int* iw2l2 = atom2.iw2l.data();
            const int* iw2n2 = atom2.iw2n.data();
            const int* iw2m2 = atom2.iw2m.data();

            auto row_indexes = pv->get_indexes_row(iat1);
            auto col_indexes = pv->get_indexes_col(iat2);

            const ModuleBase::Vector3<double>& tau1 = ucell.get_tau(iat1);
            // std::cout << "tau1: " << tau1 << " tau2: " << GlobalC::ucell.get_tau(iat2) << " r_index: " << r_index
            //           << std::endl;
            const ModuleBase::Vector3<double> tau2 = tau1 + dtau;
            for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
            {
                const int iw1 = row_indexes[iw1l] / npol;
                const int L1 = iw2l1[iw1];
                const int N1 = iw2n1[iw1];
                const int m1 = iw2m1[iw1];

                for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
                {
                    const int iw2 = col_indexes[iw2l] / npol;
                    const int L2 = iw2l2[iw2];
                    const int N2 = iw2n2[iw2];
                    const int m2 = iw2m2[iw2];

                    // std::cout<<"L1: "<<L1<<" L2: "<<L2<<" N1: "<<N1<<" N2: "<<N2<<" m1: "<<m1<<" m2:
                    // "<<m2<<std::endl;
                    ModuleBase::Vector3<double> tmp_r = r_calculator.get_psi_r_psi(tau1 * ucell.lat0,
                                                                                   T1,
                                                                                   L1,
                                                                                   m1,
                                                                                   N1,
                                                                                   tau2 * ucell.lat0,
                                                                                   T2,
                                                                                   L2,
                                                                                   m2,
                                                                                   N2);
                    for (size_t i_alpha = 0; i_alpha != 3; ++i_alpha)
                    {
                        hamilt::BaseMatrix<double>* HlocR = rR[i_alpha]->find_matrix(iat1, iat2, r_index);
                        if(HlocR!=nullptr)HlocR->add_element(iw1, iw2, tmp_r[i_alpha]);
                        // if (i_alpha == 2)
                        // {
                        //     std::cout << "iw1: " << iw1 << " iw2: " << iw2 << " i_alpha: " << i_alpha
                        //               << " tmp_r: " << tmp_r[i_alpha] << std::endl;
                        // }
                    }
                }
            }
        }
    }
    ModuleBase::TITLE("ModuleIO", "set_rR_from_sR");
}
template <typename TR>
void ModuleIO::sum_HR(
    const UnitCell& ucell,
    const Parallel_Orbitals& pv,
    const K_Vectors& kv,
    const hamilt::HContainer<TR>* hR,
    hamilt::HContainer<std::complex<double>>* full_hR,
    const Exx_NAO<std::complex<double>>& exx_nao
)
{
    ModuleBase::TITLE("ModuleIO", "sum_HR");
    ModuleBase::timer::tick("ModuleIO", "sum_HR");

    // init complex full_hR
    init_from_hR(hR, full_hR);
#ifdef __EXX
    const bool use_cell_nearest = (ModuleBase::Vector3<double>(std::fmod(kv.get_koffset(0), 1.0),
                                                               std::fmod(kv.get_koffset(1), 1.0),
                                                               std::fmod(kv.get_koffset(2), 1.0)).norm()< 1e-10);
    RI::Cell_Nearest<int, int, 3, double, 3> cell_nearest;
    // reallocate full_hR for BvK used in EXX
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        const std::array<int, 3> Rs_period = {kv.nmp[0], kv.nmp[1], kv.nmp[2]};
        if (use_cell_nearest)
        {
            // set cell_nearest
            std::map<int, std::array<double, 3>> atoms_pos;
            for (int iat = 0; iat < ucell.nat; ++iat)
            {
                atoms_pos[iat] = RI_Util::Vector3_to_array3(
                    ucell.atoms[ucell.iat2it[iat]].tau[ucell.iat2ia[iat]]);
            }
            const std::array<std::array<double, 3>, 3> latvec = {RI_Util::Vector3_to_array3(ucell.a1),
                                                                 RI_Util::Vector3_to_array3(ucell.a2),
                                                                 RI_Util::Vector3_to_array3(ucell.a3)};
            cell_nearest.init(atoms_pos, latvec, Rs_period);
            hamilt::reallocate_hcontainer(ucell.nat, full_hR, Rs_period, &cell_nearest);
        }
        else
            hamilt::reallocate_hcontainer(ucell.nat, full_hR, Rs_period);
    }
#endif
    // add other hR
    add_HR(hR, full_hR);
    // add velocity complex hR
    if (PARAM.inp.td_stype==1)
    {
        if (TD_info::td_vel_op == nullptr)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::write_current", "velocity gauge infos is null!");
        }
        const hamilt::HContainer<std::complex<double>>* velocity_hR = TD_info::td_vel_op->get_velocity_HR_pointer();
        add_HR(velocity_hR, full_hR);
    }
#ifdef __EXX
    // add HexxR to complex full_hR
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        for (size_t is = 0; is != PARAM.inp.nspin; ++is)
        {
            if (use_cell_nearest)
                RI_2D_Comm::add_HexxR(is, GlobalC::exx_info.info_global.hybrid_alpha, exx_nao.exc->get_Hexxs(), pv, PARAM.globalv.npol, *full_hR, &cell_nearest);
            else
                RI_2D_Comm::add_HexxR(is, GlobalC::exx_info.info_global.hybrid_alpha, exx_nao.exc->get_Hexxs(), pv, PARAM.globalv.npol, *full_hR, nullptr);
        }
    }
#endif

    ModuleBase::timer::tick("ModuleIO", "sum_HR");
}

template <typename Tadd, typename Tfull>
void ModuleIO::add_HR(const hamilt::HContainer<Tadd>* hR, hamilt::HContainer<Tfull>* full_hR)
{
    ModuleBase::TITLE("ModuleIO", "add_HR");
    ModuleBase::timer::tick("ModuleIO", "add_HR");

    for (int ipair = 0; ipair < hR->size_atom_pairs(); ++ipair)
    {
        hamilt::AtomPair<Tadd> atom_ij = hR->get_atom_pair(ipair);
        const int iat1 = atom_ij.get_atom_i();
        const int iat2 = atom_ij.get_atom_j();
        // loop R-index
        for (int iR = 0; iR < atom_ij.get_R_size(); iR++)
        {
            const ModuleBase::Vector3<int> r_index = atom_ij.get_R_index(iR);
            hamilt::BaseMatrix<Tfull>* full_HlocR = full_hR->find_matrix(iat1, iat2, r_index.x, r_index.y, r_index.z);
            const hamilt::BaseMatrix<Tadd>* HlocR = hR->find_matrix(iat1, iat2, r_index.x, r_index.y, r_index.z);

            if (full_HlocR == nullptr || HlocR == nullptr)
                ModuleBase::WARNING_QUIT("ModuleIO::add_HR", "HR cannot be nullptr!");

            for (int i = 0; i < atom_ij.get_row_size(); ++i)
            {
                for (int j = 0; j < atom_ij.get_col_size(); ++j)
                {
                    Tadd v = HlocR->get_value(i, j);
                    full_HlocR->add_element(i, j, Tfull(v));
                }
            }
        }
    }

    ModuleBase::timer::tick("ModuleIO", "add_HR");
}

// for molecule, if vacuum size is small, the number of R of Hs is smaller than SR
// which may lead to some errors
template <typename TR>
void ModuleIO::cal_velocity_basis_k(const UnitCell& ucell,
                                    const LCAO_Orbitals& orb,
                                    const Parallel_Orbitals* pv,
                                    const K_Vectors& kv,
                                    const ModuleBase::Vector3<hamilt::HContainer<double>*>& rR,
                                    const hamilt::HContainer<TR>& sR,
                                    const hamilt::HContainer<std::complex<double>>& hR,
                                    std::vector<ModuleBase::Vector3<std::complex<double>*>>& velocity_basis_k)
{
    ModuleBase::TITLE("ModuleIO", "cal_velocity_basis_k");
    ModuleBase::timer::tick("ModuleIO", "cal_velocity_basis_k");
#ifdef __MPI
    const int nlocal = PARAM.globalv.nlocal;
    const char N_char = 'N';
    const std::complex<double> one_imag = ModuleBase::IMAG_UNIT;
    const std::complex<double> neg_one_imag = ModuleBase::NEG_IMAG_UNIT;
    const std::complex<double> one_real = ModuleBase::ONE;
    const std::complex<double> neg_one_real = ModuleBase::NEG_ONE;
    const std::complex<double> zero_complex = ModuleBase::ZERO;
    std::cout<<pv->nloc<<std::endl;
    std::complex<double>* hk = new std::complex<double>[pv->nloc];
    std::complex<double>* sk = new std::complex<double>[pv->nloc];
    std::complex<double>* partial_hk = new std::complex<double>[pv->nloc];
    std::complex<double>* partial_sk = new std::complex<double>[pv->nloc];
    std::complex<double>* rk = new std::complex<double>[pv->nloc];
    std::complex<double>* h_is = new std::complex<double>[pv->nloc];
    std::complex<double>* h_is_r = new std::complex<double>[pv->nloc];
    std::complex<double>* r_is = new std::complex<double>[pv->nloc];
    std::complex<double>* r_is_h = new std::complex<double>[pv->nloc];
    std::complex<double>* h_is_ps = new std::complex<double>[pv->nloc];

    // for (size_t i_alpha = 0; i_alpha != 3; ++i_alpha)
    // {
    //     for (int i = 0; i < hR.size_atom_pairs(); ++i)
    //     {
    //         hamilt::AtomPair<double>& tmp = rR[i_alpha]->get_atom_pair(i);
    //         std::cout<<"cal_velocity_basis_k: "<<rR[i_alpha]->size_atom_pairs()<<" R_size:
    //         "<<tmp.get_R_size()<<std::endl; for(int ir = 0;ir < tmp.get_R_size(); ++ir )
    //         {
    //             const ModuleBase::Vector3<int> r_index = tmp.get_R_index(ir);
    //             std::cout<<"r_index: "<<r_index<<std::endl;
    //         }
    //     }
    // }

    for (size_t ik = 0; ik != kv.get_nks(); ++ik)
    {
        // set H(k), S(k)
        // 1.1 set H(k)
        ModuleBase::GlobalFunc::ZEROS(hk, pv->nloc);
        const int nrow = pv->get_row_size();
        if (elecstate::H_TDDFT_pw::stype == 2)
            module_rt::folding_HR_td(ucell, hR, hk, kv.kvec_d[ik], TD_info::cart_At, nrow, 1);
        else
            hamilt::folding_HR(hR, hk, kv.kvec_d[ik], nrow, 1);
        // 1.2 set S(k)
        ModuleBase::GlobalFunc::ZEROS(sk, pv->nloc);
        if (elecstate::H_TDDFT_pw::stype == 2)
             module_rt::folding_HR_td(ucell, sR, sk, kv.kvec_d[ik], TD_info::cart_At, nrow, 1);
        else
            hamilt::folding_HR(sR, sk, kv.kvec_d[ik], nrow, 1);
        // for (int ir = 0; ir < pv->nrow; ir++)
        // {
        //     const int iwt1 = pv->local2global_row(ir);
        //     const int iat1 = GlobalC::ucell.iwt2iat[iwt1];
        //     for (int ic = 0; ic < pv->ncol; ic++)
        //     {
        //         const int iwt2 = pv->local2global_col(ic);
        //         const int iat2 = GlobalC::ucell.iwt2iat[iwt2];
        //         const int irc = ic * pv->nrow + ir;
        //         std::cout << "ik: " << ik << " iat1:" << iat1 << " iat2:" << iat2 << " iwt1: " << iwt1
        //                   << " iwt2: " << iwt2 << " hk: " << hk[irc] << std::endl;
        //     }
        // }
        // 2. set inverse S(k) -> sk will be changed to sk_inv
        int* ipiv = new int[pv->nloc];
        int info = 0;
        // 2.1 compute ipiv
        ScalapackConnector::getrf(nlocal, nlocal, sk, 1, 1, pv->desc, ipiv, &info);
        int lwork = -1;
        int liwotk = -1;
        std::vector<std::complex<double>> work(1, 0);
        std::vector<int> iwork(1, 0);
        // 2.2 compute work
        ScalapackConnector::getri(nlocal, sk, 1, 1, pv->desc, ipiv, work.data(), &lwork, iwork.data(), &liwotk, &info);
        lwork = work[0].real();
        work.resize(lwork, 0);
        liwotk = iwork[0];
        iwork.resize(liwotk, 0);
        // 2.3 compute inverse matrix of Sk
        ScalapackConnector::getri(nlocal,
                                  sk, // return sk^-1
                                  1,
                                  1,
                                  pv->desc,
                                  ipiv,
                                  work.data(),
                                  &lwork,
                                  iwork.data(),
                                  &liwotk,
                                  &info);
        delete[] ipiv;
        assert(0 == info);
        for (size_t i_alpha = 0; i_alpha != 3; ++i_alpha)
        {
            // 3. set partial_H(k), partial_S(k) and r(k)
            // 3.1 set partial_H(k)
            ModuleBase::GlobalFunc::ZEROS(partial_hk, pv->nloc);
            if (elecstate::H_TDDFT_pw::stype == 2)
                 module_rt::folding_partial_HR_td(ucell, hR, partial_hk, kv.kvec_d[ik], TD_info::cart_At, i_alpha, nrow, 1);
            else
                 module_rt::folding_partial_HR(ucell, hR, partial_hk, kv.kvec_d[ik], i_alpha, nrow, 1);
            // 3.2 set partial S(k)
            ModuleBase::GlobalFunc::ZEROS(partial_sk, pv->nloc);
            if (elecstate::H_TDDFT_pw::stype == 2)
                 module_rt::folding_partial_HR_td(ucell, sR, partial_sk, kv.kvec_d[ik], TD_info::cart_At, i_alpha, nrow, 1);
            else
                 module_rt::folding_partial_HR(ucell, sR, partial_sk, kv.kvec_d[ik], i_alpha, nrow, 1);
            //  if(i_alpha == 2)
            // {
            //     for(int ir=0;ir< pv->nrow; ir++)
            //     {
            //         const int iwt1 = pv->local2global_row(ir);
            //         const int iat1 = GlobalC::ucell.iwt2iat[iwt1];
            //         for(int ic=0;ic< pv->ncol; ic++)
            //         {
            //             const int iwt2 = pv->local2global_col(ic);
            //             const int iat2 = GlobalC::ucell.iwt2iat[iwt2];
            //             const int irc=ic*pv->nrow + ir;
            //             std::cout<<"ik: "<<ik<<" i_alpha: "<<i_alpha<<" iat1:"<<iat1<<" iat2:"<<iat2<<" iwt1:
            //             "<<iwt1<<" iwt2: "<<iwt2<<" partial_sk: "<<partial_sk[irc]<<std::endl;
            //         }
            //     }
            // }
            // 3.3 set r(k)
            // std::cout << "set r(k): " << "i_alpha: " << i_alpha << std::endl;
            ModuleBase::GlobalFunc::ZEROS(rk, pv->nloc);
            // folding_rR(rR[i_alpha], partial_sk, rk, pv, kv.kvec_d[ik], nrow, 1);
            if (elecstate::H_TDDFT_pw::stype == 2)
                module_rt::folding_HR_td(ucell, *rR[i_alpha], rk, kv.kvec_d[ik], TD_info::cart_At, nrow, 1);
            else
                hamilt::folding_HR(*rR[i_alpha], rk, kv.kvec_d[ik], nrow, 1); // set r(k)
            // if (i_alpha == 2)
            // {
            //     std::cout << "ik: " << ik << " i_alpha: " << i_alpha << std::endl;
            //     for (int ir = 0; ir < pv->nrow; ir++)
            //     {
            //         const int iwt1 = pv->local2global_row(ir);
            //         const int iat1 = GlobalC::ucell.iwt2iat[iwt1];
            //         for (int ic = 0; ic < pv->ncol; ic++)
            //         {
            //             const int iwt2 = pv->local2global_col(ic);
            //             const int iat2 = GlobalC::ucell.iwt2iat[iwt2];
            //             const int irc = ic * pv->nrow + ir;
            //             std::cout << " iat1: " << iat1 << " iat2: " << iat2 << " iw1: " <<
            //             GlobalC::ucell.iwt2iw[iwt1]
            //                       << " iw2: " << GlobalC::ucell.iwt2iw[iwt2] << " rk: " << rk[irc] << std::endl;
            //         }
            //     }
            // }
            // 4. calculate <\vu,k|v_a|\mu,k> = partial_Hk + IMAG_UNIT * (Hk * Sk_inv * rk) - IMAG_UNIT * (rk * Sk_inv *
            // Hk) - Hk * Sk_inv * partial_Sk
            // 4.1.1 Hk * Sk_inv (note 2.)
            ModuleBase::GlobalFunc::ZEROS(h_is, pv->nloc);
            ScalapackConnector::gemm(N_char,
                                     N_char,
                                     nlocal,
                                     nlocal,
                                     nlocal,
                                     one_real,
                                     hk,
                                     1,
                                     1,
                                     pv->desc,
                                     sk,
                                     1,
                                     1,
                                     pv->desc,
                                     zero_complex,
                                     h_is,
                                     1,
                                     1,
                                     pv->desc);
            // 4.1.2 (Hk * Sk_inv) * rk
            ModuleBase::GlobalFunc::ZEROS(h_is_r, pv->nloc);
            ScalapackConnector::gemm(N_char,
                                     N_char,
                                     nlocal,
                                     nlocal,
                                     nlocal,
                                     one_real,
                                     h_is,
                                     1,
                                     1,
                                     pv->desc,
                                     rk,
                                     1,
                                     1,
                                     pv->desc,
                                     zero_complex,
                                     h_is_r,
                                     1,
                                     1,
                                     pv->desc);
            // 4.2.1 rk * Sk_inv (note 2.)
            ModuleBase::GlobalFunc::ZEROS(r_is, pv->nloc);
            ScalapackConnector::gemm(N_char,
                                     N_char,
                                     nlocal,
                                     nlocal,
                                     nlocal,
                                     one_real,
                                     rk,
                                     1,
                                     1,
                                     pv->desc,
                                     sk,
                                     1,
                                     1,
                                     pv->desc,
                                     zero_complex,
                                     r_is,
                                     1,
                                     1,
                                     pv->desc);
            // 4.2.2 (rk * Sk_inv) * Hk
            ModuleBase::GlobalFunc::ZEROS(r_is_h, pv->nloc);
            ScalapackConnector::gemm(N_char,
                                     N_char,
                                     nlocal,
                                     nlocal,
                                     nlocal,
                                     one_real,
                                     r_is,
                                     1,
                                     1,
                                     pv->desc,
                                     hk,
                                     1,
                                     1,
                                     pv->desc,
                                     zero_complex,
                                     r_is_h,
                                     1,
                                     1,
                                     pv->desc);
            // 4.3.1 (Hk * Sk_inv) * partial_Sk
            ModuleBase::GlobalFunc::ZEROS(h_is_ps, pv->nloc);
            ScalapackConnector::gemm(N_char,
                                     N_char,
                                     nlocal,
                                     nlocal,
                                     nlocal,
                                     one_real,
                                     h_is,
                                     1,
                                     1,
                                     pv->desc,
                                     partial_sk,
                                     1,
                                     1,
                                     pv->desc,
                                     zero_complex,
                                     h_is_ps,
                                     1,
                                     1,
                                     pv->desc);
            // 4.4 h_is_r will be changed to partial_Hk + IMAG_UNIT * (Hk * Sk_inv * rk)
            ScalapackConnector::geadd('N',
                                      nlocal,
                                      nlocal,
                                      one_real,
                                      partial_hk,
                                      1,
                                      1,
                                      pv->desc,
                                      one_imag,
                                      h_is_r,
                                      1,
                                      1,
                                      pv->desc);
            // 4.5 r_is_h will be changed to h_is_r - IMAG_UNIT * (rk * Sk_inv * Hk)
            ScalapackConnector::geadd('N',
                                      nlocal,
                                      nlocal,
                                      one_real,
                                      h_is_r,
                                      1,
                                      1,
                                      pv->desc,
                                      neg_one_imag,
                                      r_is_h,
                                      1,
                                      1,
                                      pv->desc);
            // 4.6 h_is_ps will be changed to r_is_h - Hk * Sk_inv * partial_Sk
            ScalapackConnector::geadd('N',
                                      nlocal,
                                      nlocal,
                                      one_real,
                                      r_is_h,
                                      1,
                                      1,
                                      pv->desc,
                                      neg_one_real,
                                      h_is_ps,
                                      1,
                                      1,
                                      pv->desc);
            // 5. copy h_is_ps to velocity_basis_k[ik][i_alpha]
            BlasConnector::copy(pv->nloc, h_is_ps, 1, velocity_basis_k[ik][i_alpha], 1);
            // if(i_alpha == 2)
            // {
            //     for(int ir=0;ir< pv->nrow; ir++)
            //     {
            //         const int iwt1 = pv->local2global_row(ir);
            //         const int iat1 = GlobalC::ucell.iwt2iat[iwt1];
            //         for(int ic=0;ic< pv->ncol; ic++)
            //         {
            //             const int iwt2 = pv->local2global_col(ic);
            //             const int iat2 = GlobalC::ucell.iwt2iat[iwt2];
            //             const int irc=ic*pv->nrow + ir;
            //             std::cout<<"ik: "<<ik<<" i_alpha: "<<i_alpha<<" iat1:"<<iat1<<" iat2:"<<iat2<<" v_basis_k:
            //             "<<velocity_basis_k[ik][i_alpha][irc]<<std::endl;
            //         }
            //     }
            // }
        }
    }

    delete[] hk;
    delete[] sk;
    delete[] partial_hk;
    delete[] partial_sk;
    delete[] rk;
    delete[] h_is;
    delete[] h_is_r;
    delete[] r_is;
    delete[] r_is_h;
    delete[] h_is_ps;
#endif //__MPI
    ModuleBase::timer::tick("ModuleIO", "cal_velocity_basis_k");
}

void ModuleIO::cal_velocity_matrix(const psi::Psi<std::complex<double>>* psi,
                                   const Parallel_Orbitals* pv,
                                   const K_Vectors& kv,
                                   const std::vector<ModuleBase::Vector3<std::complex<double>*>>& velocity_basis_k,
                                   std::vector<std::array<ModuleBase::ComplexMatrix, 3>>& velocity_k)
{
    ModuleBase::TITLE("ModuleIO", "cal_velocity_matrix");
    ModuleBase::timer::tick("ModuleIO", "cal_velocity_matrix");
#ifdef __MPI
    const char N_char = 'N';
    const char C_char = 'C';
    const std::complex<double> one_real = ModuleBase::ONE;
    const std::complex<double> zero_complex = ModuleBase::ZERO;
    const double zero_double = 0.0;
    const int nlocal = PARAM.globalv.nlocal;
    const int nbands = PARAM.inp.nbands;
    std::complex<double>* vk_c = new std::complex<double>[pv->ncol_bands * pv->nrow_bands]; // local one
    std::complex<double>* v_c = new std::complex<double>[pv->nloc_wfc];

    for (int ik = 0; ik < kv.get_nks(); ik++)
    {
        // 1. set C
        psi->fix_k(ik);
        // 2. set <\Psi_{n,\mu}|v_{\mu,\nu}|\Psi_{m,\nu}> = C^\dagger_{n,\mu} * v_{\mu,\nu} * C_{\nu,m}
        for (size_t i_alpha = 0; i_alpha != 3; ++i_alpha)
        {
            ModuleBase::GlobalFunc::ZEROS(vk_c, pv->ncol_bands * pv->nrow_bands);
            ModuleBase::GlobalFunc::ZEROS(v_c, pv->nloc_wfc);
            // v_c_{\mu,m} = v_{\mu,\nu} * C_{\nu,m}
            ScalapackConnector::gemm(N_char,
                                     N_char,
                                     nlocal,
                                     nbands,
                                     nlocal,
                                     one_real,
                                     velocity_basis_k[ik][i_alpha],
                                     1,
                                     1,
                                     pv->desc,
                                     psi[0].get_pointer(),
                                     1,
                                     1,
                                     pv->desc_wfc,
                                     zero_complex,
                                     v_c,
                                     1,
                                     1,
                                     pv->desc_wfc);
            // velocity_k_{n,m} = C^\dagger_{n,\mu} * v_c_{\mu,m}
            ScalapackConnector::gemm(C_char,
                                     N_char,
                                     nbands,
                                     nbands,
                                     nlocal,
                                     one_real,
                                     psi[0].get_pointer(),
                                     1,
                                     1,
                                     pv->desc_wfc,
                                     v_c,
                                     1,
                                     1,
                                     pv->desc_wfc,
                                     zero_complex,
                                     vk_c,
                                     1,
                                     1,
                                     pv->desc_Eij);

            for (int ir = 0; ir < PARAM.inp.nbands; ++ir)
            {
                // const int iwt1 = pv->local2global_row(ir);
                // const int iat1 = GlobalC::ucell.iwt2iat[iwt1];
                for (int ic = 0; ic < PARAM.inp.nbands; ++ic)
                {
                    const int irc = ic * pv->nrow + ir;
                    if (pv->in_this_processor(ir, ic))
                    {
                        // const int iwt2 = pv->local2global_col(ic);
                        // const int iat2 = GlobalC::ucell.iwt2iat[iwt2];
                        velocity_k[ik][i_alpha](ir, ic) = vk_c[irc];
                        // if (i_alpha == 0)
                        // {
                        //     std::cout<<"ik: "<<ik<<" i_alpha: "<<i_alpha<<" iat1:"<<iat1<<" iat2:"<<iat2<<" vk:
                        //     "<<velocity_basis_k[ik][i_alpha][irc]<<std::endl;
                        // }
                    }
                }
            }
        }
    }

    delete[] vk_c;
    delete[] v_c;
#endif //__MPI
    ModuleBase::timer::tick("ModuleIO", "cal_velocity_matrix");
}
template <typename TR>
void ModuleIO::cal_current_comm_k(const UnitCell& ucell,
                                  const Grid_Driver& GridD,
                                  const LCAO_Orbitals& orb,
                                  const Parallel_Orbitals* pv,
                                  const K_Vectors& kv,
                                  cal_r_overlap_R& r_calculator,
                                  const hamilt::HContainer<TR>& sR,
                                  const hamilt::HContainer<std::complex<double>>& hR,
                                  const psi::Psi<std::complex<double>>* psi,
                                  const elecstate::ElecState* pelec,
                                  std::vector<ModuleBase::Vector3<double>>& current_k)
{
    ModuleBase::TITLE("ModuleIO", "cal_current_exx");
    ModuleBase::timer::tick("ModuleIO", "cal_current_exx");

    const int nlocal = PARAM.globalv.nlocal;
    const int nbands = PARAM.inp.nbands;
    // init
    ModuleBase::Vector3<hamilt::HContainer<double>*> rR(nullptr, nullptr, nullptr);
    std::vector<ModuleBase::Vector3<std::complex<double>*>> velocity_basis_k;
    std::vector<std::array<ModuleBase::ComplexMatrix, 3>> velocity_k;
    velocity_basis_k.resize(kv.get_nks());
    velocity_k.resize(kv.get_nks());
    for (size_t i_alpha = 0; i_alpha != 3; ++i_alpha)
    {
        rR[i_alpha] = new hamilt::HContainer<double>(pv);
        for (int ik = 0; ik < kv.get_nks(); ik++)
        {
            velocity_basis_k[ik][i_alpha] = new std::complex<double>[pv->nloc];
            ModuleBase::GlobalFunc::ZEROS(velocity_basis_k[ik][i_alpha], pv->nloc);
            velocity_k[ik][i_alpha].create(nbands, nbands);
        }
    }
    // set rR
    set_rR_from_hR(ucell, GridD, orb, pv, r_calculator, &hR, rR);
    // set velocity_basis_k
    cal_velocity_basis_k(ucell, orb, pv, kv, rR, sR, hR, velocity_basis_k);
    // set velocity_k
    cal_velocity_matrix(psi, pv, kv, velocity_basis_k, velocity_k);

    // sum n and m for current_k
    for (size_t ik = 0; ik != kv.get_nks(); ++ik)
        for (size_t i_alpha = 0; i_alpha != 3; ++i_alpha)
        {
            for (size_t ib = 0; ib != PARAM.inp.nbands; ++ib)
                current_k[ik][i_alpha] -= pelec->wg(ik, ib) * velocity_k[ik][i_alpha](ib, ib).real() / 2.0; // for unit
        }
    for (size_t i_alpha = 0; i_alpha < 3; ++i_alpha)
    {
        delete rR[i_alpha];
        for (int ik = 0; ik < kv.get_nks(); ik++)
            delete[] velocity_basis_k[ik][i_alpha];
    }

    ModuleBase::TITLE("ModuleIO", "cal_current_exx");
}
template <typename TR>
void ModuleIO::write_current(
    const UnitCell& ucell,
    const Grid_Driver& GridD,
    const int istep,
    const psi::Psi<std::complex<double>>* psi,
    const elecstate::ElecState* pelec,
    const K_Vectors& kv,
    const Parallel_Orbitals* pv,
    const LCAO_Orbitals& orb,
    cal_r_overlap_R& r_calculator,
    const hamilt::HContainer<TR>* sR,
    const hamilt::HContainer<TR>* hR,
    const Exx_NAO<std::complex<double>>& exx_nao
)
{
    ModuleBase::TITLE("ModuleIO", "write_current");
    ModuleBase::timer::tick("ModuleIO", "write_current");
    double omega = ucell.omega;

    std::vector<ModuleBase::Vector3<double>> current_k;
    hamilt::HContainer<std::complex<double>>* full_hR;
    full_hR = new hamilt::HContainer<std::complex<double>>(pv);
    current_k.resize(kv.get_nks());
    sum_HR(ucell, *pv, kv, hR, full_hR, exx_nao);
    cal_current_comm_k(ucell, GridD, orb, pv, kv, r_calculator, *sR, *full_hR, psi, pelec, current_k);
    delete full_hR;

    int nspin0 = 1;
    if (PARAM.inp.nspin == 2)
    {
        nspin0 = 2;
    }
    for (int is = 0; is < nspin0; ++is)
    {
        for (int ik = 0; ik < kv.get_nks(); ik++)
        {
            if (is == kv.isk[ik])
            {
                if (GlobalV::MY_RANK == 0 && TD_info::out_current_k)
                {
                    std::string filename = PARAM.globalv.global_out_dir + "currents" + std::to_string(is)
                                           + "k" + std::to_string(ik) + "comm.txt";
                    std::ofstream fout;
                    fout.open(filename, std::ios::app);
                    fout << std::setprecision(16);
                    fout << std::scientific;
                    fout << istep << " " << current_k[ik][0] / omega << " " << current_k[ik][1] / omega << " "
                         << current_k[ik][2] / omega << std::endl;
                    fout.close();
                }
            }
        }
    }

    ModuleBase::Vector3<double> current_total;
    for (int dir = 0; dir < 3; dir++)
        for (int ik = 0; ik < kv.get_nks(); ik++)
            current_total[dir] += current_k[ik][dir];
    if (GlobalV::MY_RANK == 0)
    {
        std::string filename = PARAM.globalv.global_out_dir + "current_tot_comm.txt";
        std::ofstream fout;
        fout.open(filename, std::ios::app);
        fout << std::setprecision(16);
        fout << std::scientific;
        fout << istep << " " << current_total[0] / omega << " " << current_total[1] / omega << " "
             << current_total[2] / omega << std::endl;
        fout.close();
    }

    ModuleBase::timer::tick("ModuleIO", "write_current");
}
template 
void ModuleIO::write_current<double>(
        const UnitCell& ucell,
        const Grid_Driver& GridD,
        const int istep,
        const psi::Psi<std::complex<double>>* psi,
        const elecstate::ElecState* pelec,
        const K_Vectors& kv,
        const Parallel_Orbitals* pv,
        const LCAO_Orbitals& orb,
        cal_r_overlap_R& r_calculator,
        const hamilt::HContainer<double>* sR,
        const hamilt::HContainer<double>* hR,
        const Exx_NAO<std::complex<double>>& exx_nao
);

template 
void ModuleIO::write_current<std::complex<double>>(
        const UnitCell& ucell,
        const Grid_Driver& GridD,
        const int istep,
        const psi::Psi<std::complex<double>>* psi,
        const elecstate::ElecState* pelec,
        const K_Vectors& kv,
        const Parallel_Orbitals* pv,
        const LCAO_Orbitals& orb,
        cal_r_overlap_R& r_calculator,
        const hamilt::HContainer<std::complex<double>>* sR,
        const hamilt::HContainer<std::complex<double>>* hR,
        const Exx_NAO<std::complex<double>>& exx_nao
);
#endif //__LCAO
