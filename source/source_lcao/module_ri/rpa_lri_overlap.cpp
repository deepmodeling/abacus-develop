#include "exx_lri.h"
#include "rpa_lri.h"
#include "rpa_lri_detail.h"
#include "source_base/global_function.h"
#include "source_basis/module_ao/elem_basis_idx_orb.h"
#include "source_estate/elecstate_lcao.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_lcao/module_ri/module_exx_symmetry/symm_rotation.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::cal_large_Cs(const UnitCell& ucell, const LCAO_Orbitals& orb, const K_Vectors& kv)
{
    ModuleBase::TITLE("RPA_LRI", "cal_large_Cs");
    ModuleBase::timer::start("RPA_LRI", "cal_large_Cs");
    if (!exx_cut_coulomb)
    {
        Exx_Info_RI local_info = this->info;
        local_info.ccp_rmesh_times = this->ccp_rmesh_times_cut;
        Exx_LRI<double>* new_exx = new Exx_LRI<double>(local_info);
        exx_cut_coulomb.reset(new_exx);
    }
    exx_cut_coulomb->init_spencer(this->mpi_comm, ucell, kv, orb);
    ModuleBase::GlobalFunc::DONE(this->runtime.log, "exx_cut_coulomb->init");
    this->abfs = exx_cut_coulomb->abfs;
    this->MGT = exx_cut_coulomb->MGT;
    std::vector<TA> atoms(ucell.nat);
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        atoms[iat] = iat;
    }
    std::map<TA, TatomR> atoms_pos;
    for (int iat = 0; iat < ucell.nat; ++iat)
        atoms_pos[iat] = RI_Util::Vector3_to_array3(ucell.atoms[ucell.iat2it[iat]].tau[ucell.iat2ia[iat]]);
    const std::array<TatomR, Ndim> latvec = {RI_Util::Vector3_to_array3(ucell.a1),
                                             RI_Util::Vector3_to_array3(ucell.a2),
                                             RI_Util::Vector3_to_array3(ucell.a3)};
    const std::array<Tcell, Ndim> period = {p_kv->nmp[0], p_kv->nmp[1], p_kv->nmp[2]};
    this->exx_cut_coulomb->exx_lri.set_parallel(this->mpi_comm, atoms_pos, latvec, period);
    auto center2_obj_it = this->exx_cut_coulomb->exx_objs.find(Conv_Coulomb_Pot_K::Coulomb_Method::Center2);
    if (center2_obj_it == this->exx_cut_coulomb->exx_objs.end())
    {
        throw std::invalid_argument("RPA_LRI::cal_large_Cs expected a Center2 cut-Coulomb object after init_spencer.");
    }
    center2_obj_it->second.cv.set_orbitals(ucell,
                                           orb,
                                           this->lcaos,
                                           exx_cut_coulomb->abfs,
                                           center2_obj_it->second.abfs_ccp,
                                           this->info.kmesh_times,
                                           this->MGT, // get MGT from exx_cut_coulomb and used in `cal_abfs_overlap`
                                           true);

    const std::array<Tcell, Ndim> period_Vs
        = LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->ccp_rmesh_times_cut, ucell, orb_cutoff_);
    std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>> list_As_Vs
        = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Vs, 2, false);
    ModuleBase::GlobalFunc::DONE(this->runtime.log, "cal_large_Vs start");
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_cut_IJR
        = center2_obj_it->second.cv.cal_Vs(ucell, list_As_Vs.first, list_As_Vs.second[0], {{"writable_Vws", true}});
    ModuleBase::GlobalFunc::DONE(this->runtime.log, "cal_large_Vs end");

    const std::array<Tcell, Ndim> period_Cs = LRI_CV_Tools::cal_latvec_range<Tcell>(2, ucell, orb_cutoff_);
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>> list_As_Cs
        = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Cs, 2, false);
    std::pair<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>,
              std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, 3>>>>
        Cs_dCs = center2_obj_it->second.cv.cal_Cs_dCs(ucell,
                                                      list_As_Cs.first,
                                                      list_As_Cs.second[0],
                                                      {{"cal_dC", false},
                                                       {"writable_Cws", true},
                                                       {"writable_dCws", true},
                                                       {"writable_Vws", false},
                                                       {"writable_dVws", false}});
    ModuleBase::GlobalFunc::DONE(this->runtime.log, "cal_large_Cs");

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> tmp;
    if (this->runtime.input.out_unshrinked_v)
    {
        this->Vs_period = RI::RI_Tools::cal_period(Vs_cut_IJR, period);
        Vs_cut_IJR.clear();
        ModuleBase::GlobalFunc::DONE(this->runtime.log, "Vs_period");
        // MPI: {ia0, {ia1, R}} to {ia0, ia1}
        const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, TC>>>> list_As_Vs_atoms
            = RI::Distribute_Equally::distribute_atoms(this->mpi_comm, atoms, period_Vs, 2, false);
        const auto list_A0_pair_R = list_As_Vs_atoms.first;
        const auto list_A1_pair_R = list_As_Vs_atoms.second[0];
        std::set<TA> atoms00;
        std::set<TA> atoms01;
        for (const auto& I: list_A0_pair_R)
        {
            atoms00.insert(I);
        }
        for (const auto& JR: list_A1_pair_R)
        {
            atoms01.insert(JR.first);
        }

        this->Vs_period = RI_2D_Comm::comm_map2_first(this->mpi_comm, this->Vs_period, atoms00, atoms01);
        ModuleBase::GlobalFunc::DONE(this->runtime.log, "Vs_period_comm");

        Exx_LRI<double>* cut_coulomb = exx_cut_coulomb.get();
        this->out_coulomb_k(ucell, this->Vs_period, "coulomb_unshrinked_cut_", cut_coulomb);
        ModuleBase::GlobalFunc::DONE(this->runtime.log, "out_large_Vs");
        this->Vs_period.clear();
        this->Vs_period.swap(tmp);
    }

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs = std::get<0>(Cs_dCs);
    this->Cs_period = RI::RI_Tools::cal_period(Cs, period);
    this->Cs_period = exx_cut_coulomb->exx_lri.post_2D.set_tensors_map2(this->Cs_period);
    if (this->runtime.input.out_librpa_ver == 1)
    {
        this->out_librpa_basis_v1(ucell, exx_cut_coulomb.get(), "aux_basis.txt", "basis_map.txt");
        this->out_Cs_v1(ucell, this->Cs_period, "Cs_");
    }
    else
    {
        this->out_Cs(ucell, this->Cs_period, "Cs_data_");
    }
    ModuleBase::GlobalFunc::DONE(this->runtime.log, "out_large_Cs");
    this->Cs_period.clear();
    this->Cs_period.swap(tmp);
    exx_cut_coulomb.reset();
    RpaLriDetail::trim_malloc_cache();

    ModuleBase::timer::end("RPA_LRI", "cal_large_Cs");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::cal_abfs_overlap(const UnitCell& ucell, const LCAO_Orbitals& orb, const K_Vectors& kv)
{
    ModuleBase::TITLE("DFT_RPA_interface", "cal_abfs_overlap");
    const auto& abfs_s = this->abfs_shrink;

    // <smaller abfs|smaller abfs>
    Matrix_Orbs11 m_abfs_abfs;
    // <smaller abfs|larger abfs>
    Matrix_Orbs11 m_abfs_abf;

    m_abfs_abf.MGT = this->MGT;
    m_abfs_abf.init(abfs_s, this->abfs, ucell, orb, this->info.kmesh_times);
    m_abfs_abf.init_radial_table();

    m_abfs_abfs.MGT = this->MGT;
    m_abfs_abfs.init(abfs_s, abfs_s, ucell, orb, this->info.kmesh_times);
    m_abfs_abfs.init_radial_table();
    // get Rlist
    const std::array<Tcell, Ndim> period = RI_Util::get_Born_vonKarmen_period(kv);
    const auto R_period = RI_Util::get_Born_von_Karmen_cells(period);
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> overlap_abfs_abfs;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> overlap_abfs_abf;

    // index of smaller abfs
    const ModuleBase::Element_Basis_Index::Range range_abfs_s
        = ModuleBase::Element_Basis_Index::construct_range(abfs_s);
    const ModuleBase::Element_Basis_Index::IndexLNM index_abfs_s
        = ModuleBase::Element_Basis_Index::construct_index(range_abfs_s);
    // index of larger abfs
    const ModuleBase::Element_Basis_Index::Range range_abfs
        = ModuleBase::Element_Basis_Index::construct_range(this->abfs);
    const ModuleBase::Element_Basis_Index::IndexLNM index_abfs
        = ModuleBase::Element_Basis_Index::construct_index(range_abfs);

    auto orb_cutoff_ = orb.cutoffs();
    const std::array<Tcell, Ndim> period_Vs = LRI_CV_Tools::cal_latvec_range<Tcell>(2, ucell, orb_cutoff_);
    std::vector<TA> atoms(ucell.nat);
    for (int iat = 0; iat < ucell.nat; ++iat)
        atoms[iat] = iat;
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>> list_As_Vs
        = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Vs, 2, false);

// Huanjing Gong debug
// std::stringstream ss;
//  ss << "IJR_" << this->runtime.rank << ".txt";
// std::ofstream ofs;
// ofs.open(ss.str().c_str(), std::ios::out);
// for (size_t iA = 0; iA < list_As_Vs.first.size(); ++iA)
// {
//     const auto& A = list_As_Vs.first[iA];
//     for (const auto& BR: list_As_Vs.second[0])
//     {
//         const auto& B = BR.first;
//         const auto& R = BR.second;
//         ofs << "ABR: " << A << B << "," << R.at(0) << R.at(1) << R.at(2) << std::endl;
//     }
// }
// ofs.close();
#pragma omp parallel
    {
        using LocalMapType = std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>;
        std::map<int, LocalMapType> overlap_abfs_abfs_local;
        std::map<int, LocalMapType> overlap_abfs_abf_local;

#pragma omp for schedule(dynamic) nowait
        for (size_t iA = 0; iA < list_As_Vs.first.size(); ++iA)
        {
            const auto& A = list_As_Vs.first[iA];
            for (const auto& BR: list_As_Vs.second[0])
            {
                const auto& B = BR.first;
                const auto& R = BR.second;

                const size_t TA = ucell.iat2it[A];
                const size_t IA = ucell.iat2ia[A];
                const auto& tauA = ucell.atoms[TA].tau[IA];
                const size_t TB = ucell.iat2it[B];
                const size_t IB = ucell.iat2ia[B];
                const auto& tauB = ucell.atoms[TB].tau[IB];

                const ModuleBase::Vector3<double> tauB_shift = tauB + (RI_Util::array3_to_Vector3(R) * ucell.latvec);
                const ModuleBase::Vector3<double> tau_delta = tauB_shift - tauA;
                static const ModuleBase::Vector3<double> tau0(0.0, 0.0, 0.0);

                auto& local_abfs_abfs = overlap_abfs_abfs_local[A];
                local_abfs_abfs[{B, R}]
                    = m_abfs_abfs.template cal_overlap_matrix<double>(TA,
                                                                      TB,
                                                                      tau0,
                                                                      tau_delta,
                                                                      index_abfs_s,
                                                                      index_abfs_s,
                                                                      Matrix_Orbs11::Matrix_Order::AB);

                auto& local_abfs_abf = overlap_abfs_abf_local[A];
                local_abfs_abf[{B, R}]
                    = m_abfs_abf.template cal_overlap_matrix<double>(TA,
                                                                     TB,
                                                                     tau0,
                                                                     tau_delta,
                                                                     index_abfs_s,
                                                                     index_abfs,
                                                                     Matrix_Orbs11::Matrix_Order::AB);
            }
        }

#pragma omp critical(RPA_LRI_merge)
        {
            for (auto& aPair: overlap_abfs_abfs_local)
            {
                auto& aKey = aPair.first;
                auto& aSubMap = aPair.second;
                for (auto& subPair: aSubMap)
                {
                    auto& key = subPair.first;
                    auto& value = subPair.second;
                    overlap_abfs_abfs[aKey][key] = std::move(value);
                }
            }
            for (auto& aPair: overlap_abfs_abf_local)
            {
                auto& aKey = aPair.first;
                auto& aSubMap = aPair.second;
                for (auto& subPair: aSubMap)
                {
                    auto& key = subPair.first;
                    auto& value = subPair.second;
                    overlap_abfs_abf[aKey][key] = std::move(value);
                }
            }
        }
    }
    // MPI: {ia0, {ia1, R}} to {ia0, ia1}
    const std::array<Tcell, Ndim> period_Vs_IJ = LRI_CV_Tools::cal_latvec_range<Tcell>(2, ucell, orb_cutoff_);
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, TC>>>> list_As_Vs_atoms
        = RI::Distribute_Equally::distribute_atoms(this->mpi_comm, atoms, period_Vs, 2, false);
    const auto list_A0_pair_R = list_As_Vs_atoms.first;
    const auto list_A1_pair_R = list_As_Vs_atoms.second[0];
    std::set<TA> atoms00;
    std::set<TA> atoms01;
    for (const auto& I: list_A0_pair_R)
    {
        atoms00.insert(I);
    }
    for (const auto& JR: list_A1_pair_R)
    {
        atoms01.insert(JR.first);
    }
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> overlap_abfs_abfs_IJ
        = RI_2D_Comm::comm_map2_first(mpi_comm, overlap_abfs_abfs, atoms00, atoms01);
    overlap_abfs_abfs.clear();
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> overlap_abfs_abf_IJ
        = RI_2D_Comm::comm_map2_first(mpi_comm, overlap_abfs_abf, atoms00, atoms01);
    overlap_abfs_abf.clear();

    if (this->runtime.input.out_librpa_ver == 1)
    {
        out_abfs_overlap_v1(ucell, overlap_abfs_abfs_IJ, overlap_abfs_abf_IJ, "sinvS_", index_abfs_s, index_abfs);
    }
    else
    {
        out_abfs_overlap(ucell, overlap_abfs_abfs_IJ, overlap_abfs_abf_IJ, "shrink_sinvS_", index_abfs_s, index_abfs);
    }
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::inverse_olp(const UnitCell& ucell,
                                    std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>>& overlap_abfs_abfs,
                                    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s)
{
    ModuleBase::TITLE("RPA_LRI", "inverse_olp");
    ModuleBase::timer::start("RPA_LRI", "inverse_olp");
    const int nks_tot = this->runtime.input.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    size_t all_mu_s = 0;
    std::vector<int> mu_s_shift(ucell.nat);
    for (int I = 0; I != ucell.nat; I++)
    {
        mu_s_shift[I] = all_mu_s;
        all_mu_s += index_abfs_s[ucell.iat2it[I]].count_size;
    }
    RI::Tensor<std::complex<double>> olp_all = RI::Tensor<std::complex<double>>({all_mu_s, all_mu_s});
    for (int ik = 0; ik < nks_tot; ik++)
    {
        for (auto& Ip: overlap_abfs_abfs)
        {
            auto I = Ip.first;
            size_t mu_s_I = index_abfs_s[ucell.iat2it[I]].count_size;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;
                auto q = JPp.first.second;
                if (q != RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]))
                    continue;
                // std::cout << "IJ: " << I << "," << J << std::endl;
                auto mu_s_J = index_abfs_s[ucell.iat2it[J]].count_size;
                for (int ir = 0; ir < mu_s_I; ir++)
                {
                    for (int ic = 0; ic < mu_s_J; ic++)
                    {
                        olp_all(mu_s_shift[I] + ir, mu_s_shift[J] + ic) = JPp.second(ir, ic);
                    }
                }
            }
        }
        // for multi-mpi
        // for (int ir = 0; ir < all_mu_s; ir++)
        // {
        //     for (int ic = 0; ic < all_mu_s; ic++)
        //     {
        //         Parallel_Reduce::reduce_all<std::complex<double>>(olp_all(ir, ic));
        //     }
        // }

        // check Hermitian
        for (int ir = 0; ir < all_mu_s; ir++)
        {
            for (int ic = ir; ic < all_mu_s; ic++)
            {
                auto delta = std::abs(olp_all(ir, ic) - std::conj(olp_all(ic, ir)));
                if (delta > 1e-10)
                {
                    std::cout << "Warning: olp_all is not Hermitian!" << std::endl;
                    std::cout << "ik,ir,ic: " << ik << "," << ir << "," << ic << std::endl;
                    std::cout << "delta(ir, ic): " << delta << std::endl;
                }
            }
        }
        // out_pure_ri_tensor("olp_all.txt", olp_all, 0.);
        auto olp_inv = LRI_CV_Tools::cal_I(olp_all,
                                           Inverse_Matrix<std::complex<double>>::Method::syev,
                                           this->info.shrink_LU_inv_thr);
        for (int ir = 0; ir < all_mu_s; ir++)
        {
            for (int ic = ir; ic < all_mu_s; ic++)
            {
                olp_inv(ic, ir) = std::conj(olp_inv(ir, ic));
            }
        }
        // out_pure_ri_tensor("olp_inv.txt", olp_inv, 0.);
        for (auto& Ip: overlap_abfs_abfs)
        {
            auto I = Ip.first;
            size_t mu_s_I = index_abfs_s[ucell.iat2it[I]].count_size;
            for (auto& JPp: Ip.second)
            {
                auto q = JPp.first.second;
                if (q != RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]))
                    continue;
                auto J = JPp.first.first;
                auto mu_s_J = index_abfs_s[ucell.iat2it[J]].count_size;

                for (int ir = 0; ir < mu_s_I; ir++)
                {
                    for (int ic = 0; ic < mu_s_J; ic++)
                        JPp.second(ir, ic) = olp_inv(mu_s_shift[I] + ir, mu_s_shift[J] + ic);
                }
            }
        }
    }
    ModuleBase::timer::end("RPA_LRI", "inverse_olp");
}

template void RPA_LRI<double, double>::cal_large_Cs(const UnitCell& ucell,
                                                    const LCAO_Orbitals& orb,
                                                    const K_Vectors& kv);
template void RPA_LRI<std::complex<double>, double>::cal_large_Cs(const UnitCell& ucell,
                                                                  const LCAO_Orbitals& orb,
                                                                  const K_Vectors& kv);
template void RPA_LRI<double, double>::cal_abfs_overlap(const UnitCell& ucell,
                                                        const LCAO_Orbitals& orb,
                                                        const K_Vectors& kv);
template void RPA_LRI<std::complex<double>, double>::cal_abfs_overlap(const UnitCell& ucell,
                                                                      const LCAO_Orbitals& orb,
                                                                      const K_Vectors& kv);
template void RPA_LRI<double, double>::inverse_olp(
    const UnitCell& ucell,
    std::map<int, std::map<std::pair<int, std::array<double, 3>>, RI::Tensor<std::complex<double>>>>& overlap_abfs_abfs,
    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s);
template void RPA_LRI<std::complex<double>, double>::inverse_olp(
    const UnitCell& ucell,
    std::map<int, std::map<std::pair<int, std::array<double, 3>>, RI::Tensor<std::complex<double>>>>& overlap_abfs_abfs,
    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s);
