//=======================
// AUTHOR : Rong Shi
// DATE :   2022-12-09
//=======================

#ifndef RPA_LRI_HPP
#define RPA_LRI_HPP
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>
#include "source_lcao/module_ri/module_exx_symmetry/symm_rotation.h"

#include "rpa_lri.h"
#include "source_basis/module_ao/elem_basis_idx_orb.h"
#include "source_base/global_function.h"
#include "source_estate/elecstate_lcao.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_lr/utils/spectrum_mo.hpp"

#if defined(__GLIBC__)
#include <malloc.h>
#endif

namespace RpaLriDetail
{
inline void trim_malloc_cache()
{
#if defined(__GLIBC__)
    malloc_trim(0);
#endif
}
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::postSCF(const UnitCell& ucell,
                                const MPI_Comm& mpi_comm_in,
                                const elecstate::DensityMatrix<T, Tdata>& dm,
                                const elecstate::ElecState* pelec,
                                const K_Vectors& kv,
                                const LCAO_Orbitals& orb,
                                const Parallel_Orbitals& parav,
                                const psi::Psi<T>& psi)
{
    ModuleBase::TITLE("RPA_LRI", "postSCF");
    ModuleBase::timer::start("RPA_LRI", "postSCF");
    ModuleBase::GlobalFunc::MAKE_DIR(outdir);

    this->cal_postSCF_exx(dm, mpi_comm_in, ucell, kv, orb);
    this->init(mpi_comm_in, kv, orb.cutoffs());
    this->out_bands(pelec);
    this->out_eigen_vector(parav, psi);
    this->out_struc(ucell);

    std::cout << "rpa_pca_threshold: " << this->info.pca_threshold << std::endl;
    std::cout << "rpa_ccp_rmesh_times: " << this->info.ccp_rmesh_times << std::endl;
    std::cout << "rpa_lcao_exx(Ha): " << std::fixed << std::setprecision(15) << exx_cut_coulomb->Eexx / 2.0 << std::endl;

    std::cout << "etxc(Ha): " << std::fixed << std::setprecision(15) << pelec->f_en.etxc / 2.0 << std::endl;
    std::cout << "etot(Ha): " << std::fixed << std::setprecision(15) << pelec->f_en.etot / 2.0 << std::endl;
    std::cout << "Etot_without_rpa(Ha): " << std::fixed << std::setprecision(15)
              << (pelec->f_en.etot - pelec->f_en.etxc + exx_cut_coulomb->Eexx) / 2.0 << std::endl;
    delete exx_cut_coulomb;
    exx_cut_coulomb = nullptr;
    RpaLriDetail::trim_malloc_cache();

    if (this->info.shrink_abfs_pca_thr >= 0.0)
    {
        cal_large_Cs(ucell, orb, kv);
        cal_abfs_overlap(ucell, orb, kv);
        RpaLriDetail::trim_malloc_cache();
    }
    this->output_ewald_coulomb(ucell, kv, orb);

    ModuleBase::timer::end("RPA_LRI", "postSCF");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::init(const MPI_Comm& mpi_comm_in, const K_Vectors& kv_in, const std::vector<double>& orb_cutoff)
{
    ModuleBase::TITLE("RPA_LRI", "init");
    ModuleBase::timer::start("RPA_LRI", "init");
    this->mpi_comm = mpi_comm_in;
    this->orb_cutoff_ = orb_cutoff;
    this->lcaos = exx_cut_coulomb->lcaos;
    this->p_kv = &kv_in;
    this->MGT = exx_cut_coulomb->MGT;

    if (this->info.shrink_abfs_pca_thr >= 0.0)
    {
        this->abfs_shrink = exx_cut_coulomb->abfs;
    }
    else
    {
        this->abfs = exx_cut_coulomb->abfs;
    }
    //	this->cv = std::move(exx_lri_rpa.cv);
    //    exx_lri_rpa.cv = exx_lri_rpa.cv;
    ModuleBase::timer::end("RPA_LRI", "init");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::cal_postSCF_exx(const elecstate::DensityMatrix<T, Tdata>& dm,
                                        const MPI_Comm& mpi_comm_in,
                                        const UnitCell& ucell,
                                        const K_Vectors& kv,
                                        const LCAO_Orbitals& orb)
{
    ModuleBase::TITLE("RPA_LRI", "cal_postSCF_exx");
    ModuleBase::timer::start("RPA_LRI", "cal_postSCF_exx");

    this->mpi_comm = mpi_comm_in;
    this->p_kv = &kv;
    this->orb_cutoff_ = orb.cutoffs();

    Mix_DMk_2D<T> mix_DMk_2D;
    bool exx_spacegroup_symmetry = (PARAM.inp.nspin < 4 && ModuleSymmetry::Symmetry::symm_flag == 1);
    if (exx_spacegroup_symmetry)
        {mix_DMk_2D.set_nks(kv.get_nkstot_full() * (PARAM.inp.nspin == 2 ? 2 : 1));}
    else
        {mix_DMk_2D.set_nks(kv.get_nks());}
        
    mix_DMk_2D.set_mixing_plain(1.0);

    // --------------------------------------------------------------------------------------
    // NOTE: ABFs are constructed here BEFORE symmetry processing to calculate the correct
    // abfs_Lmax for symrot.set_abfs_Lmax(). Previously, abfs_Lmax was obtained either from
    // exx_cut_coulomb->abfs_Lmax() (which was nullptr at that point) or from this->info.abfs_Lmax
    // (which defaults to 0). This caused abfs_Lmax to degrade to 0 in RPA + symmetry path,
    // leading to incorrect rotation matrices for ABFs with higher angular momentum.
    // --------------------------------------------------------------------------------------
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs_for_lmax;
    if (this->info.shrink_abfs_pca_thr >= 0.0)
    {
        this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs(orb, this->info.kmesh_times);
        abfs_for_lmax = Exx_Abfs::Construct_Orbs::abfs_same_atom(ucell, orb, this->lcaos, this->info.kmesh_times, this->info.shrink_abfs_pca_thr);
        if (!this->info.files_shrink_abfs.empty())
        {
            abfs_for_lmax = Exx_Abfs::IO::construct_abfs(abfs_for_lmax, orb, this->info.files_shrink_abfs, this->info.kmesh_times);
        }
    }
    else
    {
        this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs(orb, this->info.kmesh_times);
        abfs_for_lmax = Exx_Abfs::Construct_Orbs::abfs_same_atom(ucell, orb, this->lcaos, this->info.kmesh_times, this->info.pca_threshold);
        if (!this->info.files_abfs.empty())
        {
            abfs_for_lmax = Exx_Abfs::IO::construct_abfs(abfs_for_lmax, orb, this->info.files_abfs, this->info.kmesh_times);
        }
    }

    ModuleSymmetry::Symmetry_rotation symrot;
    if (exx_spacegroup_symmetry)
    {
        const std::array<Tcell, Ndim> period = RI_Util::get_Born_vonKarmen_period(kv);
        const auto& Rs = RI_Util::get_Born_von_Karmen_cells(period);
        symrot.find_irreducible_sector(ucell.symm, ucell.atoms, ucell.st, Rs, period, ucell.lat);
        // set Lmax of the rotation matrices to max(l_ao, l_abf), to support rotation under ABF
        // NOTE: Using Exx_Abfs::Construct_Orbs::get_Lmax() to compute Lmax from the actual ABFs
        // instead of relying on exx_cut_coulomb->abfs_Lmax() (not yet initialized) or
        // this->info.abfs_Lmax (defaults to 0). This ensures correct Lmax for symmetry rotation.
        symrot.set_abfs_Lmax(Exx_Abfs::Construct_Orbs::get_Lmax(abfs_for_lmax));
        symrot.cal_Ms(kv, ucell, *dm.get_paraV_pointer());
        // output Ts (symrot_R.txt) and Ms (symrot_k.txt)
        ModuleSymmetry::print_symrot_info_R(symrot, ucell.symm, ucell.lmax, Rs);
        ModuleSymmetry::print_symrot_info_k(symrot, kv, ucell);
        mix_DMk_2D.mix(symrot.restore_dm(kv, dm.get_DMK_vector(), *dm.get_paraV_pointer()), true);
    }
    else { mix_DMk_2D.mix(dm.get_DMK_vector(), true); }
    
    const std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>>
        Ds = RI_2D_Comm::split_m2D_ktoR<Tdata>(
            ucell,
            kv,
            mix_DMk_2D.get_DMk_out(),
            *dm.get_paraV_pointer(),
            PARAM.inp.nspin,
            exx_spacegroup_symmetry);
    
    // reserve exx_ccp_rmesh_times to calculate full Coulomb
    // Note: ccp_type=Hf and hybrid_alpha=1 were previously set on the global Exx_Info
    // and sync_from_global() was called, but this->info (value copy) already has the correct
    // coulomb_param from construction time, so those writes are redundant and removed.
    this->ccp_rmesh_times_ewald = this->info.ccp_rmesh_times;
    // Using rpa_ccp_rmesh_times to calculate cut Coulomb this->Vs_period
    Exx_Info_RI local_info = this->info;
    local_info.ccp_rmesh_times = PARAM.inp.rpa_ccp_rmesh_times;
    if (!exx_cut_coulomb)
        exx_cut_coulomb = new Exx_LRI<double>(local_info);

    if (this->info.shrink_abfs_pca_thr >= 0.0)
    {
        // NOTE: Reuse abfs_for_lmax constructed earlier to avoid redundant ABFs construction.
        // This ensures consistency between the ABFs used for Lmax calculation and the ABFs
        // used for actual EXX computation.
        this->abfs_shrink = abfs_for_lmax;
        Exx_Abfs::Construct_Orbs::print_orbs_size(ucell, abfs_shrink, GlobalV::ofs_running);
        exx_cut_coulomb->init_spencer(mpi_comm_in, ucell, kv, orb, abfs_shrink);
    }
    else
        // NOTE: Reuse abfs_for_lmax constructed earlier to avoid redundant ABFs construction.
        exx_cut_coulomb->init_spencer(mpi_comm_in, ucell, kv, orb, abfs_for_lmax);
    // cal C and V for exx
    this->output_cut_coulomb_cs(ucell, exx_cut_coulomb);
    // cal CVCD
    if (exx_spacegroup_symmetry && PARAM.inp.exx_symmetry_realspace)
    {
        exx_cut_coulomb->cal_exx_elec(Ds, ucell, *dm.get_paraV_pointer(), &symrot);
    }
    else
    {
        exx_cut_coulomb->cal_exx_elec(Ds, ucell, *dm.get_paraV_pointer());
    }
    // cout<<"postSCF_Eexx: "<<exx_lri_rpa.Eexx<<endl;
    ModuleBase::timer::end("RPA_LRI", "cal_postSCF_exx");
}

// if use shrink, output Coulomb and Cs_data in small abfs
// otherwise, output in normal abfs
template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::output_cut_coulomb_cs(const UnitCell& ucell, Exx_LRI<double>* exx_lri_rpa)
{
    ModuleBase::TITLE("RPA_LRI", "output_cut_coulomb_cs");
    ModuleBase::timer::start("RPA_LRI", "output_cut_coulomb_cs");

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_cut_IJR;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Cs;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> tmp;
    std::cout << "Use rpa_ccp_rmesh_times=" << this->info.ccp_rmesh_times << " to calculate cut Coulomb" << std::endl;
    // Shrink_ABFS_ORBITAL cannot exceed this angular momentum of MGT
    exx_lri_rpa->cal_cut_coulomb_cs(Vs_cut_IJR, Cs, ucell, PARAM.inp.out_ri_cv);
    // MPI: {ia0, {ia1, R}} to {ia0, ia1}
    std::vector<TA> atoms(ucell.nat);
    for (int iat = 0; iat < ucell.nat; ++iat)
        atoms[iat] = iat;
    const std::array<Tcell, Ndim> period_Vs
        = LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->info.ccp_rmesh_times, ucell, this->orb_cutoff_);
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
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_cut_IJ
        = RI_2D_Comm::comm_map2_first(this->mpi_comm, Vs_cut_IJR, atoms00, atoms01);
    Vs_cut_IJR.clear();
    const std::array<Tcell, Ndim> period = {p_kv->nmp[0], p_kv->nmp[1], p_kv->nmp[2]};
    this->Vs_period = RI::RI_Tools::cal_period(Vs_cut_IJ, period);
    this->out_coulomb_k(ucell, this->Vs_period, "coulomb_cut_", exx_lri_rpa);
    Vs_period.clear();
    Vs_period.swap(tmp);

    this->Cs_period = RI::RI_Tools::cal_period(Cs, period);
    this->Cs_period = exx_lri_rpa->exx_lri.post_2D.set_tensors_map2(this->Cs_period);

    if (this->info.shrink_abfs_pca_thr >= 0.0)
        this->out_Cs(ucell, this->Cs_period, "Cs_shrinked_data_");
    else
        this->out_Cs(ucell, this->Cs_period, "Cs_data_");
    Cs_period.clear();
    Cs_period.swap(tmp);

    ModuleBase::timer::end("RPA_LRI", "output_cut_coulomb_cs");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::output_ewald_coulomb(const UnitCell& ucell, const K_Vectors& kv, const LCAO_Orbitals& orb)
{
    ModuleBase::TITLE("RPA_LRI", "output_ewald_coulomb");
    ModuleBase::timer::start("RPA_LRI", "output_ewald_coulomb");

    Exx_Info_RI local_info = this->info;
    local_info.ccp_rmesh_times = this->ccp_rmesh_times_ewald;
    if (!exx_full_coulomb)
        exx_full_coulomb = new Exx_LRI<double>(local_info);

    if (this->info.shrink_abfs_pca_thr >= 0.0)
        exx_full_coulomb->init(mpi_comm, ucell, kv, orb, this->abfs_shrink);
    else
        exx_full_coulomb->init(mpi_comm, ucell, kv, orb, this->abfs);
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_full_IJR;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Cs;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> tmp;
    exx_full_coulomb->cal_ewald_coulomb(Vs_full_IJR, Cs, ucell, PARAM.inp.out_ri_cv);
    // MPI: {ia0, {ia1, R}} to {ia0, ia1}
    std::vector<TA> atoms(ucell.nat);
    for (int iat = 0; iat < ucell.nat; ++iat)
        atoms[iat] = iat;
    const std::array<Tcell, Ndim> period_Vs
        = LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->ccp_rmesh_times_ewald, ucell, this->orb_cutoff_);
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, TC>>>> list_As_Vs_atoms
        = RI::Distribute_Equally::distribute_atoms(mpi_comm, atoms, period_Vs, 2, false);
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
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_full_IJ
        = RI_2D_Comm::comm_map2_first(mpi_comm, Vs_full_IJR, atoms00, atoms01);
    Vs_full_IJR.clear();

    const std::array<Tcell, Ndim> period = {p_kv->nmp[0], p_kv->nmp[1], p_kv->nmp[2]};
    this->Vs_period = RI::RI_Tools::cal_period(Vs_full_IJ, period);
    this->out_coulomb_k(ucell, this->Vs_period, "coulomb_mat_", exx_full_coulomb);
    Vs_period.clear();
    Vs_period.swap(tmp);
    Cs.clear();
    Cs.swap(tmp);

    delete exx_full_coulomb;
    exx_full_coulomb = nullptr;
    RpaLriDetail::trim_malloc_cache();

    ModuleBase::timer::end("RPA_LRI", "output_ewald_coulomb");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::cal_large_Cs(const UnitCell& ucell, const LCAO_Orbitals& orb, const K_Vectors& kv)
{
    ModuleBase::TITLE("RPA_LRI", "cal_large_Cs");
    ModuleBase::timer::start("RPA_LRI", "cal_large_Cs");
    if (!exx_cut_coulomb)
        exx_cut_coulomb = new Exx_LRI<double>(this->info);
    exx_cut_coulomb->init_spencer(this->mpi_comm, ucell, kv, orb);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "exx_cut_coulomb->init");
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
        = LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->info.ccp_rmesh_times, ucell, orb_cutoff_);
    std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>> list_As_Vs
        = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Vs, 2, false);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal_large_Vs start");
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_cut_IJR
        = center2_obj_it->second.cv.cal_Vs(ucell, list_As_Vs.first, list_As_Vs.second[0], {{"writable_Vws", true}});
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal_large_Vs end");

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
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal_large_Cs");

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> tmp;
    if (PARAM.inp.out_unshrinked_v)
    {
        this->Vs_period = RI::RI_Tools::cal_period(Vs_cut_IJR, period);
        Vs_cut_IJR.clear();
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "Vs_period");
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
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "Vs_period_comm");

        this->out_coulomb_k(ucell, this->Vs_period, "coulomb_unshrinked_cut_", exx_cut_coulomb);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "out_large_Vs");
        this->Vs_period.clear();
        this->Vs_period.swap(tmp);
    }

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs = std::get<0>(Cs_dCs);
    this->Cs_period = RI::RI_Tools::cal_period(Cs, period);
    this->Cs_period = exx_cut_coulomb->exx_lri.post_2D.set_tensors_map2(this->Cs_period);
    this->out_Cs(ucell, this->Cs_period, "Cs_data_");
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "out_large_Cs");
    this->Cs_period.clear();
    this->Cs_period.swap(tmp);
    delete exx_cut_coulomb;
    exx_cut_coulomb = nullptr;
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
    const ModuleBase::Element_Basis_Index::Range range_abfs_s = ModuleBase::Element_Basis_Index::construct_range(abfs_s);
    const ModuleBase::Element_Basis_Index::IndexLNM index_abfs_s
        = ModuleBase::Element_Basis_Index::construct_index(range_abfs_s);
    // index of larger abfs
    const ModuleBase::Element_Basis_Index::Range range_abfs = ModuleBase::Element_Basis_Index::construct_range(this->abfs);
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
//  ss << "IJR_" << GlobalV::MY_RANK << ".txt";
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

                const ModuleBase::Vector3<double> tauB_shift
                    = tauB + (RI_Util::array3_to_Vector3(R) * ucell.latvec);
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

    out_abfs_overlap(ucell, overlap_abfs_abfs_IJ, overlap_abfs_abf_IJ, "shrink_sinvS_", index_abfs_s, index_abfs);
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_abfs_overlap(const UnitCell& ucell,
                                         std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abfs,
                                         std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abf,
                                         std::string filename,
                                         const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s,
                                         const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs)
{
    ModuleBase::TITLE("RPA_LRI", "out_abfs_overlap");
    ModuleBase::timer::start("RPA_LRI", "out_abfs_overlap");
    const double threshold = 1e-15;
    const auto format = std::scientific;
    int prec = 15;

    int all_mu_s = 0;
    int all_mu = 0;
    std::vector<int> mu_s_shift(ucell.nat);
    std::vector<int> mu_shift(ucell.nat);
    for (int I = 0; I != ucell.nat; I++)
    {
        mu_s_shift[I] = all_mu_s;
        mu_shift[I] = all_mu;
        all_mu_s += index_abfs_s[ucell.iat2it[I]].count_size;
        all_mu += index_abfs[ucell.iat2it[I]].count_size;
    }
    const int nks_tot = PARAM.inp.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    std::stringstream ss;
    ss << filename << GlobalV::MY_RANK << ".txt";

    std::ofstream ofs;
    ofs.open(outdir + ss.str().c_str(), std::ios::out);

    ofs << nks_tot << std::endl;

    // Fourier of ss(R->k), s(R->k)
    std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>> olp_q_ss;
    std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>> olp_q_s;
    for (int ik = 0; ik != nks_tot; ik++)
    {
        for (auto& Ip: overlap_abfs_abfs)
        {
            auto I = Ip.first;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;
                auto R = JPp.first.second;
                auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
                RI::Tensor<std::complex<double>> tmp_olp_ss
                    = RI::Global_Func::convert<std::complex<double>>(JPp.second);
                RI::Tensor<std::complex<double>> tmp_olp_s
                    = RI::Global_Func::convert<std::complex<double>>(overlap_abfs_abf[I][{J, R}]);
                if (olp_q_ss[I][{J, q}].empty())
                {
                    olp_q_ss[I][{J, q}] = RI::Tensor<std::complex<double>>({tmp_olp_ss.shape[0], tmp_olp_ss.shape[1]});
                    olp_q_s[I][{J, q}] = RI::Tensor<std::complex<double>>({tmp_olp_s.shape[0], tmp_olp_s.shape[1]});
                }
                const double arg = 1 * (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec))
                                   * ModuleBase::TWO_PI; // latvec
                const std::complex<double> kphase = std::complex<double>(cos(arg), sin(arg));

                olp_q_ss[I][{J, q}] = olp_q_ss[I][{J, q}] + tmp_olp_ss * kphase;
                olp_q_s[I][{J, q}] = olp_q_s[I][{J, q}] + tmp_olp_s * kphase;
            }
        }
    }
    // for multi-mpi
    for (int I = 0; I != ucell.nat; I++)
    {
        for (int J = 0; J != ucell.nat; J++)
        {
            for (int ik = 0; ik != nks_tot; ik++)
            {
                auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
                if (olp_q_ss[I][{J, q}].empty())
                {
                    auto mu = index_abfs_s[ucell.iat2it[I]].count_size;
                    auto nu = index_abfs_s[ucell.iat2it[J]].count_size;
                    olp_q_ss[I][{J, q}] = RI::Tensor<std::complex<double>>({mu, nu});
                }
                if (olp_q_s[I][{J, q}].empty())
                {
                    auto mu = index_abfs_s[ucell.iat2it[I]].count_size;
                    auto nu = index_abfs[ucell.iat2it[J]].count_size;
                    olp_q_s[I][{J, q}] = RI::Tensor<std::complex<double>>({mu, nu});
                }
                for (int ir = 0; ir < olp_q_ss[I][{J, q}].shape[0]; ir++)
                {
                    for (int ic = 0; ic < olp_q_ss[I][{J, q}].shape[1]; ic++)
                    {
                        Parallel_Reduce::reduce_all<std::complex<double>>(olp_q_ss[I][{J, q}](ir, ic));
                    }
                    for (int ic = 0; ic < olp_q_s[I][{J, q}].shape[1]; ic++)
                    {
                        Parallel_Reduce::reduce_all<std::complex<double>>(olp_q_s[I][{J, q}](ir, ic));
                    }
                }
            }
        }
    }

    // out_ri_tensor("olp_ss.dat", olp_q_ss, 0.);
    // Inverse of overlap(q)
    inverse_olp(ucell, olp_q_ss, index_abfs_s);
    // out_ri_tensor("olp_ss_inv.dat", olp_q_ss, 0.);
    // out_ri_tensor("olp_s.dat", olp_q_s, 0.);
    for (auto& Ip: overlap_abfs_abf)
    {
        auto I = Ip.first;
        size_t mu_num_s = index_abfs_s[ucell.iat2it[I]].count_size;
        size_t mu_num = index_abfs[ucell.iat2it[I]].count_size;

        for (int ik = 0; ik != nks_tot; ik++)
        {
            std::map<size_t, RI::Tensor<std::complex<double>>> sinvS;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;
                auto R = JPp.first.second;
                if (sinvS[J].empty())
                {
                    sinvS[J] = RI::Tensor<std::complex<double>>(
                        {overlap_abfs_abfs[I][{J, R}].shape[0], overlap_abfs_abf[I][{J, R}].shape[1]});
                }
            }
            for (const auto& pair: sinvS)
            {
                auto J = pair.first;
                auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
                for (int K = 0; K != ucell.nat; K++)
                {
                    sinvS[J] += olp_q_ss.at(I).at({K, q}) * olp_q_s.at(K).at({J, q});
                }
            }
            for (auto& iJU: sinvS)
            {
                auto iJ = iJU.first;
                auto& vq_J = iJU.second;
                size_t nu_num = index_abfs[ucell.iat2it[iJ]].count_size;
                ofs << all_mu_s << "   " << all_mu << "   " << mu_s_shift[I] + 1 << "   " << mu_s_shift[I] + mu_num_s
                    << "  " << mu_shift[iJ] + 1 << "   " << mu_shift[iJ] + nu_num << std::endl;
                ofs << ik + 1 << "  " << p_kv->wk[ik] / 2.0 * PARAM.inp.nspin << std::endl;
                for (int i = 0; i != vq_J.data->size(); i++)
                {
                    // ofs << std::setw(25) << std::fixed << std::setprecision(15) << (*vq_J.data)[i].real()
                    //     << std::setw(25) << std::fixed << std::setprecision(15) << (*vq_J.data)[i].imag() <<
                    //     std::endl;
                    // if (fabs((*vq_J.data)[i].real()) > threshold || fabs((*vq_J.data)[i].imag()) > threshold)
                    ofs << std::showpoint << format << std::setprecision(prec) << (*vq_J.data)[i].real() << " "
                        << std::showpoint << format << std::setprecision(prec) << (*vq_J.data)[i].imag() << "\n";
                    // else
                    //     ofs << std::showpoint << format << std::setprecision(prec) << 0.0 << " " << std::showpoint
                    //         << format << std::setprecision(prec) << 0.0 << "\n";
                }
            }
        }
    }
    ofs.close();
    ModuleBase::timer::end("RPA_LRI", "out_abfs_overlap");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::inverse_olp(const UnitCell& ucell,
                                    std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>>& overlap_abfs_abfs,
                                    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s)
{
    ModuleBase::TITLE("RPA_LRI", "inverse_olp");
    ModuleBase::timer::start("RPA_LRI", "inverse_olp");
    const int nks_tot = PARAM.inp.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
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
        // out_pure_ri_tensor("olp_all.dat", olp_all, 0.);
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
        // out_pure_ri_tensor("olp_inv.dat", olp_inv, 0.);
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

// debug function
// template <typename T, typename Tdata>
// void RPA_LRI<T, Tdata>::out_pure_ri_tensor(const std::string fn,
//                                            RI::Tensor<std::complex<double>>& olp,
//                                            const double threshold)
// {
//     std::ofstream fs;
//     auto format = std::scientific;
//     int prec = 15;
//     fs.open(fn);
//     int nr = olp.shape[0];
//     int nc = olp.shape[1];
//     size_t nnz = nr * nc;
//     fs << "%%MatrixMarket matrix coordinate complex general" << std::endl;
//     fs << "%" << std::endl;

//     fs << nr << " " << nc << " " << nnz << std::endl;

//     for (int j = 0; j < nc; j++)
//     {
//         for (int i = 0; i < nr; i++)
//         {
//             auto v = olp(i, j);
//             if (fabs(v.real()) > threshold || fabs(v.imag()) > threshold)
//                 fs << i + 1 << " " << j + 1 << " " << std::showpoint << format << std::setprecision(prec) << v.real()
//                    << " " << std::showpoint << format << std::setprecision(prec) << v.imag() << "\n";
//         }
//     }

//     fs.close();
// }

// template <typename T, typename Tdata>
// void RPA_LRI<T, Tdata>::out_pure_ri_tensor(const std::string fn, RI::Tensor<double>& olp, const double threshold)
// {
//     std::ofstream fs;
//     auto format = std::scientific;
//     int prec = 15;
//     fs.open(fn);
//     int nr = olp.shape[0];
//     int nc = olp.shape[1];
//     size_t nnz = nr * nc;
//     fs << "%%MatrixMarket matrix coordinate complex general" << std::endl;
//     fs << "%" << std::endl;

//     fs << nr << " " << nc << " " << nnz << std::endl;

//     for (int j = 0; j < nc; j++)
//     {
//         for (int i = 0; i < nr; i++)
//         {
//             auto v = olp(i, j);
//             if (fabs(v) > threshold)
//                 fs << i + 1 << " " << j + 1 << " " << std::showpoint << format << std::setprecision(prec) << v << "\n";
//         }
//     }

//     fs.close();
// }

// template <typename T, typename Tdata>
// void RPA_LRI<T, Tdata>::out_ri_tensor(const std::string fn,
//                                       std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>>& olp,
//                                       const double threshold)
// {
//     std::ofstream fs;
//     auto format = std::scientific;
//     int prec = 15;
//     fs.open(fn);
//     for (auto& IJq: olp)
//     {
//         int I = IJq.first;
//         for (auto& Jq: IJq.second)
//         {
//             int J = Jq.first.first;
//             auto q = Jq.first.second;
//             auto mat = Jq.second;
//             int nr = mat.shape[0];
//             int nc = mat.shape[1];
//             size_t nnz = nr * nc;
//             fs << "%%MatrixMarket matrix coordinate complex general" << std::endl;
//             fs << I << " " << J << " " << q.at(0) << " " << q.at(1) << " " << q.at(2) << std::endl;
//             fs << "%" << std::endl;

//             fs << nr << " " << nc << " " << nnz << std::endl;

//             for (int j = 0; j < nc; j++)
//             {
//                 for (int i = 0; i < nr; i++)
//                 {
//                     auto v = mat(i, j);
//                     if (fabs(v.real()) > threshold || fabs(v.imag()) > threshold)
//                         fs << i + 1 << " " << j + 1 << " " << std::showpoint << format << std::setprecision(prec)
//                            << v.real() << " " << std::showpoint << format << std::setprecision(prec) << v.imag()
//                            << "\n";
//                 }
//             }
//         }
//     }

//     fs.close();
// }

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_eigen_vector(const Parallel_Orbitals& parav, const psi::Psi<T>& psi)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_eigen_vector");

    const int nks_tot = PARAM.inp.nspin == 2 ? p_kv->get_nks() / 2 : p_kv->get_nks();
    const int npsin_tmp = PARAM.inp.nspin == 2 ? 2 : 1;
    const int nbands = parav.get_wfc_global_nbands();
    const int nlocal = parav.get_wfc_global_nbasis();
    const std::size_t values_per_iw = static_cast<std::size_t>(nbands) * npsin_tmp;
    const std::complex<double> zero(0.0, 0.0);

#ifdef __MPI
    const MPI_Comm mpi_comm = parav.comm();
    const int mpi_rank = parav.get_coord_row() * parav.get_dim1() + parav.get_coord_col();
    const int mpi_size = parav.get_dim0() * parav.get_dim1();
    if (mpi_comm == MPI_COMM_NULL)
    {
        throw std::runtime_error("RPA eigenvector MPI-IO received MPI_COMM_NULL.");
    }
    const auto check_mpi = [mpi_comm](const int local_error, const std::string& context) {
        const int local_failed = local_error == MPI_SUCCESS ? 0 : 1;
        int any_failed = 0;
        if (MPI_Allreduce(&local_failed, &any_failed, 1, MPI_INT, MPI_MAX, mpi_comm) != MPI_SUCCESS
            || any_failed != 0)
        {
            throw std::runtime_error(context);
        }
    };
    int communicator_relation = MPI_UNEQUAL;
    check_mpi(MPI_Comm_compare(mpi_comm, this->mpi_comm, &communicator_relation),
              "Failed to compare RPA eigenvector MPI communicators.");
    if (communicator_relation != MPI_IDENT && communicator_relation != MPI_CONGRUENT)
    {
        throw std::runtime_error("RPA and wavefunction MPI communicators are inconsistent.");
    }
    if (mpi_size <= 0 || mpi_rank < 0 || mpi_rank >= mpi_size)
    {
        throw std::runtime_error("Invalid RPA eigenvector BLACS process grid.");
    }

    // Assign a contiguous orbital interval to each rank. Rank order then also
    // gives the order of the corresponding text blocks in the output file.
    std::vector<int> recv_counts(mpi_size, nlocal / mpi_size);
    for (int ip = 0; ip < nlocal % mpi_size; ++ip)
    {
        ++recv_counts[ip];
    }
    const int local_nw = recv_counts[mpi_rank];
#else
    const int mpi_rank = 0;
    const int local_nw = nlocal;
#endif

    for (int ik = 0; ik < nks_tot; ik++)
    {
        std::stringstream ss;
        ss << "KS_eigenvector_" << ik << ".dat";

        const std::size_t local_size = static_cast<std::size_t>(local_nw) * values_per_iw;
        std::vector<std::complex<double>> local_wfc(local_size, zero);
        std::vector<std::complex<double>> wfc_iks(nlocal, zero);
#ifdef __MPI
        std::vector<std::complex<double>> reduced_wfc(local_nw);
        std::complex<double> dummy;
#endif

        for (int is = 0; is < npsin_tmp; is++)
        {
            for (int ib_global = 0; ib_global < nbands; ++ib_global)
            {
                std::fill(wfc_iks.begin(), wfc_iks.end(), zero);

                const int ib_local = parav.global2local_col(ib_global);
                if (ib_local >= 0)
                {
                    for (int ir = 0; ir < psi.get_nbasis(); ir++)
                    {
                        wfc_iks[parav.local2global_row(ir)] = psi(ik + nks_tot * is, ib_local, ir);
                    }
                }

#ifdef __MPI
                MPI_Reduce_scatter(wfc_iks.data(),
                                   local_nw > 0 ? reduced_wfc.data() : &dummy,
                                   recv_counts.data(),
                                   MPI_DOUBLE_COMPLEX,
                                   MPI_SUM,
                                   mpi_comm);
#else
                const std::vector<std::complex<double>>& reduced_wfc = wfc_iks;
#endif
                for (int iw = 0; iw < local_nw; ++iw)
                {
                    const std::size_t index
                        = static_cast<std::size_t>(iw) * values_per_iw + ib_global * npsin_tmp + is;
                    local_wfc[index] = reduced_wfc[iw];
                }
            } // ib
        } // is

        std::ostringstream output;
        output << std::fixed << std::setprecision(15);
        if (mpi_rank == 0)
        {
            output << ik + 1 << '\n';
        }
        for (int iw = 0; iw < local_nw; ++iw)
        {
            for (int ib = 0; ib < nbands; ++ib)
            {
                for (int is = 0; is < npsin_tmp; is++)
                {
                    const std::complex<double>& value
                        = local_wfc[static_cast<std::size_t>(iw) * values_per_iw + ib * npsin_tmp + is];
                    output << std::setw(30) << value.real() << std::setw(30) << value.imag() << '\n';
                }
            }
        }

        const std::string output_buffer = output.str();
        const std::string filename = outdir + ss.str();
#ifdef __MPI
        const unsigned long long local_bytes
            = static_cast<unsigned long long>(output_buffer.size());
        const bool local_size_valid
            = static_cast<std::size_t>(local_bytes) == output_buffer.size();
        if (!local_size_valid)
        {
            throw std::runtime_error("RPA eigenvector local output exceeds uint64_t range.");
        }
        std::vector<unsigned long long> byte_counts(mpi_size, 0);
        check_mpi(MPI_Allgather(&local_bytes,
                                1,
                                MPI_UNSIGNED_LONG_LONG,
                                byte_counts.data(),
                                1,
                                MPI_UNSIGNED_LONG_LONG,
                                mpi_comm),
                  "Failed to collect RPA eigenvector output sizes.");

        unsigned long long file_offset_bytes = 0;
        unsigned long long total_file_bytes = 0;
        unsigned long long max_local_bytes = 0;
        const unsigned long long max_mpi_offset
            = static_cast<unsigned long long>(std::numeric_limits<MPI_Offset>::max());
        bool file_layout_valid = true;
        for (int ip = 0; ip < mpi_size; ++ip)
        {
            if (ip == mpi_rank)
            {
                file_offset_bytes = total_file_bytes;
            }
            if (byte_counts[ip] > max_mpi_offset - total_file_bytes)
            {
                file_layout_valid = false;
                break;
            }
            total_file_bytes += byte_counts[ip];
            max_local_bytes = std::max(max_local_bytes, byte_counts[ip]);
        }
        if (!file_layout_valid)
        {
            throw std::runtime_error("RPA eigenvector file layout exceeds MPI_Offset range.");
        }

        const std::string temporary_filename = filename + ".tmp";
        MPI_File file = MPI_FILE_NULL;
        bool file_opened = false;
        const char dummy_buffer = '\0';
        try
        {
            int mpi_error = MPI_File_open(mpi_comm,
                                          temporary_filename.c_str(),
                                          MPI_MODE_CREATE | MPI_MODE_WRONLY,
                                          MPI_INFO_NULL,
                                          &file);
            check_mpi(mpi_error, "Failed to open " + temporary_filename + ".");
            file_opened = true;
            check_mpi(MPI_File_set_errhandler(file, MPI_ERRORS_RETURN),
                      "Failed to set the RPA eigenvector MPI file error handler.");
            check_mpi(MPI_File_set_size(file, static_cast<MPI_Offset>(total_file_bytes)),
                      "Failed to set the RPA eigenvector file size.");

            const unsigned long long max_write_count
                = static_cast<unsigned long long>(std::numeric_limits<int>::max());
            for (unsigned long long buffer_offset = 0;
                 buffer_offset < max_local_bytes;
                 buffer_offset += max_write_count)
            {
                const unsigned long long remaining
                    = local_bytes > buffer_offset ? local_bytes - buffer_offset : 0;
                const int write_count
                    = static_cast<int>(std::min(remaining, max_write_count));
                const char* write_buffer
                    = write_count > 0
                          ? output_buffer.data() + static_cast<std::size_t>(buffer_offset)
                          : &dummy_buffer;
                MPI_Status status;
                const unsigned long long write_offset
                    = write_count > 0 ? file_offset_bytes + buffer_offset : file_offset_bytes;
                mpi_error = MPI_File_write_at_all(
                    file,
                    static_cast<MPI_Offset>(write_offset),
                    write_buffer,
                    write_count,
                    MPI_CHAR,
                    &status);
                check_mpi(mpi_error, "Failed to write " + temporary_filename + ".");
            }

            check_mpi(MPI_File_sync(file), "Failed to sync " + temporary_filename + ".");
            mpi_error = MPI_File_close(&file);
            file_opened = false;
            check_mpi(mpi_error, "Failed to close " + temporary_filename + ".");

            const int rename_error
                = mpi_rank == 0
                          && std::rename(temporary_filename.c_str(), filename.c_str()) != 0
                      ? MPI_ERR_IO
                      : MPI_SUCCESS;
            check_mpi(rename_error, "Failed to publish " + filename + ".");
        }
        catch (...)
        {
            if (file_opened)
            {
                MPI_File_close(&file);
            }
            MPI_Barrier(mpi_comm);
            if (mpi_rank == 0)
            {
                std::remove(temporary_filename.c_str());
            }
            MPI_Barrier(mpi_comm);
            throw;
        }
#else
        std::ofstream ofs(filename.c_str(), std::ios::out);
        ofs << output_buffer;
        ofs.close();
#endif
    } // ik
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_struc(const UnitCell& ucell)
{
    if (GlobalV::MY_RANK != 0)
    {
        return;
    }
    ModuleBase::TITLE("DFT_RPA_interface", "out_struc");
    double TWOPI_Bohr2A = ModuleBase::TWO_PI * ModuleBase::BOHR_TO_A;
    const int nks_tot = PARAM.inp.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    ModuleBase::Matrix3 lat = ucell.latvec / ModuleBase::BOHR_TO_A;
    ModuleBase::Matrix3 G_RPA = ucell.G * TWOPI_Bohr2A;
    std::stringstream ss;
    ss << "stru_out";
    std::ofstream ofs;
    ofs.open(outdir + ss.str().c_str(), std::ios::out);
    ofs << std::fixed << std::setprecision(9);
    ofs << lat.e11 << std::setw(15) << lat.e12 << std::setw(15) << lat.e13 << std::endl;
    ofs << lat.e21 << std::setw(15) << lat.e22 << std::setw(15) << lat.e23 << std::endl;
    ofs << lat.e31 << std::setw(15) << lat.e32 << std::setw(15) << lat.e33 << std::endl;

    ofs << G_RPA.e11 << std::setw(15) << G_RPA.e12 << std::setw(15) << G_RPA.e13 << std::endl;
    ofs << G_RPA.e21 << std::setw(15) << G_RPA.e22 << std::setw(15) << G_RPA.e23 << std::endl;
    ofs << G_RPA.e31 << std::setw(15) << G_RPA.e32 << std::setw(15) << G_RPA.e33 << std::endl;

    ofs << ucell.nat << std::endl;
    std::string& Coordinate = ucell.Coordinate;
    bool direct = (Coordinate == "Direct");
    // Only consider Direct or Cartesian
    for (int it = 0; it < ucell.ntype; it++)
    {
        Atom* atom = &ucell.atoms[it];
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const double& x = direct ? ucell.atoms[it].tau[ia].x * ucell.lat0
                                     : ucell.atoms[it].tau[ia].x;
            const double& y = direct ? ucell.atoms[it].tau[ia].y * ucell.lat0
                                     : ucell.atoms[it].tau[ia].y;
            const double& z = direct ? ucell.atoms[it].tau[ia].z * ucell.lat0
                                     : ucell.atoms[it].tau[ia].z;
            ofs << std::setw(15) << x << std::setw(15) << y << std::setw(15) << z
                << std::setw(15) << it + 1 << std::endl;
        }
    }

    ofs << p_kv->nmp[0] << std::setw(6) << p_kv->nmp[1] << std::setw(6) << p_kv->nmp[2] << std::setw(6) << std::endl;

    for (int ik = 0; ik != nks_tot; ik++)
    {
        ofs << std::setw(15) << p_kv->kvec_c[ik].x * TWOPI_Bohr2A << std::setw(15)
            << p_kv->kvec_c[ik].y * TWOPI_Bohr2A << std::setw(15)
            << p_kv->kvec_c[ik].z * TWOPI_Bohr2A << std::endl;
    }
    // added for BZ to IBZ (actually LibRPA interface only support BZ by 2025/03/30)
    if (PARAM.inp.symmetry == "-1")
    {
        for (int ik = 0; ik != nks_tot; ++ik)
        {
            ofs << (ik + 1) << std::endl;
        }
    }
    ofs.close();
    return;
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_bands(const elecstate::ElecState* pelec)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_bands");
    if (GlobalV::MY_RANK != 0)
    {
        return;
    }
    const int nks_tot = PARAM.inp.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    const int nspin_tmp = PARAM.inp.nspin == 2 ? 2 : 1;
    std::stringstream ss;
    ss << "band_out";
    std::ofstream ofs;
    ofs.open(outdir + ss.str().c_str(), std::ios::out);
    ofs << std::fixed << std::setprecision(15);
    ofs << nks_tot << std::endl;
    ofs << nspin_tmp << std::endl;
    ofs << PARAM.inp.nbands << std::endl;
    ofs << PARAM.globalv.nlocal << std::endl;
    ofs << (pelec->eferm.ef / 2.0) << std::endl;

    for (int ik = 0; ik != nks_tot; ik++)
    {
        for (int is = 0; is != nspin_tmp; is++)
        {
            ofs << std::setw(6) << ik + 1 << std::setw(6) << is + 1 << std::endl;
            for (int ib = 0; ib != PARAM.inp.nbands; ib++)
            {
                ofs << std::setw(5) << ib + 1 << "   " << std::setw(8) << pelec->wg(ik + is * nks_tot, ib) * nks_tot
                    << std::setw(25) << pelec->ekb(ik + is * nks_tot, ib) / 2.0
                    << std::setw(25)
                    << pelec->ekb(ik + is * nks_tot, ib) * ModuleBase::Ry_to_eV << std::endl;
            }
        }
    }
    ofs.close();
    return;
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_Cs(const UnitCell& ucell, std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs_in, std::string filename)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_Cs");
    ModuleBase::timer::start("RPA_LRI", "out_Cs");

    std::stringstream ss;
    ss << filename << GlobalV::MY_RANK << ".txt";
    std::ofstream ofs;
    ofs.open(outdir + ss.str().c_str(), std::ios::out);
    ofs << ucell.nat << "    " << 0 << std::endl;
    ofs << std::fixed << std::setprecision(15);
    for (auto& Ip: Cs_in)
    {
        size_t I = Ip.first;
        size_t i_num = ucell.atoms[ucell.iat2it[I]].nw;
        for (auto& JPp: Ip.second)
        {
            size_t J = JPp.first.first;
            auto R = JPp.first.second;
            auto& tmp_Cs = JPp.second;
            size_t j_num = ucell.atoms[ucell.iat2it[J]].nw;

            ofs << I + 1 << "   " << J + 1 << "   " << R[0] << "   " << R[1] << "   " << R[2] << "   " << i_num
                << std::endl;
            ofs << j_num << "   " << tmp_Cs.shape[0] << std::endl;
            for (int i = 0; i != i_num; i++)
            {
                for (int j = 0; j != j_num; j++)
                {
                    for (int mu = 0; mu != tmp_Cs.shape[0]; mu++)
                    {
                        ofs << std::setw(30) << tmp_Cs(mu, i, j) << std::endl;
                    }
                }
            }
        }
    }
    ofs.close();
    ModuleBase::timer::end("RPA_LRI", "out_Cs");
    return;
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_coulomb_k(const UnitCell& ucell,
                                      std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs,
                                      std::string filename,
                                      Exx_LRI<double>* exx_lri)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_coulomb_k");
    ModuleBase::timer::start("RPA_LRI", "out_coulomb_k");

    int all_mu = 0;
    std::vector<int> mu_shift(ucell.nat);
    if (exx_lri->exx_objs.empty())
    {
        throw std::invalid_argument(std::string(__FILE__) + " line " + std::to_string(__LINE__));
    }
    const auto basis_method = exx_lri->exx_objs.count(Conv_Coulomb_Pot_K::Coulomb_Method::Center2)
        ? Conv_Coulomb_Pot_K::Coulomb_Method::Center2
        : exx_lri->exx_objs.begin()->first;
    for (int I = 0; I != ucell.nat; I++)
    {
        mu_shift[I] = all_mu;
        all_mu += exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);
    }
    const int nks_tot = PARAM.inp.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    std::stringstream ss;
    ss << filename << GlobalV::MY_RANK << ".txt";

    std::ofstream ofs;
    ofs.open(outdir + ss.str().c_str(), std::ios::out);

    ofs << nks_tot << std::endl;
    ofs << std::fixed << std::setprecision(15);
    for (auto& Ip: Vs)
    {
        auto I = Ip.first;
        size_t mu_num = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);

        for (int ik = 0; ik != nks_tot; ik++)
        {
            std::map<size_t, RI::Tensor<std::complex<double>>> Vq_k_IJ;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;

                auto R = JPp.first.second;
                if (J < I)
                {
                    continue;
                }
                RI::Tensor<std::complex<double>> tmp_VR = RI::Global_Func::convert<std::complex<double>>(JPp.second);
                const double arg = 1 * (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec))
                                   * ModuleBase::TWO_PI; // latvec
                const std::complex<double> kphase = std::complex<double>(cos(arg), sin(arg));
                if (Vq_k_IJ[J].empty())
                {
                    Vq_k_IJ[J] = RI::Tensor<std::complex<double>>({tmp_VR.shape[0], tmp_VR.shape[1]});
                }
                Vq_k_IJ[J] = Vq_k_IJ[J] + tmp_VR * kphase;
            }
            for (auto& vq_Jp: Vq_k_IJ)
            {
                auto iJ = vq_Jp.first;
                auto& vq_J = vq_Jp.second;
                size_t nu_num = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[iJ]);
                ofs << all_mu << "   " << mu_shift[I] + 1 << "   " << mu_shift[I] + mu_num << "  " << mu_shift[iJ] + 1
                    << "   " << mu_shift[iJ] + nu_num << std::endl;
                ofs << ik + 1 << "  " << p_kv->wk[ik] / 2.0 * PARAM.inp.nspin << std::endl;
                for (int i = 0; i != vq_J.data->size(); i++)
                {
                    ofs << std::setw(25) << (*vq_J.data)[i].real()
                        << std::setw(25) << (*vq_J.data)[i].imag() << std::endl;
                }
            }
        }
    }
    ofs.close();
    ModuleBase::timer::end("RPA_LRI", "out_coulomb_k");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_velocity(const UnitCell &ucell,
                                    const Grid_Driver &gd,
                                    const TwoCenterBundle &two_center_bundle,
                                    const Parallel_Orbitals &parav,/*nbasis×nbasis*/
                                    const psi::Psi<T> &psi,
                                    const elecstate::ElecState* pelec)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_velocity");
    ModuleBase::timer::start("RPA_LRI", "out_velocity");

    Parallel_2D parac;/*nbasis×nbands*/
    LR_Util::setup_2d_division(parac, parav.get_block_size(), PARAM.globalv.nlocal, PARAM.inp.nbands
        #ifdef __MPI
            , parav.blacs_ctxt
        #endif
            );

    const int nk = PARAM.inp.nspin == 2 ? p_kv->get_nks() / 2 : p_kv->get_nks();
    const int nspin_tmp = PARAM.inp.nspin == 2 ? 2 : 1;
    const int nbands = parav.get_wfc_global_nbands();
    const int nbasis = parav.get_wfc_global_nbasis();

    // nocc and nvirt dosen't matter, their sum is actually used
    std::vector<int> nocc(2, nbands);
    std::vector<int> nvirt(2, 0);

    std::vector<std::complex<double>> velocity_mo = LR_Util::cal_velocity_mo(ucell, gd, two_center_bundle,
        parav, parac, *this->p_kv, psi, nk, nspin_tmp, PARAM.globalv.nlocal, nocc, nvirt);
    if (GlobalV::MY_RANK == 0){
        // for librpa readable
        LR_Util::output_spectrum_mo_librpa(velocity_mo, outdir + "velocity_matrix",
            nk, nspin_tmp, nbands, nbasis, nbands, *this->p_kv);
        // for human readable
        // LR_Util::output_spectrum_mo(velocity_mo, PARAM.globalv.global_out_dir + "velocity_matrix_rpa.dat",
        //     pelec->ekb.c, nk, nspin_tmp, PARAM.inp.nbands, *this->p_kv);
    }
    ModuleBase::timer::end("RPA_LRI", "out_velocity");
}

// template<typename Tdata>
// void RPA_LRI<T, Tdata>::init(const MPI_Comm &mpi_comm_in)
// {
// 	if(this->info == this->exx.info)
// 	{
// 		this->lcaos = this->exx.lcaos;
// 		this->abfs = this->exx.abfs;
// 		this->abfs_ccp = this->exx.abfs_ccp;

// 		exx_lri_rpa.cv = std::move(this->exx.cv);
// 	}
// 	else
// 	{
// 		this->lcaos = ...
// 		this->abfs = ...
// 		this->abfs_ccp = ...

// 		exx_lri_rpa.cv.set_orbitals(
// 			this->lcaos, this->abfs, this->abfs_ccp,
// 			this->info.kmesh_times, this->info.ccp_rmesh_times );
// 	}

// }



// template<typename Tdata>
// void RPA_LRI<T, Tdata>::cal_rpa_ions()
// {
// 	// this->rpa_lri.set_parallel(this->mpi_comm, atoms_pos, latvec, period);

// 	if(this->info == this->exx.info)
// 		exx_lri_rpa.cv.Vws = std::move(this->exx.cv.Vws);

// 	const std::array<Tcell,Ndim> period_Vs =
// LRI_CV_Tools::cal_latvec_range<Tcell>(1+this->info.ccp_rmesh_times); const
// std::pair<std::vector<TA>,
// std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>> 		list_As_Vs
// = RI::Distribute_Equally::distribute_atoms(this->mpi_comm, atoms, period_Vs,
// 2, false);

// 	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>
// 		Vs = exx_lri_rpa.cv.cal_Vs(
// 			list_As_Vs.first, list_As_Vs.second[0],
// 			{{"writable_Vws",true}});

// 	// Vs[iat0][{iat1,cell1}]	按 (iat0,iat1) 分进程，每个进程有所有 cell1
// 	Vqs = FFT(Vs);
// 	out_Vs(Vqs);

// 	if(this->info == this->exx.info)
// 		exx_lri_rpa.cv.Cws = std::move(this->exx.cv.Cws);

// 	const std::array<Tcell,Ndim> period_Cs =
// LRI_CV_Tools::cal_latvec_range<Tcell>(2); 	const std::pair<std::vector<TA>,
// std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>> 		list_As_Cs
// = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms,
// period_Cs, 2, false);

// 	std::pair<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,
// std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>> 		Cs_dCs =
// exx_lri_rpa.cv.cal_Cs_dCs( 			list_As_Cs.first, list_As_Cs.second[0],
// 			{{"cal_dC",false},
// 			 {"writable_Cws",true}, {"writable_dCws",true},
// {"writable_Vws",false},
// {"writable_dVws",false}}); 	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>> &Cs
// = std::get<0>(Cs_dCs);

// 	out_Cs(Cs);

// 	// rpa_lri.set_Cs(Cs);
// }

#endif
