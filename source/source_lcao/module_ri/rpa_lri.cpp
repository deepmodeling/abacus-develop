#include "rpa_lri.h"

#include "exx_lri.h"
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
RPA_LRI<T, Tdata>::~RPA_LRI()
{
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::postSCF(const UnitCell& ucell,
                                const MPI_Comm& mpi_comm_in,
                                const module_dm::DensityMatrix<T, Tdata>& dm,
                                const elecstate::ElecState* pelec,
                                const K_Vectors& kv,
                                const LCAO_Orbitals& orb,
                                const Parallel_Orbitals& parav,
                                const psi::Psi<T>& psi)
try
{
    ModuleBase::TITLE("RPA_LRI", "postSCF");
    ModuleBase::timer::start("RPA_LRI", "postSCF");
    ModuleBase::GlobalFunc::MAKE_DIR(this->runtime.input.rpa_outdir);
    this->ccp_rmesh_times_cut = this->runtime.input.rpa_ccp_rmesh_times;
    this->ccp_rmesh_times_ewald = this->info.ccp_rmesh_times; // should be `exx_ccp_rmesh_times`

    this->cal_postSCF_exx(dm, mpi_comm_in, ucell, kv, orb, parav);
    this->init(mpi_comm_in, kv, orb.cutoffs());
    this->out_bands(pelec);
    this->out_eigen_vector(parav, psi);
    this->out_struc(ucell);
    this->out_bz_sampling(ucell);

    std::cout << "rpa_pca_threshold: " << this->info.pca_threshold << std::endl;
    std::cout << "rpa_ccp_rmesh_times_cut: " << this->ccp_rmesh_times_cut << std::endl;
    std::cout << "rpa_ccp_rmesh_times_ewald: " << this->ccp_rmesh_times_ewald << std::endl;
    std::cout << "rpa_lcao_exx(Ha): " << std::fixed << std::setprecision(15) << exx_cut_coulomb->Eexx / 2.0
              << std::endl;

    std::cout << "etxc(Ha): " << std::fixed << std::setprecision(15) << pelec->f_en.etxc / 2.0 << std::endl;
    std::cout << "etot(Ha): " << std::fixed << std::setprecision(15) << pelec->f_en.etot / 2.0 << std::endl;
    std::cout << "Etot_without_rpa(Ha): " << std::fixed << std::setprecision(15)
              << (pelec->f_en.etot - pelec->f_en.etxc + exx_cut_coulomb->Eexx) / 2.0 << std::endl;
    exx_cut_coulomb.reset();
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
catch (const std::exception& error)
{
    RpaLriDetail::abort_output(mpi_comm_in, error.what());
}
catch (...)
{
    RpaLriDetail::abort_output(mpi_comm_in, "Unknown exception in postSCF");
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
void RPA_LRI<T, Tdata>::cal_postSCF_exx(const module_dm::DensityMatrix<T, Tdata>& dm,
                                        const MPI_Comm& mpi_comm_in,
                                        const UnitCell& ucell,
                                        const K_Vectors& kv,
                                        const LCAO_Orbitals& orb,
                                        const Parallel_Orbitals& parav)
{
    ModuleBase::TITLE("RPA_LRI", "cal_postSCF_exx");
    ModuleBase::timer::start("RPA_LRI", "cal_postSCF_exx");

    this->mpi_comm = mpi_comm_in;
    this->p_kv = &kv;
    this->orb_cutoff_ = orb.cutoffs();

    Mix_DMk_2D<T> mix_DMk_2D;
    bool exx_spacegroup_symmetry = (this->runtime.input.nspin < 4 && ModuleSymmetry::Symmetry::symm_flag == 1);
    if (exx_spacegroup_symmetry)
    {
        mix_DMk_2D.set_nks(kv.get_nkstot_nospin() * (this->runtime.input.nspin == 2 ? 2 : 1));
    }
    else
    {
        mix_DMk_2D.set_nks(kv.get_nks());
    }

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
        abfs_for_lmax = Exx_Abfs::Construct_Orbs::abfs_same_atom(ucell,
                                                                 orb,
                                                                 this->lcaos,
                                                                 this->info.kmesh_times,
                                                                 this->info.shrink_abfs_pca_thr);
        if (!this->info.files_shrink_abfs.empty())
        {
            abfs_for_lmax = Exx_Abfs::IO::construct_abfs(abfs_for_lmax,
                                                         orb,
                                                         this->info.files_shrink_abfs,
                                                         this->info.kmesh_times);
        }
    }
    else
    {
        this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs(orb, this->info.kmesh_times);
        abfs_for_lmax = Exx_Abfs::Construct_Orbs::abfs_same_atom(ucell,
                                                                 orb,
                                                                 this->lcaos,
                                                                 this->info.kmesh_times,
                                                                 this->info.pca_threshold);
        if (!this->info.files_abfs.empty())
        {
            abfs_for_lmax
                = Exx_Abfs::IO::construct_abfs(abfs_for_lmax, orb, this->info.files_abfs, this->info.kmesh_times);
        }
    }

    ModuleSymmetry::Symmetry_rotation symrot;
    if (exx_spacegroup_symmetry)
    {
        const std::array<Tcell, Ndim> period = RI_Util::get_Born_vonKarmen_period(kv);
        const auto& Rs = RI_Util::get_Born_von_Karmen_cells(period);
        symrot.find_irreducible_sector(ucell.symm,
                                       ucell.atoms,
                                       ucell.st,
                                       Rs,
                                       period,
                                       ucell.lat,
                                       this->runtime.global_out_dir);
        // set Lmax of the rotation matrices to max(l_ao, l_abf), to support rotation under ABF
        // NOTE: Using Exx_Abfs::Construct_Orbs::get_Lmax() to compute Lmax from the actual ABFs
        // instead of relying on exx_cut_coulomb->abfs_Lmax() (not yet initialized) or
        // this->info.abfs_Lmax (defaults to 0). This ensures correct Lmax for symmetry rotation.
        symrot.set_abfs_Lmax(Exx_Abfs::Construct_Orbs::get_Lmax(abfs_for_lmax));
        symrot.cal_Ms(kv, ucell, parav, this->runtime.input.nspin);
        // output Ts (symrot_R.txt) and Ms (symrot_k.txt)
        ModuleSymmetry::print_symrot_info_R(symrot, ucell.symm, ucell.lmax, Rs);
        ModuleSymmetry::print_symrot_info_k(symrot, kv, ucell);
        mix_DMk_2D.mix(symrot.restore_dm(kv, dm.get_dmk_vec(), parav), true);
    }
    else
    {
        mix_DMk_2D.mix(dm.get_dmk_vec(), true);
    }

    const std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>> Ds
        = RI_2D_Comm::split_m2D_ktoR<Tdata>(ucell,
                                            kv,
                                            mix_DMk_2D.get_DMk_out(),
                                            parav,
                                            this->runtime.input.nspin,
                                            exx_spacegroup_symmetry);

    if (!exx_cut_coulomb)
    {
        Exx_Info_RI local_info = this->info;
        local_info.ccp_rmesh_times = this->ccp_rmesh_times_cut;
        Exx_LRI<double>* new_exx = new Exx_LRI<double>(local_info);
        exx_cut_coulomb.reset(new_exx);
    }

    if (this->info.shrink_abfs_pca_thr >= 0.0)
    {
        // NOTE: Reuse abfs_for_lmax constructed earlier to avoid redundant ABFs construction.
        // This ensures consistency between the ABFs used for Lmax calculation and the ABFs
        // used for actual EXX computation.
        this->abfs_shrink = abfs_for_lmax;
        Exx_Abfs::Construct_Orbs::print_orbs_size(ucell, abfs_shrink, this->runtime.log);
        exx_cut_coulomb->init_spencer(mpi_comm_in, ucell, kv, orb, abfs_shrink);
    }
    else
        // NOTE: Reuse abfs_for_lmax constructed earlier to avoid redundant ABFs construction.
        exx_cut_coulomb->init_spencer(mpi_comm_in, ucell, kv, orb, abfs_for_lmax);
    // cal C and V for exx
    Exx_LRI<double>* cut_coulomb = exx_cut_coulomb.get();
    this->output_cut_coulomb_cs(ucell, cut_coulomb);
    // cal CVCD
    if (exx_spacegroup_symmetry && this->runtime.input.exx_symmetry_realspace)
    {
        exx_cut_coulomb->cal_exx_elec(Ds, ucell, parav, &symrot);
    }
    else
    {
        exx_cut_coulomb->cal_exx_elec(Ds, ucell, parav);
    }
    // cout<<"postSCF_Eexx: "<<exx_lri_rpa.Eexx<<endl;
    ModuleBase::timer::end("RPA_LRI", "cal_postSCF_exx");
}

template void RPA_LRI<double, double>::postSCF(const UnitCell& ucell,
                                               const MPI_Comm& mpi_comm_in,
                                               const module_dm::DensityMatrix<double, double>& dm,
                                               const elecstate::ElecState* pelec,
                                               const K_Vectors& kv,
                                               const LCAO_Orbitals& orb,
                                               const Parallel_Orbitals& parav,
                                               const psi::Psi<double>& psi);
template void RPA_LRI<std::complex<double>, double>::postSCF(
    const UnitCell& ucell,
    const MPI_Comm& mpi_comm_in,
    const module_dm::DensityMatrix<std::complex<double>, double>& dm,
    const elecstate::ElecState* pelec,
    const K_Vectors& kv,
    const LCAO_Orbitals& orb,
    const Parallel_Orbitals& parav,
    const psi::Psi<std::complex<double>>& psi);
template void RPA_LRI<double, double>::init(const MPI_Comm& mpi_comm_in,
                                            const K_Vectors& kv_in,
                                            const std::vector<double>& orb_cutoff);
template void RPA_LRI<std::complex<double>, double>::init(const MPI_Comm& mpi_comm_in,
                                                          const K_Vectors& kv_in,
                                                          const std::vector<double>& orb_cutoff);
template void RPA_LRI<double, double>::cal_postSCF_exx(const module_dm::DensityMatrix<double, double>& dm,
                                                       const MPI_Comm& mpi_comm_in,
                                                       const UnitCell& ucell,
                                                       const K_Vectors& kv,
                                                       const LCAO_Orbitals& orb,
                                                       const Parallel_Orbitals& parav);
template void RPA_LRI<std::complex<double>, double>::cal_postSCF_exx(
    const module_dm::DensityMatrix<std::complex<double>, double>& dm,
    const MPI_Comm& mpi_comm_in,
    const UnitCell& ucell,
    const K_Vectors& kv,
    const LCAO_Orbitals& orb,
    const Parallel_Orbitals& parav);
template RPA_LRI<double, double>::~RPA_LRI();
template RPA_LRI<std::complex<double>, double>::~RPA_LRI();

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_eigen_vector(const Parallel_Orbitals& parav, const psi::Psi<T>& psi)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_eigen_vector");
    if (this->runtime.input.out_librpa_ver == 1)
    {
        this->out_eigen_vector_v1(parav, psi);
    }
    else
    {
        this->out_eigen_vector_legacy(parav, psi);
    }
}

template void RPA_LRI<double, double>::out_eigen_vector(const Parallel_Orbitals& parav, const psi::Psi<double>& psi);
template void RPA_LRI<std::complex<double>, double>::out_eigen_vector(const Parallel_Orbitals& parav,
                                                                      const psi::Psi<std::complex<double>>& psi);
