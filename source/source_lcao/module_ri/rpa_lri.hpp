//=======================
// AUTHOR : Rong Shi
// DATE :   2022-12-09
//=======================

#ifndef RPA_LRI_HPP
#define RPA_LRI_HPP
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include "source_lcao/module_ri/module_exx_symmetry/symm_rotation.h"

#include "rpa_lri.h"
#include "source_basis/module_ao/elem_basis_idx_orb.h"
#include "source_base/global_function.h"
#include "source_estate/elecstate_lcao.h"
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_restart/restart_exx_csr.h"

#if defined(__GLIBC__)
#include <malloc.h>
#endif

namespace RpaLriDetail
{
constexpr int LIBRPA_COULOMB_V1_MARKER = -20129433;
constexpr int LIBRPA_LRICOEF_V1_MARKER = -10267453;
constexpr int LIBRPA_COULOMB_V1_COMPLEX_FLAG = 1;

static_assert(sizeof(std::complex<double>) == 2 * sizeof(double),
              "LibRPA v1 Coulomb output expects complex<double> as two doubles.");

inline void trim_malloc_cache()
{
#if defined(__GLIBC__)
    malloc_trim(0);
#endif
}

inline bool debug_dump_exx_ao_enabled()
{
    const char* env = std::getenv("ABACUS_DUMP_EXX_AO");
    if (env == nullptr)
    {
        return false;
    }
    const std::string value(env);
    return !(value.empty() || value == "0" || value == "f" || value == "F"
             || value == "false" || value == "FALSE");
}

inline std::size_t coulomb_atom_pair_index(const std::size_t I, const std::size_t J, const std::size_t natoms)
{
    if (I > J)
    {
        throw std::runtime_error("LibRPA v1 Coulomb output expects upper-triangular atom pairs.");
    }
    return I * natoms - I * (I - 1) / 2 + (J - I);
}

inline void checked_write(std::ofstream& ofs, const void* data, const std::size_t bytes, const std::string& filename)
{
    const char* ptr = reinterpret_cast<const char*>(data);
    std::size_t bytes_left = bytes;
    const std::size_t max_chunk = static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max());
    while (bytes_left > 0)
    {
        const std::size_t chunk = std::min(bytes_left, max_chunk);
        ofs.write(ptr, static_cast<std::streamsize>(chunk));
        if (!ofs.good())
        {
            throw std::runtime_error("Failed to write " + filename);
        }
        ptr += chunk;
        bytes_left -= chunk;
    }
}

template <typename Value>
inline void write_scalar(std::ofstream& ofs, const Value& value, const std::string& filename)
{
    checked_write(ofs, &value, sizeof(Value), filename);
}

inline unsigned long long checked_mul_u64(const unsigned long long lhs,
                                          const unsigned long long rhs,
                                          const std::string& context)
{
    if (lhs != 0 && rhs > std::numeric_limits<unsigned long long>::max() / lhs)
    {
        throw std::runtime_error(context + " exceeds uint64_t range.");
    }
    return lhs * rhs;
}

inline std::int64_t checked_i64_from_u64(const unsigned long long value, const std::string& context)
{
    if (value > static_cast<unsigned long long>(std::numeric_limits<std::int64_t>::max()))
    {
        throw std::runtime_error(context + " exceeds int64_t range.");
    }
    return static_cast<std::int64_t>(value);
}

inline int sum_int_vector(const std::vector<int>& values)
{
    int sum = 0;
    for (const int value: values)
    {
        if (value > std::numeric_limits<int>::max() - sum)
        {
            throw std::runtime_error("Integer overflow while summing LibRPA v1 basis sizes.");
        }
        sum += value;
    }
    return sum;
}

template <typename Value>
inline double real_as_double(const Value& value)
{
    return static_cast<double>(value);
}

inline double real_as_double(const std::complex<double>& value)
{
    return value.real();
}

inline std::vector<std::vector<int>>
collect_abfs_l_nchi(const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs)
{
    std::vector<std::vector<int>> abfs_l_nchi;
    abfs_l_nchi.reserve(abfs.size());
    for (const auto& abfs_type : abfs)
    {
        std::vector<int> shell_counts;
        shell_counts.reserve(abfs_type.size());
        for (const auto& abfs_l : abfs_type)
        {
            shell_counts.push_back(static_cast<int>(abfs_l.size()));
        }
        abfs_l_nchi.push_back(std::move(shell_counts));
    }
    return abfs_l_nchi;
}

inline void append_unique_abfs_layout_candidates(
    std::vector<std::vector<std::vector<int>>>& candidates_by_type,
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs)
{
    if (abfs.empty())
    {
        return;
    }

    const auto shell_counts_by_type = collect_abfs_l_nchi(abfs);
    if (candidates_by_type.empty())
    {
        candidates_by_type.resize(shell_counts_by_type.size());
    }
    else if (candidates_by_type.size() != shell_counts_by_type.size())
    {
        throw std::runtime_error("ABF shell-layout candidates are inconsistent with the atom-type count.");
    }

    for (std::size_t itype = 0; itype < shell_counts_by_type.size(); ++itype)
    {
        const auto& shell_counts = shell_counts_by_type[itype];
        auto& layouts = candidates_by_type[itype];
        if (std::find(layouts.begin(), layouts.end(), shell_counts) == layouts.end())
        {
            layouts.push_back(shell_counts);
        }
    }
}

inline std::vector<std::string> collect_atom_type_labels(const UnitCell& ucell)
{
    std::vector<std::string> labels(static_cast<std::size_t>(ucell.ntype));
    for (int itype = 0; itype < ucell.ntype; ++itype)
    {
        labels[static_cast<std::size_t>(itype)] = ucell.atom_label[itype];
    }
    return labels;
}

inline int max_layout_lmax(const std::vector<std::vector<std::vector<int>>>& candidates_by_type)
{
    int lmax = -1;
    for (const auto& type_candidates : candidates_by_type)
    {
        for (const auto& shell_counts : type_candidates)
        {
            lmax = std::max(lmax, static_cast<int>(shell_counts.size()) - 1);
        }
    }
    return lmax;
}

template<typename Tdata>
inline bool has_valid_matrix_shape(const RI::Tensor<Tdata>& tensor)
{
    return tensor.shape.size() == 2 && tensor.shape[0] > 0 && tensor.shape[1] > 0;
}

template<typename Tdata>
inline std::map<RI_2D_Comm::TA, std::map<RI_2D_Comm::TAC, RI::Tensor<Tdata>>>
collect_local_irreducible_abf_blocks(
    const std::map<RI_2D_Comm::TA, std::map<RI_2D_Comm::TAC, RI::Tensor<Tdata>>>& period_blocks,
    const std::map<ModuleSymmetry::Tap, std::set<ModuleSymmetry::TC>>& irreducible_sector,
    std::size_t& n_skipped_irreducible_blocks)
{
    std::map<RI_2D_Comm::TA, std::map<RI_2D_Comm::TAC, RI::Tensor<Tdata>>> irreducible_blocks;
    n_skipped_irreducible_blocks = 0;
    for (const auto& irap_Rs: irreducible_sector)
    {
        const auto period_iter = period_blocks.find(irap_Rs.first.first);
        for (const auto& irR: irap_Rs.second)
        {
            const RI_2D_Comm::TAC ir_key = {irap_Rs.first.second, irR};
            if (period_iter == period_blocks.end())
            {
                ++n_skipped_irreducible_blocks;
                continue;
            }
            const auto block_iter = period_iter->second.find(ir_key);
            if (block_iter == period_iter->second.end()
                || !has_valid_matrix_shape(block_iter->second))
            {
                ++n_skipped_irreducible_blocks;
                continue;
            }
            irreducible_blocks[irap_Rs.first.first][ir_key] = block_iter->second;
        }
    }
    return irreducible_blocks;
}

inline std::size_t sum_skipped_irreducible_blocks(const MPI_Comm& mpi_comm,
                                                  const std::size_t local_count)
{
    unsigned long long global_count = static_cast<unsigned long long>(local_count);
    unsigned long long reduced_count = global_count;
    MPI_Allreduce(&global_count, &reduced_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm);
    return static_cast<std::size_t>(reduced_count);
}

}

template <typename T, typename Tdata>
RPA_LRI<T, Tdata>::~RPA_LRI() = default;

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
    if (RpaLriDetail::debug_dump_exx_ao_enabled())
    {
        const std::string file_name_exx
            = PARAM.globalv.global_out_dir + "HexxR" + std::to_string(GlobalV::MY_RANK);
        ModuleIO::write_Hexxs_csr(file_name_exx, ucell, exx_cut_coulomb->Hexxs);
    }
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
    exx_cut_coulomb.reset();
    RpaLriDetail::trim_malloc_cache();

    if (this->info.shrink_abfs_pca_thr >= 0.0)
    {
        cal_large_Cs(ucell, orb, kv);
        cal_abfs_overlap(ucell, orb, kv);
        RpaLriDetail::trim_malloc_cache();
    }
    this->output_symmetry_sidecars(ucell, kv, dm);
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

    Mix_DMk_2D mix_DMk_2D;
    this->use_spacegroup_symmetry_ = (PARAM.inp.nspin < 4 && ModuleSymmetry::Symmetry::symm_flag == 1);
    if (this->use_spacegroup_symmetry_)
        {mix_DMk_2D.set_nks(kv.get_nkstot_full() * (PARAM.inp.nspin == 2 ? 2 : 1), PARAM.globalv.gamma_only_local);}
    else
        {mix_DMk_2D.set_nks(kv.get_nks());}
        
    mix_DMk_2D.set_mixing(nullptr);
    if (this->use_spacegroup_symmetry_)
    {
        const std::array<Tcell, Ndim> period = RI_Util::get_Born_vonKarmen_period(kv);
        const auto& Rs = RI_Util::get_Born_von_Karmen_cells(period);
        this->symmetry_rotation_.find_irreducible_sector(ucell.symm, ucell.atoms, ucell.st, Rs, period, ucell.lat);
        // set Lmax of the rotation matrices to max(l_ao, l_abf), to support rotation under ABF
        this->symmetry_rotation_.set_abfs_Lmax(GlobalC::exx_info.info_ri.abfs_Lmax);
        this->symmetry_rotation_.cal_Ms(kv, ucell, *dm.get_paraV_pointer());
        mix_DMk_2D.mix(this->symmetry_rotation_.restore_dm(kv, dm.get_DMK_vector(), *dm.get_paraV_pointer()), true);
    }
    else { mix_DMk_2D.mix(dm.get_DMK_vector(), true); }
    
    const std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>>
		Ds = PARAM.globalv.gamma_only_local
        ? RI_2D_Comm::split_m2D_ktoR<Tdata>(ucell,kv, mix_DMk_2D.get_DMk_gamma_out(), *dm.get_paraV_pointer(), PARAM.inp.nspin)
        : RI_2D_Comm::split_m2D_ktoR<Tdata>(ucell,kv, mix_DMk_2D.get_DMk_k_out(), *dm.get_paraV_pointer(), PARAM.inp.nspin, this->use_spacegroup_symmetry_);
    
    // reserve exx_ccp_rmesh_times to calculate full Coulomb
    // Note: ccp_type=Hf and hybrid_alpha=1 were previously set on the global Exx_Info
    // and sync_from_global() was called, but this->info (value copy) already has the correct
    // coulomb_param from construction time, so those writes are redundant and removed.
    this->ccp_rmesh_times_ewald = this->info.ccp_rmesh_times;
    // Using rpa_ccp_rmesh_times to calculate cut Coulomb this->Vs_period
    Exx_Info_RI local_info = this->info;
    local_info.ccp_rmesh_times = PARAM.inp.rpa_ccp_rmesh_times;
    if (!exx_cut_coulomb)
        exx_cut_coulomb.reset(new Exx_LRI<double>(local_info));

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
        exx_cut_coulomb->init_spencer(mpi_comm_in, ucell, kv, orb);

    if (this->use_spacegroup_symmetry_)
    {
        // Refresh the ABF-side spherical-harmonic rotation matrices after the auxiliary basis
        // is finalized by `init_spencer()`. The earlier `cal_Ms()` call only guaranteed the AO
        // rotation blocks needed for density-matrix restoration.
        this->symmetry_rotation_.set_Cs_rotation(exx_cut_coulomb->get_abfs_nchis());
        this->symmetry_rotation_.cal_Ms(kv, ucell, *dm.get_paraV_pointer());
    }

    // cal C and V for exx
    this->output_cut_coulomb_cs(ucell, exx_cut_coulomb.get());
    // cal CVCD
    if (this->use_spacegroup_symmetry_ && PARAM.inp.exx_symmetry_realspace)
    {
        exx_cut_coulomb->cal_exx_elec(Ds, ucell, *dm.get_paraV_pointer(), &this->symmetry_rotation_);
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
    if (PARAM.inp.out_librpa_reader_version == 1)
    {
        this->out_librpa_basis_v1(ucell, exx_lri_rpa);
        this->out_coulomb_k_v1(ucell, this->Vs_period, "v1_coulomb_cut_iq_", exx_lri_rpa);
    }
    else
    {
        this->out_coulomb_k(ucell, this->Vs_period, "coulomb_cut_", exx_lri_rpa);
    }
    Vs_period.clear();
    Vs_period.swap(tmp);

    this->Cs_period = RI::RI_Tools::cal_period(Cs, period);
    this->Cs_period = exx_lri_rpa->exx_lri.post_2D.set_tensors_map2(this->Cs_period);

    if (PARAM.inp.out_librpa_reader_version == 1)
    {
        this->out_Cs_v1(ucell, this->Cs_period, "v1_Cs_data_");
    }
    else
    {
        if (GlobalC::exx_info.info_ri.shrink_abfs_pca_thr >= 0.0)
            this->out_Cs(ucell, this->Cs_period, "Cs_shrinked_data_");
        else
            this->out_Cs(ucell, this->Cs_period, "Cs_data_");
    }
    Cs_period.clear();
    Cs_period.swap(tmp);

    ModuleBase::timer::end("RPA_LRI", "output_cut_coulomb_cs");
}

template <typename T, typename Tdata>
Conv_Coulomb_Pot_K::Coulomb_Method RPA_LRI<T, Tdata>::select_coulomb_basis_method_(Exx_LRI<double>* exx_lri) const
{
    if (exx_lri == nullptr || exx_lri->exx_objs.empty())
    {
        throw std::invalid_argument("Cannot select Coulomb basis method from an empty Exx_LRI object.");
    }
    return exx_lri->exx_objs.count(Conv_Coulomb_Pot_K::Coulomb_Method::Center2)
        ? Conv_Coulomb_Pot_K::Coulomb_Method::Center2
        : exx_lri->exx_objs.begin()->first;
}

template <typename T, typename Tdata>
std::vector<int> RPA_LRI<T, Tdata>::collect_atom_naux_(const UnitCell& ucell, Exx_LRI<double>* exx_lri) const
{
    const auto basis_method = this->select_coulomb_basis_method_(exx_lri);
    std::vector<int> atom_naux(static_cast<std::size_t>(ucell.nat), 0);
    for (int I = 0; I != ucell.nat; ++I)
    {
        atom_naux[static_cast<std::size_t>(I)]
            = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);
        if (atom_naux[static_cast<std::size_t>(I)] <= 0)
        {
            throw std::runtime_error("LibRPA v1 output found a non-positive per-atom auxiliary size.");
        }
    }
    return atom_naux;
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::output_ewald_coulomb(const UnitCell& ucell, const K_Vectors& kv, const LCAO_Orbitals& orb)
{
    ModuleBase::TITLE("RPA_LRI", "output_ewald_coulomb");
    ModuleBase::timer::start("RPA_LRI", "output_ewald_coulomb");

    Exx_Info_RI local_info = this->info;
    local_info.ccp_rmesh_times = this->ccp_rmesh_times_ewald;
    if (!exx_full_coulomb)
        exx_full_coulomb.reset(new Exx_LRI<double>(local_info));

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
    if (this->use_spacegroup_symmetry_)
    {
        // Rebuild the full periodic bare Coulomb from the irreducible real-space sector so that
        // the exported `coulomb_mat_` uses the same ABF rotation convention as the sidecar files.
        this->symmetry_rotation_.set_Cs_rotation(exx_full_coulomb->get_abfs_nchis());

        std::size_t n_skipped_irreducible_blocks_local = 0;
        auto Vs_irreducible = RpaLriDetail::collect_local_irreducible_abf_blocks(
            this->Vs_period,
            this->symmetry_rotation_.get_irreducible_sector(),
            n_skipped_irreducible_blocks_local);
        if (GlobalV::NPROC > 1)
        {
            const std::set<TA> all_atoms_set(atoms.begin(), atoms.end());
            Vs_irreducible = RI_2D_Comm::comm_map2_first(
                this->mpi_comm, Vs_irreducible, all_atoms_set, all_atoms_set);
        }
        const std::size_t n_skipped_irreducible_blocks
            = RpaLriDetail::sum_skipped_irreducible_blocks(
                this->mpi_comm, n_skipped_irreducible_blocks_local);
        if (n_skipped_irreducible_blocks != 0 && GlobalV::MY_RANK == 0)
        {
            std::cout << "Warning: skipped " << n_skipped_irreducible_blocks
                      << " missing irreducible ABF bare-Coulomb blocks during symmetry restoration"
                      << std::endl;
        }
        this->Vs_period
            = this->symmetry_rotation_.restore_HR_abf(ucell.symm, ucell.atoms, ucell.st, Vs_irreducible);
        if (GlobalV::NPROC > 1)
        {
            this->Vs_period = RI_2D_Comm::comm_map2_first(this->mpi_comm, this->Vs_period, atoms00, atoms01);
        }
    }
    if (PARAM.inp.out_librpa_reader_version == 1)
    {
        this->out_librpa_basis_v1(ucell, exx_full_coulomb.get());
        this->out_coulomb_k_v1(ucell, this->Vs_period, "v1_coulomb_full_iq_", exx_full_coulomb.get());
    }
    else
    {
        this->out_coulomb_k(ucell, this->Vs_period, "coulomb_mat_", exx_full_coulomb.get());
    }
    Vs_period.clear();
    Vs_period.swap(tmp);
    Cs.clear();
    Cs.swap(tmp);

    exx_full_coulomb.reset();
    RpaLriDetail::trim_malloc_cache();

    ModuleBase::timer::end("RPA_LRI", "output_ewald_coulomb");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::cal_large_Cs(const UnitCell& ucell, const LCAO_Orbitals& orb, const K_Vectors& kv)
{
    ModuleBase::TITLE("RPA_LRI", "cal_large_Cs");
    ModuleBase::timer::start("RPA_LRI", "cal_large_Cs");
    if (!exx_cut_coulomb)
        exx_cut_coulomb.reset(new Exx_LRI<double>(this->info));
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

        this->out_coulomb_k(ucell, this->Vs_period, "coulomb_unshrinked_cut_", exx_cut_coulomb.get());
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
void RPA_LRI<T, Tdata>::output_symmetry_sidecars(const UnitCell& ucell,
                                                 const K_Vectors& kv,
                                                 const elecstate::DensityMatrix<T, Tdata>& dm)
{
    const bool exx_spacegroup_symmetry =
        (PARAM.inp.nspin < 4 && ModuleSymmetry::Symmetry::symm_flag == 1);
    if (!exx_spacegroup_symmetry)
    {
        return;
    }
    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    std::vector<std::vector<std::vector<int>>> abf_layout_candidates;
    if (GlobalC::exx_info.info_ri.shrink_abfs_pca_thr >= 0.0)
    {
        RpaLriDetail::append_unique_abfs_layout_candidates(abf_layout_candidates, this->abfs_shrink);
        RpaLriDetail::append_unique_abfs_layout_candidates(abf_layout_candidates, this->abfs);
    }
    else
    {
        RpaLriDetail::append_unique_abfs_layout_candidates(abf_layout_candidates, this->abfs);
    }

    ModuleSymmetry::Symmetry_rotation symrot;
    const std::array<Tcell, Ndim> period = RI_Util::get_Born_vonKarmen_period(kv);
    const auto& Rs = RI_Util::get_Born_von_Karmen_cells(period);
    symrot.find_irreducible_sector(ucell.symm, ucell.atoms, ucell.st, Rs, period, ucell.lat);

    const int abf_lmax = RpaLriDetail::max_layout_lmax(abf_layout_candidates);
    if (abf_lmax >= 0)
    {
        symrot.set_abfs_Lmax(abf_lmax);
    }
    else
    {
        symrot.set_abfs_Lmax(GlobalC::exx_info.info_ri.abfs_Lmax);
    }
    symrot.cal_Ms(kv, ucell, *dm.get_paraV_pointer());

    ModuleSymmetry::print_symrot_info_R(symrot, ucell.symm, ucell.lmax, Rs);
    ModuleSymmetry::print_symrot_info_k(symrot, kv, ucell);
    ModuleSymmetry::print_symrot_info_abf_k(
        symrot, kv, ucell, RpaLriDetail::collect_atom_type_labels(ucell), abf_layout_candidates);
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
    ss << filename << (GlobalV::MY_RANK + 1) << ".txt";

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

    // out_ri_tensor("olp_ss.txt", olp_q_ss, 0.);
    // Inverse of overlap(q)
    inverse_olp(ucell, olp_q_ss, index_abfs_s);
    // out_ri_tensor("olp_ss_inv.txt", olp_q_ss, 0.);
    // out_ri_tensor("olp_s.txt", olp_q_s, 0.);
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
    const int nbasis = parav.get_wfc_global_nbasis();
    const std::size_t values_per_iw = static_cast<std::size_t>(nbands) * npsin_tmp;

#ifdef __MPI
    const MPI_Comm mpi_comm = parav.comm();
    const int mpi_rank = parav.get_coord_row() * parav.get_dim1() + parav.get_coord_col();
    const int mpi_size = parav.get_dim0() * parav.get_dim1();
    const auto check_mpi = [mpi_comm](const int local_error, const std::string& context) {
        const int local_failed = local_error == MPI_SUCCESS ? 0 : 1;
        int any_failed = 0;
        if (MPI_Allreduce(&local_failed, &any_failed, 1, MPI_INT, MPI_MAX, mpi_comm) != MPI_SUCCESS
            || any_failed != 0)
        {
            throw std::runtime_error(context);
        }
    };
    // 1. Assign a contiguous basis interval to each rank. Rank order also
    // gives the order of the corresponding text blocks in the output file.
    std::vector<int> output_basis_counts(mpi_size, nbasis / mpi_size);
    for (int ip = 0; ip < nbasis % mpi_size; ++ip)
    {
        ++output_basis_counts[ip];
    }
    std::vector<int> output_basis_offsets(mpi_size + 1, 0);
    std::vector<int> output_basis_owner(nbasis);
    for (int ip = 0; ip < mpi_size; ++ip)
    {
        output_basis_offsets[ip + 1]
            = output_basis_offsets[ip] + output_basis_counts[ip];
        std::fill(output_basis_owner.begin() + output_basis_offsets[ip],
                  output_basis_owner.begin() + output_basis_offsets[ip + 1],
                  ip);
    }
    const int local_nw = output_basis_counts[mpi_rank];
#else
    const int mpi_rank = 0;
    const int local_nw = nbasis;
#endif
    const std::size_t local_size = static_cast<std::size_t>(local_nw) * values_per_iw;
    std::vector<std::complex<double>> local_wfc(local_size);
#ifdef __MPI
    struct WfcPackIndex
    {
        int spin;
        int band;
        int basis;
    };

    // 2. Build the 2D block-cyclic to contiguous-basis redistribution map.
    // ncol_bands is the number of wavefunction bands stored on this rank.
    const int local_band_count = parav.ncol_bands;
    const unsigned long long send_size_wide
        = static_cast<unsigned long long>(local_band_count) * psi.get_nbasis() * npsin_tmp;
    const unsigned long long recv_size_wide
        = static_cast<unsigned long long>(local_nw) * values_per_iw;
    const unsigned long long max_alltoallv_count = std::numeric_limits<int>::max();
    check_mpi(send_size_wide <= max_alltoallv_count && recv_size_wide <= max_alltoallv_count
                  ? MPI_SUCCESS
                  : MPI_ERR_COUNT,
              "RPA eigenvector buffer exceeds MPI_Alltoallv count range.");

    // Group every locally owned (band, basis, spin) coefficient by destination rank.
    std::vector<int> send_counts(mpi_size, 0);
    for (int ib = 0; ib < nbands; ++ib)
    {
        if (parav.global2local_col(ib) < 0)
        {
            continue;
        }
        for (int ir = 0; ir < psi.get_nbasis(); ++ir)
        {
            send_counts[output_basis_owner[parav.local2global_row(ir)]] += npsin_tmp;
        }
    }
    std::vector<int> send_displacements(mpi_size, 0);
    std::vector<int> recv_counts(mpi_size);
    std::vector<int> recv_displacements(mpi_size, 0);
    for (int ip = 1; ip < mpi_size; ++ip)
    {
        send_displacements[ip] = send_displacements[ip - 1] + send_counts[ip - 1];
    }
    check_mpi(MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, mpi_comm),
              "Failed to exchange RPA eigenvector redistribution counts.");
    for (int ip = 1; ip < mpi_size; ++ip)
    {
        recv_displacements[ip] = recv_displacements[ip - 1] + recv_counts[ip - 1];
    }
    const int send_size = static_cast<int>(send_size_wide);
    const int recv_size = static_cast<int>(recv_size_wide);
    check_mpi(recv_displacements.back() + recv_counts.back() == recv_size
                  ? MPI_SUCCESS
                  : MPI_ERR_COUNT,
              "RPA eigenvector Alltoallv redistribution map is inconsistent.");
    // Store source-local coordinates and destination-local flat indices in the
    // same packing order. The map is k-independent, so exchange target indices once.
    std::vector<WfcPackIndex> pack_indices(send_size);
    std::vector<unsigned long long> send_targets(send_size);
    std::vector<unsigned long long> recv_targets(recv_size);
    std::vector<int> send_positions = send_displacements;
    for (int ib = 0; ib < nbands; ++ib)
    {
        const int ib_local = parav.global2local_col(ib);
        if (ib_local < 0)
        {
            continue;
        }
        for (int ir = 0; ir < psi.get_nbasis(); ++ir)
        {
            const int iw = parav.local2global_row(ir);
            const int destination = output_basis_owner[iw];
            for (int is = 0; is < npsin_tmp; ++is)
            {
                const int position = send_positions[destination]++;
                pack_indices[position] = {is, ib_local, ir};
                send_targets[position]
                    = static_cast<unsigned long long>(iw - output_basis_offsets[destination])
                          * values_per_iw
                      + ib * npsin_tmp + is;
            }
        }
    }
    std::complex<double> dummy(0.0, 0.0);
    unsigned long long dummy_target = 0;
    check_mpi(MPI_Alltoallv(send_size > 0 ? send_targets.data() : &dummy_target,
                            send_counts.data(),
                            send_displacements.data(),
                            MPI_UNSIGNED_LONG_LONG,
                            recv_size > 0 ? recv_targets.data() : &dummy_target,
                            recv_counts.data(),
                            recv_displacements.data(),
                            MPI_UNSIGNED_LONG_LONG,
                            mpi_comm),
              "Failed to exchange RPA eigenvector destination indices.");
    std::vector<std::complex<double>> send_wfc(send_size);
    std::vector<std::complex<double>> recv_wfc(recv_size);
#endif
    std::string output_buffer;
    constexpr std::size_t output_line_bytes = 61; // required width: 2 * 30 columns plus '\n'
    output_buffer.reserve(local_size * output_line_bytes + 32); // 32 bytes for the k-point index and newline

    // 3. Open the unified output file once for all k points.
    const std::string filename = outdir + "KS_eigenvector.txt";
#ifdef __MPI
    MPI_File file = MPI_FILE_NULL;
    unsigned long long total_file_bytes = 0; // start offset of the current k point
    const unsigned long long max_mpi_offset
        = static_cast<unsigned long long>(std::numeric_limits<MPI_Offset>::max());
    const char dummy_buffer = '\0';
    check_mpi(MPI_File_open(mpi_comm,
                            filename.c_str(),
                            MPI_MODE_CREATE | MPI_MODE_WRONLY,
                            MPI_INFO_NULL,
                            &file),
              "Failed to open " + filename + ".");
    check_mpi(MPI_File_set_size(file, 0), "Failed to truncate " + filename + ".");
#else
    std::ofstream ofs(filename.c_str(), std::ios::out);
#endif
    // 4. Process k points in file order; all ranks participate in each iteration.
    for (int ik = 0; ik < nks_tot; ik++)
    {
#ifdef __MPI
        // 4.1 Pack and redistribute all bands and spins for this k point.
        for (int index = 0; index < send_size; ++index)
        {
            const WfcPackIndex& source = pack_indices[index];
            send_wfc[index] = psi(ik + nks_tot * source.spin, source.band, source.basis);
        }
        check_mpi(MPI_Alltoallv(send_size > 0 ? send_wfc.data() : &dummy,
                                send_counts.data(),
                                send_displacements.data(),
                                MPI_DOUBLE_COMPLEX,
                                recv_size > 0 ? recv_wfc.data() : &dummy,
                                recv_counts.data(),
                                recv_displacements.data(),
                                MPI_DOUBLE_COMPLEX,
                                mpi_comm),
                  "Failed to redistribute RPA eigenvectors with MPI_Alltoallv.");
        for (int index = 0; index < recv_size; ++index)
        {
            local_wfc[static_cast<std::size_t>(recv_targets[index])] = recv_wfc[index];
        }
#else
        for (int is = 0; is < npsin_tmp; is++)
        {
            for (int ib = 0; ib < nbands; ++ib)
            {
                for (int iw = 0; iw < psi.get_nbasis(); ++iw)
                {
                    local_wfc[iw * values_per_iw + ib * npsin_tmp + is]
                        = psi(ik + nks_tot * is, ib, iw);
                }
            }
        }
#endif
        // 4.2 Format this rank's contiguous basis interval in basis-band-spin order.
        output_buffer.clear();
        if (mpi_rank == 0)
        {
            output_buffer += std::to_string(ik + 1);
            output_buffer += '\n';
        }
        char line_buffer[128];
        for (int iw = 0; iw < local_nw; ++iw)
        {
            for (int ib = 0; ib < nbands; ++ib)
            {
                for (int is = 0; is < npsin_tmp; is++)
                {
                    const std::complex<double>& value
                        = local_wfc[static_cast<std::size_t>(iw) * values_per_iw + ib * npsin_tmp + is];
                    const int line_length = std::snprintf(line_buffer,
                                                          sizeof(line_buffer),
                                                          "%30.15f%30.15f\n",
                                                          value.real(),
                                                          value.imag());
                    if (line_length != static_cast<int>(output_line_bytes))
                    {
                        throw std::runtime_error("Failed to format an RPA eigenvector value.");
                    }
                    output_buffer.append(line_buffer, static_cast<std::size_t>(line_length));
                }
            }
        }

#ifdef __MPI
        // 4.3 Fixed-width value lines make each rank's file offset deterministic;
        // only the k-point header length varies.
        const unsigned long long local_bytes
            = static_cast<unsigned long long>(output_buffer.size());
        const unsigned long long header_bytes = std::to_string(ik + 1).size() + 1;
        const unsigned long long bytes_per_basis = values_per_iw * output_line_bytes;
        if (bytes_per_basis > (max_mpi_offset - header_bytes) / std::max(nbasis, 1)
            || header_bytes + static_cast<unsigned long long>(nbasis) * bytes_per_basis
                   > max_mpi_offset - total_file_bytes)
        {
            throw std::runtime_error("RPA eigenvector file layout exceeds MPI_Offset range.");
        }
        const unsigned long long kpoint_bytes
            = header_bytes + static_cast<unsigned long long>(nbasis) * bytes_per_basis;
        const unsigned long long rank_offset_bytes
            = mpi_rank == 0
                  ? 0
                  : header_bytes
                        + static_cast<unsigned long long>(output_basis_offsets[mpi_rank])
                              * bytes_per_basis;
        const unsigned long long file_offset_bytes = total_file_bytes + rank_offset_bytes;
        // Rank 0 has the largest basis slice and additionally owns the k-point header.
        const unsigned long long max_local_bytes
            = header_bytes
              + static_cast<unsigned long long>(output_basis_counts[0]) * bytes_per_basis;
        const unsigned long long max_write_count
            = static_cast<unsigned long long>(std::numeric_limits<int>::max());
        // MPI_File_write_at_all is collective, so every rank executes the same
        // number of chunks; ranks without data in a chunk pass a zero count.
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
            const unsigned long long write_offset
                = write_count > 0 ? file_offset_bytes + buffer_offset : file_offset_bytes;
            check_mpi(MPI_File_write_at_all(file,
                                            static_cast<MPI_Offset>(write_offset),
                                            write_buffer,
                                            write_count,
                                            MPI_CHAR,
                                            MPI_STATUS_IGNORE),
                      "Failed to write " + filename + ".");
        }
        total_file_bytes += kpoint_bytes;
#else
        ofs << output_buffer;
#endif
    } // ik
#ifdef __MPI
    check_mpi(MPI_File_sync(file), "Failed to sync " + filename + ".");
    check_mpi(MPI_File_close(&file), "Failed to close " + filename + ".");
#endif
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_struc(const UnitCell& ucell)
{
    if (GlobalV::MY_RANK != 0)
    {
        return;
    }
    ModuleBase::TITLE("DFT_RPA_interface", "out_struc");
    const int nks_tot = PARAM.inp.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    const ModuleBase::Matrix3 lat = ucell.latvec * ucell.lat0; // in unit of Bohr
    const ModuleBase::Matrix3 G_RPA = ucell.G * (ModuleBase::TWO_PI / ucell.lat0); // in unit of 1/Bohr
    std::ofstream ofs;
    ofs.open(ss.str().c_str(), std::ios::out);
    const auto write_scientific_triplet = [&ofs](const double x, const double y, const double z) {
        ofs << std::setw(24) << std::scientific << std::setprecision(15) << x
            << std::setw(24) << std::scientific << std::setprecision(15) << y
            << std::setw(24) << std::scientific << std::setprecision(15) << z << std::endl;
    };
    write_scientific_triplet(lat.e11, lat.e12, lat.e13);
    write_scientific_triplet(lat.e21, lat.e22, lat.e23);
    write_scientific_triplet(lat.e31, lat.e32, lat.e33);

    write_scientific_triplet(G_RPA.e11, G_RPA.e12, G_RPA.e13);
    write_scientific_triplet(G_RPA.e21, G_RPA.e22, G_RPA.e23);
    write_scientific_triplet(G_RPA.e31, G_RPA.e32, G_RPA.e33);

    ofs << ucell.nat << std::endl;
    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const double& x = direct ? ucell.atoms[it].tau[ia].x * ucell.lat0
                                     : ucell.atoms[it].tau[ia].x;
            const double& y = direct ? ucell.atoms[it].tau[ia].y * ucell.lat0
                                     : ucell.atoms[it].tau[ia].y;
            const double& z = direct ? ucell.atoms[it].tau[ia].z * ucell.lat0
                                     : ucell.atoms[it].tau[ia].z;
            ofs << std::setw(24) << std::scientific << std::setprecision(15) << x
                << std::setw(24) << std::scientific << std::setprecision(15) << y
                << std::setw(24) << std::scientific << std::setprecision(15) << z
                << std::setw(15) << (it + 1) << std::endl;
        }
    }

    ofs << p_kv->nmp[0] << std::setw(6) << p_kv->nmp[1] << std::setw(6) << p_kv->nmp[2] << std::setw(6) << std::endl;

    for (int ik = 0; ik != nks_tot; ik++)
    {
        write_scientific_triplet(p_kv->kvec_c[ik].x * TWOPI_Bohr2A,
                                 p_kv->kvec_c[ik].y * TWOPI_Bohr2A,
                                 p_kv->kvec_c[ik].z * TWOPI_Bohr2A);
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
    std::ofstream ofs;
    ofs.open(outdir + "band_out.txt", std::ios::out);
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
    ss << filename << (GlobalV::MY_RANK + 1) << ".txt";
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
void RPA_LRI<T, Tdata>::out_Cs_v1(const UnitCell& ucell,
                                  std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs_in,
                                  std::string filename)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_Cs_v1");
    ModuleBase::timer::tick("RPA_LRI", "out_Cs_v1");

    struct CsRecord
    {
        int ia1 = 0;
        int ia2 = 0;
        int R[3] = {0, 0, 0};
        double max_abs = 0.0;
        std::int64_t offset = 0;
        RI::Tensor<Tdata>* tensor = nullptr;
        int nw1 = 0;
        int nw2 = 0;
        int naux = 0;
    };

    std::vector<CsRecord> records;
    records.reserve(Cs_in.size());
    for (auto& Ip: Cs_in)
    {
        const int I = static_cast<int>(Ip.first);
        const int i_num = ucell.atoms[ucell.iat2it[I]].nw;
        for (auto& JPp: Ip.second)
        {
            const int J = static_cast<int>(JPp.first.first);
            auto& tmp_Cs = JPp.second;
            if (tmp_Cs.shape.size() != 3 || tmp_Cs.shape[0] <= 0 || tmp_Cs.shape[1] <= 0 || tmp_Cs.shape[2] <= 0)
            {
                continue;
            }
            const int j_num = ucell.atoms[ucell.iat2it[J]].nw;
            if (static_cast<int>(tmp_Cs.shape[1]) != i_num || static_cast<int>(tmp_Cs.shape[2]) != j_num)
            {
                throw std::runtime_error("LibRPA v1 Cs output encountered an inconsistent tensor shape.");
            }
            CsRecord record;
            record.ia1 = I + 1;
            record.ia2 = J + 1;
            record.R[0] = JPp.first.second[0];
            record.R[1] = JPp.first.second[1];
            record.R[2] = JPp.first.second[2];
            record.tensor = &tmp_Cs;
            record.nw1 = i_num;
            record.nw2 = j_num;
            record.naux = static_cast<int>(tmp_Cs.shape[0]);
            for (int i = 0; i != record.nw1; ++i)
            {
                for (int j = 0; j != record.nw2; ++j)
                {
                    for (int mu = 0; mu != record.naux; ++mu)
                    {
                        record.max_abs = std::max(record.max_abs, std::abs(RpaLriDetail::real_as_double(tmp_Cs(mu, i, j))));
                    }
                }
            }
            records.push_back(record);
        }
    }

    const std::int64_t nblocks = static_cast<std::int64_t>(records.size());
    const std::int64_t record_bytes = static_cast<std::int64_t>(5 * sizeof(std::int32_t)
        + sizeof(double) + sizeof(std::int64_t));
    std::int64_t offset = static_cast<std::int64_t>(3 * sizeof(std::int32_t) + 2 * sizeof(std::int64_t))
        + nblocks * record_bytes;
    for (auto& record: records)
    {
        record.offset = offset;
        const unsigned long long nw_product = RpaLriDetail::checked_mul_u64(
            static_cast<unsigned long long>(record.nw1),
            static_cast<unsigned long long>(record.nw2),
            "LibRPA v1 Cs block size");
        const unsigned long long values = RpaLriDetail::checked_mul_u64(
            nw_product,
            static_cast<unsigned long long>(record.naux),
            "LibRPA v1 Cs block size");
        const unsigned long long bytes = RpaLriDetail::checked_mul_u64(
            values,
            static_cast<unsigned long long>(sizeof(double)),
            "LibRPA v1 Cs block size");
        offset += RpaLriDetail::checked_i64_from_u64(bytes, "LibRPA v1 Cs block size");
    }

    std::stringstream ss;
    ss << filename << GlobalV::MY_RANK << ".txt";
    const std::string out_name = ss.str();
    std::ofstream ofs(out_name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open " + out_name);
    }

    const std::int32_t marker = RpaLriDetail::LIBRPA_LRICOEF_V1_MARKER;
    const std::int32_t natom = static_cast<std::int32_t>(ucell.nat);
    const std::int32_t ncell = 0;
    RpaLriDetail::write_scalar(ofs, marker, out_name);
    RpaLriDetail::write_scalar(ofs, natom, out_name);
    RpaLriDetail::write_scalar(ofs, ncell, out_name);
    RpaLriDetail::write_scalar(ofs, nblocks, out_name);
    RpaLriDetail::write_scalar(ofs, nblocks, out_name);
    for (const auto& record: records)
    {
        const std::int32_t ia1 = record.ia1;
        const std::int32_t ia2 = record.ia2;
        const std::int32_t r0 = record.R[0];
        const std::int32_t r1 = record.R[1];
        const std::int32_t r2 = record.R[2];
        RpaLriDetail::write_scalar(ofs, ia1, out_name);
        RpaLriDetail::write_scalar(ofs, ia2, out_name);
        RpaLriDetail::write_scalar(ofs, r0, out_name);
        RpaLriDetail::write_scalar(ofs, r1, out_name);
        RpaLriDetail::write_scalar(ofs, r2, out_name);
        RpaLriDetail::write_scalar(ofs, record.max_abs, out_name);
        RpaLriDetail::write_scalar(ofs, record.offset, out_name);
    }
    for (const auto& record: records)
    {
        for (int i = 0; i != record.nw1; ++i)
        {
            for (int j = 0; j != record.nw2; ++j)
            {
                for (int mu = 0; mu != record.naux; ++mu)
                {
                    const double value = RpaLriDetail::real_as_double((*record.tensor)(mu, i, j));
                    RpaLriDetail::write_scalar(ofs, value, out_name);
                }
            }
        }
    }
    ofs.close();
    ModuleBase::timer::tick("RPA_LRI", "out_Cs_v1");
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
    const auto basis_method = this->select_coulomb_basis_method_(exx_lri);
    for (int I = 0; I != ucell.nat; I++)
    {
        mu_shift[I] = all_mu;
        all_mu += exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);
    }
    const int nks_tot = PARAM.inp.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    std::stringstream ss;
    ss << filename << (GlobalV::MY_RANK + 1) << ".txt";

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
                if (!RpaLriDetail::has_valid_matrix_shape(JPp.second))
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
void RPA_LRI<T, Tdata>::out_librpa_basis_v1(const UnitCell& ucell, Exx_LRI<double>* exx_lri)
{
    if (GlobalV::MY_RANK != 0)
    {
        return;
    }
    ModuleBase::TITLE("DFT_RPA_interface", "out_librpa_basis_v1");

    const auto basis_method = this->select_coulomb_basis_method_(exx_lri);
    std::vector<int> type_naux(static_cast<std::size_t>(ucell.ntype), 0);
    int total_wfc = 0;
    int total_aux = 0;
    for (int it = 0; it != ucell.ntype; ++it)
    {
        type_naux[static_cast<std::size_t>(it)] = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(it);
        if (type_naux[static_cast<std::size_t>(it)] <= 0)
        {
            throw std::runtime_error("LibRPA v1 basis_out found a non-positive per-type auxiliary size.");
        }
        total_wfc += ucell.atoms[it].nw * ucell.atoms[it].na;
        total_aux += type_naux[static_cast<std::size_t>(it)] * ucell.atoms[it].na;
    }

    std::ofstream ofs("basis_out", std::ios::out | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open basis_out");
    }
    ofs << std::setw(10) << ucell.ntype
        << std::setw(10) << total_wfc
        << std::setw(10) << total_aux
        << "    fallback" << std::endl;
    for (int it = 0; it != ucell.ntype; ++it)
    {
        ofs << std::setw(10) << it + 1
            << std::setw(10) << ucell.atoms[it].nw
            << std::setw(10) << type_naux[static_cast<std::size_t>(it)]
            << std::endl;
    }
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_coulomb_k_v1(const UnitCell& ucell,
                                         std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs,
                                         std::string filename,
                                         Exx_LRI<double>* exx_lri)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_coulomb_k_v1");
    ModuleBase::timer::tick("RPA_LRI", "out_coulomb_k_v1");

    const auto basis_method = this->select_coulomb_basis_method_(exx_lri);
    const auto atom_naux = this->collect_atom_naux_(ucell, exx_lri);
    const int all_mu = RpaLriDetail::sum_int_vector(atom_naux);
    const int nks_tot = PARAM.inp.nspin == 2 ? static_cast<int>(p_kv->get_nks()) / 2 : p_kv->get_nks();
    const std::size_t natoms = static_cast<std::size_t>(ucell.nat);

    struct V1Block
    {
        int pair_index = 0;
        int I = 0;
        int J = 0;
        std::int64_t offset = 0;
        RI::Tensor<std::complex<double>> tensor;
    };

    for (int ik = 0; ik != nks_tot; ++ik)
    {
        std::vector<V1Block> blocks;
        for (auto& Ip: Vs)
        {
            const int I = static_cast<int>(Ip.first);
            const int mu_num = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);
            std::map<size_t, RI::Tensor<std::complex<double>>> Vq_k_IJ;
            for (auto& JPp: Ip.second)
            {
                const int J = static_cast<int>(JPp.first.first);
                if (J < I)
                {
                    continue;
                }
                if (!RpaLriDetail::has_valid_matrix_shape(JPp.second))
                {
                    continue;
                }
                RI::Tensor<std::complex<double>> tmp_VR = RI::Global_Func::convert<std::complex<double>>(JPp.second);
                const auto R = JPp.first.second;
                const double arg = (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec))
                    * ModuleBase::TWO_PI;
                const std::complex<double> kphase = std::complex<double>(std::cos(arg), std::sin(arg));
                if (Vq_k_IJ[J].empty())
                {
                    Vq_k_IJ[J] = RI::Tensor<std::complex<double>>({tmp_VR.shape[0], tmp_VR.shape[1]});
                }
                Vq_k_IJ[J] = Vq_k_IJ[J] + tmp_VR * kphase;
            }
            for (auto& vq_Jp: Vq_k_IJ)
            {
                const int J = static_cast<int>(vq_Jp.first);
                auto& vq_J = vq_Jp.second;
                const int nu_num = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[J]);
                if (static_cast<int>(vq_J.shape[0]) != mu_num || static_cast<int>(vq_J.shape[1]) != nu_num)
                {
                    throw std::runtime_error("LibRPA v1 Coulomb output encountered an inconsistent tensor shape.");
                }
                V1Block block;
                block.pair_index = static_cast<int>(RpaLriDetail::coulomb_atom_pair_index(
                    static_cast<std::size_t>(I), static_cast<std::size_t>(J), natoms));
                block.I = I;
                block.J = J;
                block.tensor = std::move(vq_J);
                blocks.push_back(std::move(block));
            }
        }

        std::sort(blocks.begin(), blocks.end(), [](const V1Block& lhs, const V1Block& rhs) {
            return lhs.pair_index < rhs.pair_index;
        });
        if (blocks.empty())
        {
            continue;
        }
        for (std::size_t ib = 1; ib < blocks.size(); ++ib)
        {
            if (blocks[ib - 1].pair_index == blocks[ib].pair_index)
            {
                throw std::runtime_error("LibRPA v1 Coulomb output found duplicate atom-pair blocks on one MPI rank.");
            }
        }

        const std::int64_t nblocks = static_cast<std::int64_t>(blocks.size());
        const std::int64_t header_size = static_cast<std::int64_t>(6 * sizeof(std::int32_t))
            + static_cast<std::int64_t>(ucell.nat * sizeof(std::int32_t))
            + nblocks * static_cast<std::int64_t>(sizeof(std::int32_t) + sizeof(std::int64_t));
        std::int64_t offset = header_size;
        for (auto& block: blocks)
        {
            block.offset = offset;
            const unsigned long long values = RpaLriDetail::checked_mul_u64(
                static_cast<unsigned long long>(atom_naux[static_cast<std::size_t>(block.I)]),
                static_cast<unsigned long long>(atom_naux[static_cast<std::size_t>(block.J)]),
                "LibRPA v1 Coulomb block size");
            const unsigned long long bytes = RpaLriDetail::checked_mul_u64(
                values,
                static_cast<unsigned long long>(sizeof(std::complex<double>)),
                "LibRPA v1 Coulomb block size");
            offset += RpaLriDetail::checked_i64_from_u64(bytes, "LibRPA v1 Coulomb block size");
        }

        std::stringstream ss;
        ss << filename << ik + 1 << "_rank" << GlobalV::MY_RANK << ".dat";
        const std::string out_name = ss.str();
        std::ofstream ofs(out_name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
        if (!ofs.good())
        {
            throw std::runtime_error("Failed to open " + out_name);
        }

        const std::int32_t marker = RpaLriDetail::LIBRPA_COULOMB_V1_MARKER;
        const std::int32_t iq = ik + 1;
        const std::int32_t naux = all_mu;
        const std::int32_t value_flag = RpaLriDetail::LIBRPA_COULOMB_V1_COMPLEX_FLAG;
        const std::int32_t natom = ucell.nat;
        const std::int32_t nblock_i32 = static_cast<std::int32_t>(nblocks);
        RpaLriDetail::write_scalar(ofs, marker, out_name);
        RpaLriDetail::write_scalar(ofs, iq, out_name);
        RpaLriDetail::write_scalar(ofs, naux, out_name);
        RpaLriDetail::write_scalar(ofs, value_flag, out_name);
        RpaLriDetail::write_scalar(ofs, natom, out_name);
        RpaLriDetail::write_scalar(ofs, nblock_i32, out_name);
        for (const int atom_aux: atom_naux)
        {
            const std::int32_t atom_aux_i32 = atom_aux;
            RpaLriDetail::write_scalar(ofs, atom_aux_i32, out_name);
        }
        for (const auto& block: blocks)
        {
            const std::int32_t pair_index = block.pair_index;
            RpaLriDetail::write_scalar(ofs, pair_index, out_name);
            RpaLriDetail::write_scalar(ofs, block.offset, out_name);
        }
        for (const auto& block: blocks)
        {
            const std::size_t nvalues = block.tensor.get_shape_all();
            if (nvalues == 0)
            {
                throw std::runtime_error("LibRPA v1 Coulomb output encountered an empty tensor payload.");
            }
            RpaLriDetail::checked_write(
                ofs,
                block.tensor.ptr(),
                nvalues * sizeof(std::complex<double>),
                out_name);
        }
        ofs.close();
    }

    ModuleBase::timer::tick("RPA_LRI", "out_coulomb_k_v1");
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

// 	// Vs[iat0][{iat1,cell1}]	distributed across processes by (iat0,iat1); each process holds all cell1
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
