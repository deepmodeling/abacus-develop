//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef EXX_LRI_HPP
#define EXX_LRI_HPP

#include "Exx_LRI.h"
#include "RI_2D_Comm.h"
#include "RI_Util.h"
#include "source_lcao/module_ri/exx_rotate_abfs.h"
#include "source_lcao/module_ri/exx_abfs-construct_orbs.h"
#include "source_lcao/module_ri/exx_abfs-io.h"
#include "source_lcao/module_ri/conv_coulomb_pot_k.h"
#include "source_base/tool_title.h"
#include "source_base/timer.h"
#include "source_lcao/module_ri/mix_dmk_2d.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_io/module_parameter/parameter.h"

#include <RI/distribute/Distribute_Equally.h>
#include <RI/global/Map_Operator-3.h>

#include <fstream>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

#if defined(__GLIBC__)
#include <malloc.h>
#endif

namespace ExxLriDetail
{
using CoulombParam
    = std::map<Conv_Coulomb_Pot_K::Coulomb_Type, std::vector<std::map<std::string, std::string>>>;

inline void trim_malloc_cache()
{
#if defined(__GLIBC__)
	malloc_trim(0);
#endif
}

inline double default_spencer_rcut(const UnitCell& ucell, const K_Vectors& kv)
{
    return std::pow(0.75 * kv.get_nkstot_full() * ucell.omega / (ModuleBase::PI), 1.0 / 3.0);
}

inline bool rotate_abfs_in_place_for_current_full_matrix(const Exx_Info::Exx_Info_RI&)
{
    // The rotated-basis Ewald split uses the full rotated ABFS for the
    // short-range channel and a separate N=0 truncation only for the
    // Gaussian long-range channel.
    return true;
}

inline int get_cal_hs_benchmark_repeat()
{
    const char* repeat_env = std::getenv("ABACUS_EXX_CALHS_BENCH_REPEAT");
    if (repeat_env == nullptr)
    {
        return 1;
    }
    const int repeat = std::atoi(repeat_env);
    return std::max(repeat, 1);
}

inline bool debug_parallel_exx_enabled()
{
    const char* debug_env = std::getenv("ABACUS_EXX_DEBUG_PARALLEL");
    return debug_env != nullptr && std::atoi(debug_env) != 0;
}

inline CoulombParam build_center2_cut_coulomb_param(const CoulombParam& coulomb_param,
                                                    const UnitCell& ucell,
                                                    const K_Vectors& kv,
                                                    bool* synthesized_rcut = nullptr)
{
    CoulombParam center2_param = RI_Util::update_coulomb_param(coulomb_param, ucell, &kv);
    const double fallback_rcut = default_spencer_rcut(ucell, kv);
    bool used_fallback_rcut = false;

    for (auto& param_list: center2_param)
    {
        if (param_list.first != Conv_Coulomb_Pot_K::Coulomb_Type::Fock)
        {
            continue;
        }
        for (auto& param: param_list.second)
        {
            auto rcut_it = param.find("Rcut");
            if (rcut_it == param.end() || rcut_it->second.empty())
            {
                param["Rcut"] = ModuleBase::GlobalFunc::TO_STRING(fallback_rcut);
                used_fallback_rcut = true;
            }
        }
    }

    if (synthesized_rcut != nullptr)
    {
        *synthesized_rcut = used_fallback_rcut;
    }
    return center2_param;
}

inline std::vector<std::vector<double>> build_rotation_rows(const std::vector<double>& moments,
                                                            const double threshold)
{
    const std::size_t n = moments.size();
    std::vector<std::vector<double>> rows;
    rows.reserve(n);

    if (n == 0)
    {
        return rows;
    }

    const double norm2 = std::inner_product(moments.begin(), moments.end(), moments.begin(), 0.0);
    const double norm = std::sqrt(norm2);
    const double ortho_tol = std::max(threshold, 10.0 * std::numeric_limits<double>::epsilon());
    if (norm <= threshold)
    {
        rows.assign(n, std::vector<double>(n, 0.0));
        for (std::size_t i = 0; i != n; ++i)
        {
            rows[i][i] = 1.0;
        }
        return rows;
    }

    rows.emplace_back(n, 0.0);
    for (std::size_t i = 0; i != n; ++i)
    {
        rows[0][i] = moments[i] / norm;
    }

    for (std::size_t basis = 0; basis != n && rows.size() != n; ++basis)
    {
        std::vector<double> candidate(n, 0.0);
        candidate[basis] = 1.0;

        for (const auto& row: rows)
        {
            const double dot = std::inner_product(candidate.begin(), candidate.end(), row.begin(), 0.0);
            for (std::size_t i = 0; i != n; ++i)
            {
                candidate[i] -= dot * row[i];
            }
        }

        const double candidate_norm2
            = std::inner_product(candidate.begin(), candidate.end(), candidate.begin(), 0.0);
        if (candidate_norm2 <= ortho_tol * ortho_tol)
        {
            continue;
        }

        const double candidate_norm = std::sqrt(candidate_norm2);
        for (double& value: candidate)
        {
            value /= candidate_norm;
        }
        rows.push_back(std::move(candidate));
    }

    if (rows.size() != n)
    {
        throw std::runtime_error("Failed to build a complete ABFS rotation basis.");
    }
    return rows;
}

inline Numerical_Orbital_Lm combine_orbital_block(
    const std::vector<Numerical_Orbital_Lm>& original_block,
    const std::vector<double>& coeffs,
    const int output_index)
{
    if (original_block.empty())
    {
        throw std::runtime_error("Cannot combine an empty ABFS block.");
    }
    if (original_block.size() != coeffs.size())
    {
        throw std::runtime_error("ABFS rotation coefficients do not match the orbital block size.");
    }

    const Numerical_Orbital_Lm& ref = original_block.front();
    const int nr = ref.getNr();
    std::vector<double> psi(nr, 0.0);
    for (std::size_t iorb = 0; iorb != original_block.size(); ++iorb)
    {
        const Numerical_Orbital_Lm& orb = original_block[iorb];
        if (orb.getNr() != nr || orb.getL() != ref.getL() || orb.getType() != ref.getType())
        {
            throw std::runtime_error("ABFS rotation requires a consistent radial grid within each (T, L) block.");
        }
        for (int ir = 0; ir != nr; ++ir)
        {
            psi[ir] += coeffs[iorb] * orb.getPsi(ir);
        }
    }

    const Numerical_Orbital_Lm& out_ref = original_block.at(output_index);
    Numerical_Orbital_Lm rotated;
    rotated.set_orbital_info(out_ref.getLabel(),
                             out_ref.getType(),
                             out_ref.getL(),
                             out_ref.getChi(),
                             nr,
                             ref.getRab(),
                             ref.getRadial(),
                             Numerical_Orbital_Lm::Psi_Type::Psi,
                             ModuleBase::GlobalFunc::VECTOR_TO_PTR(psi),
                             out_ref.getNk(),
                             out_ref.getDk(),
                             out_ref.getDruniform(),
                             false,
                             true,
                             PARAM.inp.cal_force);
    return rotated;
}

inline void rotate_abfs_by_multipole(std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs,
                                     const double threshold)
{
    const auto multipole = Exx_Abfs::Construct_Orbs::get_multipole(abfs);
    if (GlobalV::MY_RANK == 0)
    {
        std::cout << "\nRotated ABFS multipole summary" << std::endl;
        std::cout << "threshold for displayed zero moments = " << threshold << std::endl;
        std::cout << std::setprecision(16);
    }
    for (std::size_t T = 0; T != abfs.size(); ++T)
    {
        for (std::size_t L = 0; L != abfs[T].size(); ++L)
        {
            if (abfs[T][L].empty())
            {
                continue;
            }

            if (abfs[T][L].size() <= 1)
            {
                if (GlobalV::MY_RANK == 0)
                {
                    std::cout << "Atom type " << T << ", L " << L << ", N 0"
                              << ", multipole before rotation: " << multipole[T][L][0]
                              << ", multipole after rotation: " << multipole[T][L][0] << std::endl;
                }
                continue;
            }

            const std::vector<std::vector<double>> rows = build_rotation_rows(multipole[T][L], threshold);
            const std::vector<Numerical_Orbital_Lm> original_block = abfs[T][L];
            for (std::size_t N = 0; N != abfs[T][L].size(); ++N)
            {
                const double rotated_moment_raw
                    = std::inner_product(rows[N].begin(), rows[N].end(), multipole[T][L].begin(), 0.0);
                const double rotated_moment
                    = (std::abs(rotated_moment_raw) <= threshold) ? 0.0 : rotated_moment_raw;
                if (GlobalV::MY_RANK == 0)
                {
                    std::cout << "Atom type " << T << ", L " << L << ", N " << N
                              << ", multipole before rotation: " << multipole[T][L][N]
                              << ", multipole after rotation: " << rotated_moment << std::endl;
                }
                abfs[T][L][N] = combine_orbital_block(original_block, rows[N], N);
            }
        }
    }
}

inline void keep_only_leading_radial_channel(std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs)
{
    for (auto& abfs_T: abfs)
    {
        for (auto& abfs_L: abfs_T)
        {
            if (abfs_L.size() > 1)
            {
                abfs_L.resize(1);
            }
        }
    }
}

inline std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> make_leading_radial_channel_copy(
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs)
{
    auto abfs_n0 = abfs;
    keep_only_leading_radial_channel(abfs_n0);
    return abfs_n0;
}

struct AuxLongPrefixPermutation
{
    std::vector<std::vector<std::size_t>> old_to_new_by_type;
    std::vector<std::size_t> long_prefix_size_by_type;
};

inline AuxLongPrefixPermutation build_long_prefix_permutation(
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs)
{
    const ModuleBase::Element_Basis_Index::Range range = ModuleBase::Element_Basis_Index::construct_range(abfs);
    const ModuleBase::Element_Basis_Index::IndexLNM index = ModuleBase::Element_Basis_Index::construct_index(range);

    AuxLongPrefixPermutation permutation;
    permutation.old_to_new_by_type.resize(abfs.size());
    permutation.long_prefix_size_by_type.resize(abfs.size(), 0);

    for (std::size_t T = 0; T != abfs.size(); ++T)
    {
        const std::size_t full_size = index[T].count_size;
        permutation.old_to_new_by_type[T].assign(full_size, 0);

        std::vector<std::size_t> old_order;
        old_order.reserve(full_size);

        for (std::size_t L = 0; L != abfs[T].size(); ++L)
        {
            if (abfs[T][L].empty())
            {
                continue;
            }
            for (std::size_t M = 0; M != 2 * L + 1; ++M)
            {
                const std::size_t old_index = index[T][L][0][M];
                old_order.push_back(old_index);
            }
        }
        permutation.long_prefix_size_by_type[T] = old_order.size();

        for (std::size_t L = 0; L != abfs[T].size(); ++L)
        {
            for (std::size_t N = 1; N != abfs[T][L].size(); ++N)
            {
                for (std::size_t M = 0; M != 2 * L + 1; ++M)
                {
                    old_order.push_back(index[T][L][N][M]);
                }
            }
        }

        if (old_order.size() != full_size)
        {
            throw std::runtime_error("Failed to construct the rotated-ABFS long-prefix permutation.");
        }

        for (std::size_t new_index = 0; new_index != old_order.size(); ++new_index)
        {
            permutation.old_to_new_by_type[T][old_order[new_index]] = new_index;
        }
    }

    return permutation;
}

template<typename Tdata>
inline RI::Tensor<Tdata> permute_aux_tensor_rows(
    const RI::Tensor<Tdata>& tensor_in,
    const std::vector<std::size_t>& old_to_new)
{
    if (tensor_in.empty())
    {
        return tensor_in;
    }
    if (tensor_in.shape.empty() || tensor_in.shape[0] != old_to_new.size())
    {
        throw std::runtime_error("Auxiliary-row permutation does not match tensor shape.");
    }

    auto shape_out = tensor_in.shape;
    RI::Tensor<Tdata> tensor_out(shape_out);
    const std::size_t slice_size = tensor_in.get_shape_all() / tensor_in.shape[0];
    for (std::size_t old_index = 0; old_index != old_to_new.size(); ++old_index)
    {
        std::copy_n(tensor_in.ptr() + old_index * slice_size,
                    slice_size,
                    tensor_out.ptr() + old_to_new[old_index] * slice_size);
    }
    return tensor_out;
}

template<typename Tdata>
inline RI::Tensor<Tdata> permute_aux_tensor_matrix(
    const RI::Tensor<Tdata>& tensor_in,
    const std::vector<std::size_t>& row_old_to_new,
    const std::vector<std::size_t>& col_old_to_new)
{
    if (tensor_in.empty())
    {
        return tensor_in;
    }
    if (tensor_in.shape.size() != 2
        || tensor_in.shape[0] != row_old_to_new.size()
        || tensor_in.shape[1] != col_old_to_new.size())
    {
        throw std::runtime_error("Auxiliary-matrix permutation does not match tensor shape.");
    }

    RI::Tensor<Tdata> tensor_out({tensor_in.shape[0], tensor_in.shape[1]});
    for (std::size_t old_row = 0; old_row != row_old_to_new.size(); ++old_row)
    {
        const std::size_t new_row = row_old_to_new[old_row];
        for (std::size_t old_col = 0; old_col != col_old_to_new.size(); ++old_col)
        {
            tensor_out(new_row, col_old_to_new[old_col]) = tensor_in(old_row, old_col);
        }
    }
    return tensor_out;
}

template<typename Tdata>
inline void permute_aux_tensor_map_rows_by_type(
    const UnitCell& ucell,
    const std::vector<std::vector<std::size_t>>& old_to_new_per_type,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<Tdata>>>& tensors_io)
{
    for (auto& tensors_I: tensors_io)
    {
        const int type_I = ucell.iat2it[tensors_I.first];
        const auto& permutation = old_to_new_per_type.at(type_I);
        for (auto& tensors_JR: tensors_I.second)
        {
            tensors_JR.second = permute_aux_tensor_rows(tensors_JR.second, permutation);
        }
    }
}

template<typename Tdata>
inline void permute_aux_tensor_map_matrices_by_type(
    const UnitCell& ucell,
    const std::vector<std::vector<std::size_t>>& old_to_new_per_type,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<Tdata>>>& tensors_io)
{
    for (auto& tensors_I: tensors_io)
    {
        const int type_I = ucell.iat2it[tensors_I.first];
        const auto& row_permutation = old_to_new_per_type.at(type_I);
        for (auto& tensors_JR: tensors_I.second)
        {
            const int type_J = ucell.iat2it[tensors_JR.first.first];
            const auto& col_permutation = old_to_new_per_type.at(type_J);
            tensors_JR.second = permute_aux_tensor_matrix(tensors_JR.second, row_permutation, col_permutation);
        }
    }
}

template<typename Tdata>
inline void overwrite_ewald_far_field_with_moment(
    const Exx_Info::Exx_Info_RI& info,
    const UnitCell& ucell,
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs,
    const std::pair<std::vector<int>, std::vector<std::vector<std::pair<int, std::array<int, 3>>>>>& list_As_Vs,
    LRI_CV<Tdata>& cv,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<Tdata>>>& Vs_io)
{
    if (!info.coul_moment)
    {
        return;
    }
    if (PARAM.inp.cal_force || PARAM.inp.cal_stress)
    {
        throw std::invalid_argument("exx_coul_moment for Ewald currently supports energy/SCF only.");
    }

    Moment_abfs<Tdata> moment_abfs(info);
    moment_abfs.cal_multipole(abfs);
    // The overwrite targets V(a0b0), so the asymptotic condition must be based
    // on the support of the auxiliary basis itself rather than the NAO cutoff.
    const std::vector<double> abfs_cutoff = Exx_Abfs::Construct_Orbs::get_Rcut(abfs);
    moment_abfs.cal_VR(ucell,
                       abfs,
                       list_As_Vs,
                       abfs_cutoff,
                       0.0,
                       cv,
                       Vs_io,
                       true,
                       false,
                       false,
                       false);
}

template<typename Tdata>
inline double tensor_max_abs_value(const RI::Tensor<Tdata>& tensor)
{
    double max_abs = 0.0;
    for (int i = 0; i < tensor.get_shape_all(); ++i)
    {
        max_abs = std::max(max_abs, std::abs(tensor.ptr()[i]));
    }
    return max_abs;
}

struct TensorMapStats
{
    std::size_t block_count = 0;
    std::size_t element_count = 0;
    double sum_abs = 0.0;
    double max_abs = 0.0;
};

template<typename Tdata>
inline TensorMapStats tensor_map_stats(
    const std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<Tdata>>>& tensor_map)
{
    TensorMapStats stats;
    for (const auto& iat_pair: tensor_map)
    {
        (void)iat_pair;
        for (const auto& jr_pair: iat_pair.second)
        {
            const RI::Tensor<Tdata>& tensor = jr_pair.second;
            ++stats.block_count;
            stats.element_count += static_cast<std::size_t>(tensor.get_shape_all());
            for (int i = 0; i < tensor.get_shape_all(); ++i)
            {
                const double abs_value = std::abs(tensor.ptr()[i]);
                stats.sum_abs += abs_value;
                stats.max_abs = std::max(stats.max_abs, abs_value);
            }
        }
    }
    return stats;
}

inline std::string format_tensor_map_stats(const TensorMapStats& stats)
{
    std::ostringstream oss;
    oss << "blocks=" << stats.block_count
        << ", elements=" << stats.element_count
        << ", sum_abs=" << stats.sum_abs
        << ", max_abs=" << stats.max_abs;
    return oss.str();
}

inline TensorMapStats reduce_tensor_map_stats(const MPI_Comm& mpi_comm, const TensorMapStats& local_stats)
{
    TensorMapStats global_stats;
    const unsigned long long local_counts[2] = {static_cast<unsigned long long>(local_stats.block_count),
                                                static_cast<unsigned long long>(local_stats.element_count)};
    unsigned long long global_counts[2] = {0ULL, 0ULL};
    MPI_Allreduce(local_counts, global_counts, 2, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm);
    global_stats.block_count = static_cast<std::size_t>(global_counts[0]);
    global_stats.element_count = static_cast<std::size_t>(global_counts[1]);
    MPI_Allreduce(&local_stats.sum_abs, &global_stats.sum_abs, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);
    MPI_Allreduce(&local_stats.max_abs, &global_stats.max_abs, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
    return global_stats;
}

template<typename Tdata>
inline std::pair<double, double> max_abs_diff_for_tensor_map(
    const std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<Tdata>>>& lhs,
    const std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<Tdata>>>& rhs)
{
    double max_diff = 0.0;
    double max_ref = 0.0;

    const auto update_from_pair = [&](const auto& reference_map, const auto& other_map)
    {
        for (const auto& iat0_pair: reference_map)
        {
            const auto other_iat0 = other_map.find(iat0_pair.first);
            for (const auto& jr_pair: iat0_pair.second)
            {
                max_ref = std::max(max_ref, tensor_max_abs_value(jr_pair.second));
                if (other_iat0 == other_map.end())
                {
                    max_diff = std::max(max_diff, tensor_max_abs_value(jr_pair.second));
                    continue;
                }
                const auto other_jr = other_iat0->second.find(jr_pair.first);
                if (other_jr == other_iat0->second.end())
                {
                    max_diff = std::max(max_diff, tensor_max_abs_value(jr_pair.second));
                    continue;
                }
                const RI::Tensor<Tdata>& tensor_a = jr_pair.second;
                const RI::Tensor<Tdata>& tensor_b = other_jr->second;
                if (tensor_a.get_shape_all() != tensor_b.get_shape_all())
                {
                    max_diff = std::numeric_limits<double>::infinity();
                    return;
                }
                for (int i = 0; i < tensor_a.get_shape_all(); ++i)
                {
                    max_diff = std::max(max_diff, std::abs(tensor_a.ptr()[i] - tensor_b.ptr()[i]));
                }
            }
        }
    };

    update_from_pair(lhs, rhs);
    update_from_pair(rhs, lhs);
    return {max_diff, max_ref};
}
}

template<typename Tdata>
void Exx_LRI<Tdata>::init(const MPI_Comm &mpi_comm_in,
						  const UnitCell &ucell,
						  const K_Vectors &kv_in,
						  const LCAO_Orbitals& orb)
{
	ModuleBase::TITLE("Exx_LRI","init");
	ModuleBase::timer::tick("Exx_LRI", "init");

	this->mpi_comm = mpi_comm_in;
	this->p_kv = &kv_in;
	this->orb_cutoff_ = orb.cutoffs();

	this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs( orb, this->info.kmesh_times );
	Exx_Abfs::Construct_Orbs::filter_empty_orbs(this->lcaos);

	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		abfs_same_atom = Exx_Abfs::Construct_Orbs::abfs_same_atom(ucell, orb, this->lcaos, this->info.kmesh_times, this->info.pca_threshold );
	if(this->info.files_abfs.empty())
		{ this->abfs = abfs_same_atom;}
	else
		{ this->abfs = Exx_Abfs::IO::construct_abfs( abfs_same_atom, orb, this->info.files_abfs, this->info.kmesh_times ); 	}
	Exx_Abfs::Construct_Orbs::filter_empty_orbs(this->abfs);
    if (this->info.rotate_abfs && ExxLriDetail::rotate_abfs_in_place_for_current_full_matrix(this->info))
    {
        ExxLriDetail::rotate_abfs_by_multipole(this->abfs, this->info.multip_moments_threshold);
    }
	Exx_Abfs::Construct_Orbs::print_orbs_size(ucell, this->abfs, GlobalV::ofs_running);

	for( size_t T=0; T!=this->abfs.size(); ++T )
		{ GlobalC::exx_info.info_ri.abfs_Lmax = std::max( GlobalC::exx_info.info_ri.abfs_Lmax, static_cast<int>(this->abfs[T].size())-1 ); }

		this->exx_objs.clear();
		this->exx_objs_long.clear();
	    this->abfs_long_n0.clear();
		this->coulomb_settings = RI_Util::update_coulomb_settings(this->info.coulomb_param, ucell, this->p_kv);
	    this->use_rotated_n0_long_range
	        = this->info.rotate_abfs
          && this->info.coul_moment
          && this->coulomb_settings.find(Conv_Coulomb_Pot_K::Coulomb_Method::Ewald) != this->coulomb_settings.end();
    ModuleBase::Element_Basis_Index::IndexPermutation abfs_old_to_new_per_type;
    std::vector<std::size_t> abfs_long_prefix_size_per_type;
    if (this->use_rotated_n0_long_range)
    {
        if (PARAM.inp.cal_force || PARAM.inp.cal_stress)
	        {
	            throw std::invalid_argument(
	                "Rotated-ABFS split Ewald currently supports energy/SCF only.");
	        }
	        this->abfs_long_n0 = ExxLriDetail::make_leading_radial_channel_copy(this->abfs);
	        const auto permutation = ExxLriDetail::build_long_prefix_permutation(this->abfs);
	        abfs_old_to_new_per_type = permutation.old_to_new_by_type;
	        abfs_long_prefix_size_per_type = permutation.long_prefix_size_by_type;
	        if (GlobalV::MY_RANK == 0)
	        {
	            std::cout << "Rotated ABFS long-prefix sizes by type:";
	            for (std::size_t T = 0; T != abfs_long_prefix_size_per_type.size(); ++T)
	            {
	                std::cout << " T" << T << "=" << abfs_long_prefix_size_per_type[T];
	            }
	            std::cout << std::endl;
	        }
	    }

		this->MGT = std::make_shared<ORB_gaunt_table>();
		for(const auto &settings_list : this->coulomb_settings)
		{
		this->exx_objs[settings_list.first].abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp(this->abfs, settings_list.second.second, this->info.ccp_rmesh_times);
		this->exx_objs[settings_list.first].cv.set_orbitals(ucell, orb,
															this->lcaos, this->abfs, this->exx_objs[settings_list.first].abfs_ccp,
															this->info.kmesh_times, this->MGT, settings_list.second.first,
															abfs_old_to_new_per_type );
				if (settings_list.first == Conv_Coulomb_Pot_K::Coulomb_Method::Ewald)
				{
					this->exx_objs[settings_list.first].evq.init(ucell, orb,
																this->mpi_comm, this->p_kv, this->lcaos, this->abfs,
																settings_list.second.second, this->MGT, this->info.ccp_rmesh_times, this->info.kmesh_times,
																abfs_old_to_new_per_type);
					if (this->use_rotated_n0_long_range)
					{
						this->exx_objs_long[settings_list.first].abfs_ccp
						    = Conv_Coulomb_Pot_K::cal_orbs_ccp(this->abfs_long_n0,
						                                       settings_list.second.second,
						                                       this->info.ccp_rmesh_times);
						this->exx_objs_long[settings_list.first].cv.set_orbitals(ucell,
						                                                       orb,
						                                                       this->lcaos,
						                                                       this->abfs_long_n0,
						                                                       this->exx_objs_long[settings_list.first].abfs_ccp,
						                                                       this->info.kmesh_times,
						                                                       this->MGT,
						                                                       settings_list.second.first);
						this->exx_objs_long[settings_list.first].evq.init(ucell,
						                                                orb,
						                                                this->mpi_comm,
						                                                this->p_kv,
						                                                this->lcaos,
						                                                this->abfs_long_n0,
						                                                settings_list.second.second,
						                                                this->MGT,
						                                                this->info.ccp_rmesh_times,
						                                                this->info.kmesh_times);
					}
				}
			}

	ModuleBase::timer::tick("Exx_LRI", "init");
}

template<typename Tdata>
void Exx_LRI<Tdata>::init(const MPI_Comm &mpi_comm_in,
						  const UnitCell &ucell,
						  const K_Vectors &kv_in,
						  const LCAO_Orbitals& orb,
						  const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs_in)
{
	ModuleBase::TITLE("Exx_LRI","init");
	ModuleBase::timer::start("Exx_LRI", "init");

	this->mpi_comm = mpi_comm_in;
	this->p_kv = &kv_in;
	this->orb_cutoff_ = orb.cutoffs();

	this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs( orb, this->info.kmesh_times );
	Exx_Abfs::Construct_Orbs::filter_empty_orbs(this->lcaos);

	this->abfs = abfs_in;
	Exx_Abfs::Construct_Orbs::filter_empty_orbs(this->abfs);
    if (this->info.rotate_abfs && ExxLriDetail::rotate_abfs_in_place_for_current_full_matrix(this->info))
    {
        ExxLriDetail::rotate_abfs_by_multipole(this->abfs, this->info.multip_moments_threshold);
    }
	Exx_Abfs::Construct_Orbs::print_orbs_size(ucell, this->abfs, GlobalV::ofs_running);

	for( size_t T=0; T!=this->abfs.size(); ++T )
	{
		this->abfs_Lmax_ = std::max(this->abfs_Lmax_, static_cast<int>(this->abfs[T].size())-1);
	}

		this->exx_objs.clear();
		this->exx_objs_long.clear();
	    this->abfs_long_n0.clear();
		this->coulomb_settings = RI_Util::update_coulomb_settings(this->info.coulomb_param, ucell, this->p_kv);
	    this->use_rotated_n0_long_range
	        = this->info.rotate_abfs
          && this->info.coul_moment
          && this->coulomb_settings.find(Conv_Coulomb_Pot_K::Coulomb_Method::Ewald) != this->coulomb_settings.end();
    ModuleBase::Element_Basis_Index::IndexPermutation abfs_old_to_new_per_type;
    std::vector<std::size_t> abfs_long_prefix_size_per_type;
    if (this->use_rotated_n0_long_range)
    {
        if (PARAM.inp.cal_force || PARAM.inp.cal_stress)
	        {
	            throw std::invalid_argument(
	                "Rotated-ABFS split Ewald currently supports energy/SCF only.");
	        }
	        this->abfs_long_n0 = ExxLriDetail::make_leading_radial_channel_copy(this->abfs);
	        const auto permutation = ExxLriDetail::build_long_prefix_permutation(this->abfs);
	        abfs_old_to_new_per_type = permutation.old_to_new_by_type;
	        abfs_long_prefix_size_per_type = permutation.long_prefix_size_by_type;
	        if (GlobalV::MY_RANK == 0)
	        {
	            std::cout << "Rotated ABFS long-prefix sizes by type:";
	            for (std::size_t T = 0; T != abfs_long_prefix_size_per_type.size(); ++T)
	            {
	                std::cout << " T" << T << "=" << abfs_long_prefix_size_per_type[T];
	            }
	            std::cout << std::endl;
	        }
	    }

		this->MGT = std::make_shared<ORB_gaunt_table>();
		for(const auto &settings_list : this->coulomb_settings)
		{
		this->exx_objs[settings_list.first].abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp(this->abfs, settings_list.second.second, this->info.ccp_rmesh_times);
		this->exx_objs[settings_list.first].cv.set_orbitals(ucell, orb,
															this->lcaos, this->abfs, this->exx_objs[settings_list.first].abfs_ccp,
															this->info.kmesh_times, this->MGT, settings_list.second.first,
															abfs_old_to_new_per_type );
				if (settings_list.first == Conv_Coulomb_Pot_K::Coulomb_Method::Ewald)
				{
					this->exx_objs[settings_list.first].evq.init(ucell, orb,
																this->mpi_comm, this->p_kv, this->lcaos, this->abfs,
																settings_list.second.second, this->MGT, this->info.ccp_rmesh_times, this->info.kmesh_times,
																abfs_old_to_new_per_type);
					if (this->use_rotated_n0_long_range)
					{
						this->exx_objs_long[settings_list.first].abfs_ccp
						    = Conv_Coulomb_Pot_K::cal_orbs_ccp(this->abfs_long_n0,
						                                       settings_list.second.second,
						                                       this->info.ccp_rmesh_times);
						this->exx_objs_long[settings_list.first].cv.set_orbitals(ucell,
						                                                       orb,
						                                                       this->lcaos,
						                                                       this->abfs_long_n0,
						                                                       this->exx_objs_long[settings_list.first].abfs_ccp,
						                                                       this->info.kmesh_times,
						                                                       this->MGT,
						                                                       settings_list.second.first);
						this->exx_objs_long[settings_list.first].evq.init(ucell,
						                                                orb,
						                                                this->mpi_comm,
						                                                this->p_kv,
						                                                this->lcaos,
						                                                this->abfs_long_n0,
						                                                settings_list.second.second,
						                                                this->MGT,
						                                                this->info.ccp_rmesh_times,
						                                                this->info.kmesh_times);
					}
				}
			}

	ModuleBase::timer::end("Exx_LRI", "init");
}

template <typename Tdata>
void Exx_LRI<Tdata>::init_spencer(
	const MPI_Comm& mpi_comm_in,
	const UnitCell& ucell,
	const K_Vectors& kv_in,
	const LCAO_Orbitals& orb,
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs_in)
{
	ModuleBase::TITLE("Exx_LRI", "init_spencer");
	ModuleBase::timer::start("Exx_LRI", "init_spencer");

	this->mpi_comm = mpi_comm_in;
	this->p_kv = &kv_in;
	this->orb_cutoff_ = orb.cutoffs();

    this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs(orb, this->info.kmesh_times);
    Exx_Abfs::Construct_Orbs::filter_empty_orbs(this->lcaos);

    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs_same_atom
        = Exx_Abfs::Construct_Orbs::abfs_same_atom(ucell,
                                                   orb,
                                                   this->lcaos,
                                                   this->info.kmesh_times,
                                                   this->info.pca_threshold);
    if (this->info.files_abfs.empty())
    {
        this->abfs = abfs_same_atom;
    }
    else
    {
        this->abfs = Exx_Abfs::IO::construct_abfs(abfs_same_atom, orb, this->info.files_abfs, this->info.kmesh_times);
    }
    Exx_Abfs::Construct_Orbs::filter_empty_orbs(this->abfs);
    if (this->info.rotate_abfs && ExxLriDetail::rotate_abfs_in_place_for_current_full_matrix(this->info))
    {
        ExxLriDetail::rotate_abfs_by_multipole(this->abfs, this->info.multip_moments_threshold);
    }
    Exx_Abfs::Construct_Orbs::print_orbs_size(ucell, this->abfs, GlobalV::ofs_running);

	for (size_t T = 0; T != this->abfs.size(); ++T)
	{
		this->abfs_Lmax_ = std::max(this->abfs_Lmax_, static_cast<int>(this->abfs[T].size()) - 1);
	}

	    this->exx_objs.clear();
	    this->exx_objs_long.clear();
	    this->abfs_long_n0.clear();
	    this->use_rotated_n0_long_range = false;
    this->coulomb_settings.clear();
    this->coulomb_settings[Conv_Coulomb_Pot_K::Coulomb_Method::Center2]
        = std::make_pair(true,
                         ExxLriDetail::build_center2_cut_coulomb_param(
                             this->info.coulomb_param, ucell, kv_in));

    this->MGT = std::make_shared<ORB_gaunt_table>();
    const auto center2_settings = this->coulomb_settings.find(Conv_Coulomb_Pot_K::Coulomb_Method::Center2);
    if (center2_settings == this->coulomb_settings.end())
    {
        throw std::invalid_argument("Exx_LRI::init_spencer failed to prepare Center2 settings.");
    }

    this->exx_objs[Conv_Coulomb_Pot_K::Coulomb_Method::Center2].abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp_spencer(
        this->abfs,
        center2_settings->second.second,
        this->info.ccp_rmesh_times);
    this->exx_objs[Conv_Coulomb_Pot_K::Coulomb_Method::Center2].cv.set_orbitals(
        ucell,
        orb,
        this->lcaos,
        this->abfs,
        this->exx_objs[Conv_Coulomb_Pot_K::Coulomb_Method::Center2].abfs_ccp,
        this->info.kmesh_times,
        this->MGT,
        center2_settings->second.first);

	ModuleBase::timer::end("Exx_LRI", "init_spencer");
}

template <typename Tdata>
void Exx_LRI<Tdata>::init_spencer(const MPI_Comm& mpi_comm_in,
                                  const UnitCell& ucell,
                                  const K_Vectors& kv_in,
                                  const LCAO_Orbitals& orb,
                                  const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs_in)
{
    ModuleBase::TITLE("Exx_LRI", "init_spencer");
    ModuleBase::timer::tick("Exx_LRI", "init_spencer");

    this->mpi_comm = mpi_comm_in;
    this->p_kv = &kv_in;
    this->orb_cutoff_ = orb.cutoffs();

    this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs(orb, this->info.kmesh_times);
    Exx_Abfs::Construct_Orbs::filter_empty_orbs(this->lcaos);

    this->abfs = abfs_in;
    Exx_Abfs::Construct_Orbs::filter_empty_orbs(this->abfs);
    if (this->info.rotate_abfs && ExxLriDetail::rotate_abfs_in_place_for_current_full_matrix(this->info))
    {
        ExxLriDetail::rotate_abfs_by_multipole(this->abfs, this->info.multip_moments_threshold);
    }
    Exx_Abfs::Construct_Orbs::print_orbs_size(ucell, this->abfs, GlobalV::ofs_running);

    for (size_t T = 0; T != this->abfs.size(); ++T)
    {
        GlobalC::exx_info.info_ri.abfs_Lmax
            = std::max(GlobalC::exx_info.info_ri.abfs_Lmax, static_cast<int>(this->abfs[T].size()) - 1);
    }

	    this->exx_objs.clear();
	    this->exx_objs_long.clear();
	    this->abfs_long_n0.clear();
	    this->use_rotated_n0_long_range = false;
    this->coulomb_settings.clear();
    this->coulomb_settings[Conv_Coulomb_Pot_K::Coulomb_Method::Center2]
        = std::make_pair(true,
                         ExxLriDetail::build_center2_cut_coulomb_param(
                             this->info.coulomb_param, ucell, kv_in));

    this->MGT = std::make_shared<ORB_gaunt_table>();
    const auto center2_settings = this->coulomb_settings.find(Conv_Coulomb_Pot_K::Coulomb_Method::Center2);
    if (center2_settings == this->coulomb_settings.end())
    {
        throw std::invalid_argument("Exx_LRI::init_spencer failed to prepare Center2 settings.");
    }
    this->exx_objs[Conv_Coulomb_Pot_K::Coulomb_Method::Center2].abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp_spencer(
        this->abfs,
        center2_settings->second.second,
        this->info.ccp_rmesh_times);
    this->exx_objs[Conv_Coulomb_Pot_K::Coulomb_Method::Center2].cv.set_orbitals(
        ucell,
        orb,
        this->lcaos,
        this->abfs,
        this->exx_objs[Conv_Coulomb_Pot_K::Coulomb_Method::Center2].abfs_ccp,
        this->info.kmesh_times,
        this->MGT,
        center2_settings->second.first);

    ModuleBase::timer::tick("Exx_LRI", "init_spencer");
}

template<typename Tdata>
void Exx_LRI<Tdata>::cal_exx_ions(const UnitCell& ucell,
								  const bool write_cv)
{
	ModuleBase::TITLE("Exx_LRI","cal_exx_ions");
	ModuleBase::timer::start("Exx_LRI", "cal_exx_ions");

	// init_radial_table_ions( cal_atom_centres_core(atom_pairs_core_origin), atom_pairs_core_origin );

	// this->m_abfsabfs.init_radial_table(Rradial);
	// this->m_abfslcaos_lcaos.init_radial_table(Rradial);

	std::vector<TA> atoms(ucell.nat);
	for(int iat=0; iat<ucell.nat; ++iat)
		{ atoms[iat] = iat; }
    std::set<TA> all_atoms;
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        all_atoms.insert(iat);
    }
    const bool debug_parallel_exx = ExxLriDetail::debug_parallel_exx_enabled();
    int mpi_size = 1;
    MPI_Comm_size(this->mpi_comm, &mpi_size);
	std::map<TA,TatomR> atoms_pos;
	for(int iat=0; iat<ucell.nat; ++iat)
		{ atoms_pos[iat] = RI_Util::Vector3_to_array3( ucell.atoms[ ucell.iat2it[iat] ].tau[ ucell.iat2ia[iat] ] ); }
	const std::array<TatomR,Ndim> latvec
		= {RI_Util::Vector3_to_array3(ucell.a1),
		   RI_Util::Vector3_to_array3(ucell.a2),
		   RI_Util::Vector3_to_array3(ucell.a3)};
	const std::array<Tcell,Ndim> period = {this->p_kv->nmp[0], this->p_kv->nmp[1], this->p_kv->nmp[2]};

	this->exx_lri.set_parallel(this->mpi_comm, atoms_pos, latvec, period);

	// std::max(3) for gamma_only, list_A2 should contain cell {-1,0,1}. In the future distribute will be neighbour.
	const std::array<Tcell,Ndim> period_Vs = LRI_CV_Tools::cal_latvec_range<Tcell>(1+this->info.ccp_rmesh_times, ucell, orb_cutoff_);
	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>>
		list_As_Vs = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Vs, 2, false);

	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>> Vs;
	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>> Vs_long;
	std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>> dVs;
	for (const auto& settings_list : this->coulomb_settings)
	{
		auto Vs_temp = this->exx_objs[settings_list.first].cv.cal_Vs(
			ucell,
			list_As_Vs.first,
			list_As_Vs.second[0],
			{{"writable_Vws",true}});
		this->exx_objs[settings_list.first].cv.Vws = LRI_CV_Tools::get_CVws(ucell, Vs_temp);

		if (debug_parallel_exx && settings_list.first == Conv_Coulomb_Pot_K::Coulomb_Method::Ewald)
		{
			const auto bare_vs_local_stats = ExxLriDetail::tensor_map_stats(Vs_temp);
			const auto bare_vs_global_stats
				= ExxLriDetail::reduce_tensor_map_stats(this->mpi_comm, bare_vs_local_stats);
			if (GlobalV::MY_RANK == 0)
			{
				std::cout << "EXX debug Vs bare local(rank0): "
						  << ExxLriDetail::format_tensor_map_stats(bare_vs_local_stats) << std::endl;
				std::cout << "EXX debug Vs bare reduced     : "
						  << ExxLriDetail::format_tensor_map_stats(bare_vs_global_stats) << std::endl;
			}
		}

		if (settings_list.first == Conv_Coulomb_Pot_K::Coulomb_Method::Ewald)
		{
			ExxLriDetail::overwrite_ewald_far_field_with_moment(
				this->info,
				ucell,
				this->abfs,
				list_As_Vs,
				this->exx_objs[settings_list.first].cv,
				Vs_temp);
			if (this->info.coul_moment && GlobalV::MY_RANK == 0)
			{
				std::cout << "Overwrite Ewald far-field bare Coulomb blocks with moment tensors in the current ABFS basis."
						  << std::endl;
			}

			this->exx_objs[settings_list.first].evq.init_ions(ucell, period_Vs);
			if (this->use_rotated_n0_long_range)
			{
				this->exx_objs_long[settings_list.first].evq.init_ions(ucell, period_Vs);
			}

			std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_ewald;
			std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_ewald_long;
			for (const auto& param_list : settings_list.second.second)
			{
				std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_ewald_temp;
				switch (param_list.first)
				{
					case Conv_Coulomb_Pot_K::Coulomb_Type::Fock:
					{
						const double chi
							= this->exx_objs[settings_list.first].evq.get_singular_chi(ucell, param_list.second, 2.0);
						if (this->use_rotated_n0_long_range)
						{
							std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_ewald_long_temp;
							if (mpi_size > 1)
							{
								MPI_Barrier(this->mpi_comm);
								auto Vs_bare_root
									= RI_2D_Comm::comm_map2_first(this->mpi_comm, Vs_temp, all_atoms, all_atoms);
								MPI_Barrier(this->mpi_comm);
								if (GlobalV::MY_RANK == 0)
								{
									Vs_ewald_temp
										= this->exx_objs[settings_list.first].evq.cal_short_range_Vs_serial_full(
											ucell,
											Vs_bare_root,
											period_Vs);
									Vs_ewald_long_temp
										= this->exx_objs_long[settings_list.first].evq.cal_long_range_Vs_gauss_serial_full(
											ucell,
											chi,
											period_Vs);
								}
							}
							else
							{
								Vs_ewald_temp = this->exx_objs[settings_list.first].evq.cal_short_range_Vs(
									ucell,
									list_As_Vs.first,
									list_As_Vs.second[0],
									Vs_temp);
								Vs_ewald_long_temp
									= this->exx_objs_long[settings_list.first].evq.cal_long_range_Vs_gauss(
										ucell,
										chi);
							}
							Vs_ewald_long = Vs_ewald_long.empty() ? std::move(Vs_ewald_long_temp)
																  : LRI_CV_Tools::add(Vs_ewald_long, Vs_ewald_long_temp);
						}
						else
						{
							if (mpi_size > 1)
							{
								MPI_Barrier(this->mpi_comm);
								auto Vs_bare_root
									= RI_2D_Comm::comm_map2_first(this->mpi_comm, Vs_temp, all_atoms, all_atoms);
								MPI_Barrier(this->mpi_comm);
								if (GlobalV::MY_RANK == 0)
								{
									Vs_ewald_temp = this->exx_objs[settings_list.first].evq.cal_Vs_serial_full(
										ucell,
										chi,
										Vs_bare_root,
										period_Vs);
								}
							}
							else
							{
								Vs_ewald_temp = this->exx_objs[settings_list.first].evq.cal_Vs(ucell, chi, Vs_temp);
							}
						}
						break;
					}
					default:
						throw std::invalid_argument(std::string(__FILE__) + " line " + std::to_string(__LINE__));
				}

				Vs_ewald = Vs_ewald.empty() ? std::move(Vs_ewald_temp) : LRI_CV_Tools::add(Vs_ewald, Vs_ewald_temp);
			}

			Vs_temp = std::move(Vs_ewald);
			if (debug_parallel_exx)
			{
				const auto ewald_vs_local_stats = ExxLriDetail::tensor_map_stats(Vs_temp);
				const auto ewald_vs_global_stats
					= ExxLriDetail::reduce_tensor_map_stats(this->mpi_comm, ewald_vs_local_stats);
				if (GlobalV::MY_RANK == 0)
				{
					std::cout << "EXX debug Vs ewald local(rank0): "
							  << ExxLriDetail::format_tensor_map_stats(ewald_vs_local_stats) << std::endl;
					std::cout << "EXX debug Vs ewald reduced     : "
							  << ExxLriDetail::format_tensor_map_stats(ewald_vs_global_stats) << std::endl;
				}
			}
			if (this->use_rotated_n0_long_range)
			{
				Vs_long = Vs_long.empty() ? std::move(Vs_ewald_long) : LRI_CV_Tools::add(Vs_long, Vs_ewald_long);
			}
		}

		Vs = Vs.empty() ? std::move(Vs_temp) : LRI_CV_Tools::add(Vs, Vs_temp);

		if (PARAM.inp.cal_force || PARAM.inp.cal_stress)
		{
			auto dVs_temp = this->exx_objs[settings_list.first].cv.cal_dVs(
				ucell,
				list_As_Vs.first,
				list_As_Vs.second[0],
				{{"writable_dVws",true}});
			this->exx_objs[settings_list.first].cv.dVws = LRI_CV_Tools::get_dCVws(ucell, dVs_temp);
			dVs = dVs.empty() ? std::move(dVs_temp) : LRI_CV_Tools::add(dVs, dVs_temp);
		}
	}

	if (write_cv && GlobalV::MY_RANK == 0)
	{
		LRI_CV_Tools::write_Vs_abf(Vs, PARAM.globalv.global_out_dir + "Vs");
		if (this->use_rotated_n0_long_range)
		{
			LRI_CV_Tools::write_Vs_abf(Vs_long, PARAM.globalv.global_out_dir + "Vs_long_n0");
		}
	}
	if (debug_parallel_exx)
	{
		const auto vs_local_stats = ExxLriDetail::tensor_map_stats(Vs);
		const auto vs_global_stats = ExxLriDetail::reduce_tensor_map_stats(this->mpi_comm, vs_local_stats);
		if (GlobalV::MY_RANK == 0)
		{
			std::cout << "EXX debug Vs local(rank0): "
					  << ExxLriDetail::format_tensor_map_stats(vs_local_stats) << std::endl;
			std::cout << "EXX debug Vs reduced     : "
					  << ExxLriDetail::format_tensor_map_stats(vs_global_stats) << std::endl;
		}
	}
	if (mpi_size > 1
		&& this->coulomb_settings.find(Conv_Coulomb_Pot_K::Coulomb_Method::Ewald) != this->coulomb_settings.end())
	{
		MPI_Barrier(this->mpi_comm);
		auto Vs_root = RI_2D_Comm::comm_map2_first(this->mpi_comm, Vs, all_atoms, all_atoms);
		MPI_Barrier(this->mpi_comm);
		if (debug_parallel_exx && GlobalV::MY_RANK == 0)
		{
			const auto vs_root_stats = ExxLriDetail::tensor_map_stats(Vs_root);
			std::cout << "EXX debug Vs_root full    : "
					  << ExxLriDetail::format_tensor_map_stats(vs_root_stats) << std::endl;
		}
		if (GlobalV::MY_RANK != 0)
		{
			Vs_root.clear();
		}
		Vs = std::move(Vs_root);
		if (this->use_rotated_n0_long_range)
		{
			MPI_Barrier(this->mpi_comm);
			auto Vs_long_root = RI_2D_Comm::comm_map2_first(this->mpi_comm, Vs_long, all_atoms, all_atoms);
			MPI_Barrier(this->mpi_comm);
			if (GlobalV::MY_RANK != 0)
			{
				Vs_long_root.clear();
			}
			Vs_long = std::move(Vs_long_root);
		}
	}
	const double V_threshold_short = this->info.V_threshold;
	this->exx_lri.set_Vs(std::move(Vs), V_threshold_short, this->use_rotated_n0_long_range ? "short" : "");
	if (this->use_rotated_n0_long_range)
	{
		this->exx_lri.set_Vs(std::move(Vs_long), this->info.V_threshold_long, "long");
	}

    if (cal_dCV)
	{
		std::array<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>, Ndim>
			dVs_order = LRI_CV_Tools::change_order(std::move(dVs));
		this->exx_lri.set_dVs(std::move(dVs_order), this->info.V_grad_threshold);
		if(PARAM.inp.cal_stress)
		{
			std::array<std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>,3> dVRs = LRI_CV_Tools::cal_dMRs(ucell,dVs_order);
			this->exx_lri.set_dVRs(std::move(dVRs), this->info.V_grad_R_threshold);
		}
	}

	const std::array<Tcell,Ndim> period_Cs = LRI_CV_Tools::cal_latvec_range<Tcell>(2, ucell,orb_cutoff_);
	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>>
		list_As_Cs = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Cs, 2, false);

	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>> Cs;
	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>> Cs_long;
	std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, 3>>> dCs;
	for (const auto& settings_list : this->coulomb_settings)
	{
		if (settings_list.second.first)
		{
			auto Cs_dCs = this->exx_objs[settings_list.first].cv.cal_Cs_dCs(
				ucell,
				list_As_Cs.first,
				list_As_Cs.second[0],
				{{"cal_dC",PARAM.inp.cal_force||PARAM.inp.cal_stress},
				 {"writable_Cws",true},
				 {"writable_dCws",true},
				 {"writable_Vws",false},
				 {"writable_dVws",false}});
			auto& Cs_temp = std::get<0>(Cs_dCs);
			this->exx_objs[settings_list.first].cv.Cws = LRI_CV_Tools::get_CVws(ucell, Cs_temp);
			Cs = Cs.empty() ? Cs_temp : LRI_CV_Tools::add(Cs, Cs_temp);

			if (PARAM.inp.cal_force || PARAM.inp.cal_stress)
			{
				auto& dCs_temp = std::get<1>(Cs_dCs);
				this->exx_objs[settings_list.first].cv.dCws = LRI_CV_Tools::get_dCVws(ucell, dCs_temp);
				dCs = dCs.empty() ? dCs_temp : LRI_CV_Tools::add(dCs, dCs_temp);
			}
		}
	}
	if (this->use_rotated_n0_long_range)
	{
		const auto long_ewald_it = this->exx_objs_long.find(Conv_Coulomb_Pot_K::Coulomb_Method::Ewald);
		if (long_ewald_it == this->exx_objs_long.end())
		{
			throw std::runtime_error("Missing reduced-basis Ewald object for rotated-ABFS long-range channel.");
		}
		auto Cs_long_pair = long_ewald_it->second.cv.cal_Cs_dCs(
			ucell,
			list_As_Cs.first,
			list_As_Cs.second[0],
			{{"cal_dC", false},
			 {"writable_Cws", true},
			 {"writable_dCws", false},
			 {"writable_Vws", false},
			 {"writable_dVws", false}});
		Cs_long = std::get<0>(Cs_long_pair);
		long_ewald_it->second.cv.Cws = LRI_CV_Tools::get_CVws(ucell, Cs_long);
	}
	if (write_cv && GlobalV::MY_RANK == 0)
	{
		LRI_CV_Tools::write_Cs_ao(Cs, PARAM.globalv.global_out_dir + "Cs");
	}
	if (debug_parallel_exx)
	{
		const auto cs_local_stats = ExxLriDetail::tensor_map_stats(Cs);
		const auto cs_global_stats = ExxLriDetail::reduce_tensor_map_stats(this->mpi_comm, cs_local_stats);
		if (GlobalV::MY_RANK == 0)
		{
			std::cout << "EXX debug Cs local(rank0): "
					  << ExxLriDetail::format_tensor_map_stats(cs_local_stats) << std::endl;
			std::cout << "EXX debug Cs reduced     : "
					  << ExxLriDetail::format_tensor_map_stats(cs_global_stats) << std::endl;
		}
	}
	if (this->use_rotated_n0_long_range)
	{
		this->exx_lri.set_Cs(std::move(Cs_long), this->info.C_threshold, "long");
	}
	this->exx_lri.set_Cs(std::move(Cs), this->info.C_threshold, this->use_rotated_n0_long_range ? "short" : "");

    if (cal_dCV)
	{
		std::array<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>, Ndim>
			dCs_order = LRI_CV_Tools::change_order(std::move(dCs));
		this->exx_lri.set_dCs(std::move(dCs_order), this->info.C_grad_threshold);
		if(PARAM.inp.cal_stress)
		{
			std::array<std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>,3> dCRs = LRI_CV_Tools::cal_dMRs(ucell,dCs_order);
			this->exx_lri.set_dCRs(std::move(dCRs), this->info.C_grad_R_threshold);
		}
	}
	ModuleBase::timer::end("Exx_LRI", "cal_exx_ions");
}

	#if 0
	template <typename Tdata>
	void Exx_LRI<Tdata>::cal_cut_coulomb_cs(
		std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_cut,
		std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs,
		const UnitCell& ucell,
		const bool write_cv)
{
	ModuleBase::TITLE("Exx_LRI", "cal_cut_coulomb_cs");
	ModuleBase::timer::start("Exx_LRI", "cal_cut_coulomb_cs");

	std::vector<TA> atoms(ucell.nat);
	for(int iat=0; iat<ucell.nat; ++iat)
		atoms[iat] = iat;
	std::map<TA,TatomR> atoms_pos;
	for(int iat=0; iat<ucell.nat; ++iat)
		atoms_pos[iat] = RI_Util::Vector3_to_array3( ucell.atoms[ucell.iat2it[iat]].tau[ucell.iat2ia[iat]] );
	const std::array<TatomR,Ndim> latvec
		= {RI_Util::Vector3_to_array3(ucell.a1),
		   RI_Util::Vector3_to_array3(ucell.a2),
		   RI_Util::Vector3_to_array3(ucell.a3)};
	const std::array<Tcell,Ndim> period = {this->p_kv->nmp[0], this->p_kv->nmp[1], this->p_kv->nmp[2]};

	this->exx_lri.set_parallel(this->mpi_comm, atoms_pos, latvec, period);

	// std::max(3) for gamma_only, list_A2 should contain cell {-1,0,1}. In the future distribute will be neighbour.
	const std::array<Tcell, Ndim> period_Vs
		= LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->info.ccp_rmesh_times, ucell, orb_cutoff_);
	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>>
		list_As_Vs = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Vs, 2, false);

	std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>> dVs;
	for(const auto &settings_list : this->coulomb_settings)
	{
		if(!settings_list.second.first) continue;
		std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>
			Vs_temp = this->exx_objs[settings_list.first].cv.cal_Vs(ucell,
				list_As_Vs.first, list_As_Vs.second[0],
				{{"writable_Vws",true}});
		this->exx_objs[settings_list.first].cv.Vws = LRI_CV_Tools::get_CVws(ucell,Vs_temp);

		// ======rotate ABFs begin======
		int flag = 0;
		for (const auto& IJRc: this->exx_objs[settings_list.first].cv.Vws)
		{
			const TA& I = IJRc.first;
			const auto& JRc = IJRc.second;
			for (const auto& JRc_tensor: JRc)
			{
				const TA& J = JRc_tensor.first;
				const auto Rc = JRc_tensor.second;
				for (const auto& Rc_tensor: Rc)
				{
					const auto& R = Rc_tensor.first;
					flag += 1;
				}
			}
		}
		std::cout << "Coulomb: number of atom-pairs inside atomic overlap is " << flag << ". " << std::endl;
		if (this->info.coul_moment == true)
		{
			double hf_Rcut = std::pow(0.75 * this->p_kv->get_nkstot_nospin() * ucell.omega / (ModuleBase::PI), 1.0 / 3.0);
			// To cal Cs, we still cal all Vs(R) in r space
			// moment_abfs->cal_VR(ucell,
			//					 this->abfs,
			//					 list_As_Vs,
			//					 orb_cutoff_,
			//					 hf_Rcut,
			//					 this->exx_objs[settings_list.first].cv,
			//					 Vs_cut);
			delete moment_abfs;
			moment_abfs = nullptr;
			malloc_trim(0);
		}

		flag = 0;
		for (const auto& IJRc: this->exx_objs[settings_list.first].cv.Vws)
		{
			const auto& JRc = IJRc.second;
			for (const auto& JRc_tensor: JRc)
			{
				const auto Rc = JRc_tensor.second;
				for (const auto& Rc_tensor: Rc)
				{
					flag += 1;
				}
			}
		}
		std::cout << "Coulomb: number of all atom-pairs is " << flag << ". " << std::endl;
		// ======rotate ABFs end======
		
		Vs_cut = Vs.empty() ? Vs_temp : LRI_CV_Tools::add(Vs_cut, Vs_temp);

		if(PARAM.inp.cal_force || PARAM.inp.cal_stress)
		{
			std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>>
				dVs_temp = this->exx_objs[settings_list.first].cv.cal_dVs(ucell,
					list_As_Vs.first, list_As_Vs.second[0],
					{{"writable_dVws",true}});
			this->exx_objs[settings_list.first].cv.dVws = LRI_CV_Tools::get_dCVws(ucell,dVs_temp);
			dVs = dVs.empty() ? dVs_temp : LRI_CV_Tools::add(dVs, dVs_temp);
		}
	}

	if (write_cv && GlobalV::MY_RANK == 0)
	{
		LRI_CV_Tools::write_Vs_abf(Vs_cut, PARAM.globalv.global_out_dir + "Vs_cut");
	}
	this->exx_lri.set_Vs(std::move(Vs_cut), this->info.V_threshold);

	if(PARAM.inp.cal_force || PARAM.inp.cal_stress)
	{
		std::map<TA,std::map<TAC,std::array<RI::Tensor<Tdata>,Ndim>>> dVs
			= this->exx_objs[coulomb_method].cv.cal_dVs(ucell,
				list_As_Vs.first, list_As_Vs.second[0],
				{{"writable_dVws",true}});
		this->exx_objs[coulomb_method].cv.dVws = LRI_CV_Tools::get_dCVws(ucell,dVs);

		std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,Ndim> dVs_order
			= LRI_CV_Tools::change_order(std::move(dVs));
		this->exx_lri.set_dVs(std::move(dVs_order), this->info.V_grad_threshold);
		if(PARAM.inp.cal_stress)
		{
			std::array<std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>,3> dVRs
				= LRI_CV_Tools::cal_dMRs(ucell,dVs_order);
			this->exx_lri.set_dVRs(std::move(dVRs), this->info.V_grad_R_threshold);
		}
	}

	const std::array<Tcell,Ndim> period_Cs = LRI_CV_Tools::cal_latvec_range<Tcell>(2, ucell,orb_cutoff_);
	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>>
		list_As_Cs = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Cs, 2, false);
	std::pair<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,
			  std::map<TA,std::map<TAC,std::array<RI::Tensor<Tdata>,3>>>>
		Cs_dCs = this->exx_objs[coulomb_method].cv.cal_Cs_dCs(ucell,
			list_As_Cs.first, list_As_Cs.second[0],
			{{"cal_dC",PARAM.inp.cal_force||PARAM.inp.cal_stress},
			 {"writable_Cws",true},
			 {"writable_dCws",true},
			 {"writable_Vws",false},
			 {"writable_dVws",false}});
	Cs = std::get<0>(Cs_dCs);
	this->exx_objs[coulomb_method].cv.Cws = LRI_CV_Tools::get_CVws(ucell,Cs);
	if(write_cv && GlobalV::MY_RANK==0)
	{
		LRI_CV_Tools::write_Cs_ao(Cs, PARAM.globalv.global_out_dir + "Cs");
	}
	this->exx_lri.set_Cs(Cs, this->info.C_threshold);

	if(PARAM.inp.cal_force || PARAM.inp.cal_stress)
	{
		std::map<TA,std::map<TAC,std::array<RI::Tensor<Tdata>,3>>>& dCs = std::get<1>(Cs_dCs);
		this->exx_objs[coulomb_method].cv.dCws = LRI_CV_Tools::get_dCVws(ucell,dCs);
		std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,Ndim> dCs_order
			= LRI_CV_Tools::change_order(std::move(dCs));
		this->exx_lri.set_dCs(std::move(dCs_order), this->info.C_grad_threshold);
		if(PARAM.inp.cal_stress)
		{
			std::array<std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>,3> dCRs
				= LRI_CV_Tools::cal_dMRs(ucell,dCs_order);
			this->exx_lri.set_dCRs(std::move(dCRs), this->info.C_grad_R_threshold);
		}
	}
	ModuleBase::timer::end("Exx_LRI", "cal_cut_coulomb_cs");
}

template <typename Tdata>
void Exx_LRI<Tdata>::cal_ewald_coulomb(
	std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_full,
	std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs,
	const UnitCell& ucell,
	const bool write_cv)
{
	ModuleBase::TITLE("Exx_LRI", "cal_ewald_coulomb");
	ModuleBase::timer::start("Exx_LRI", "cal_ewald_coulomb");

	std::vector<TA> atoms(ucell.nat);
	for (int iat = 0; iat < ucell.nat; ++iat)
		atoms[iat] = iat;
	std::map<TA, TatomR> atoms_pos;
	for (int iat = 0; iat < ucell.nat; ++iat)
		atoms_pos[iat] = RI_Util::Vector3_to_array3(ucell.atoms[ucell.iat2it[iat]].tau[ucell.iat2ia[iat]]);
	const std::array<TatomR, Ndim> latvec = {
		RI_Util::Vector3_to_array3(ucell.a1),
		RI_Util::Vector3_to_array3(ucell.a2),
		RI_Util::Vector3_to_array3(ucell.a3)};
	const std::array<Tcell, Ndim> period = {this->p_kv->nmp[0], this->p_kv->nmp[1], this->p_kv->nmp[2]};

	this->exx_lri.set_parallel(this->mpi_comm, atoms_pos, latvec, period);

	// std::max(3) for gamma_only, list_A2 should contain cell {-1,0,1}. In the future distribute will be neighbour.
	const std::array<Tcell, Ndim> period_Vs
		= LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->info.ccp_rmesh_times, ucell, orb_cutoff_);
	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>> list_As_Vs
		= RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Vs, 2, false);

	std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, Ndim>>> dVs;
	for (const auto& settings_list: this->coulomb_settings)
	{
		if (!settings_list.second.first)
			continue;
		std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_temp
			= this->exx_objs[settings_list.first].cv.cal_Vs(
				ucell,
				list_As_Vs.first,
				list_As_Vs.second[0],
				{{"writable_Vws", true}});
		this->exx_objs[settings_list.first].cv.Vws = LRI_CV_Tools::get_CVws(ucell, Vs_temp);

		// ======rotate ABFs begin======
		int flag = 0;
		for (const auto& IJRc: this->exx_objs[settings_list.first].cv.Vws)
		{
			const TA& I = IJRc.first;
			const auto& JRc = IJRc.second;
			for (const auto& JRc_tensor: JRc)
			{
				const TA& J = JRc_tensor.first;
				const auto Rc = JRc_tensor.second;
				for (const auto& Rc_tensor: Rc)
				{
					const auto& R = Rc_tensor.first;
					flag += 1;
				}
			}
		}
		std::cout << "Coulomb: number of atom-pairs inside atomic overlap is " << flag << ". " << std::endl;
		if (this->info.coul_moment == true)
		{
			double hf_Rcut = std::pow(0.75 * this->p_kv->get_nkstot_nospin() * ucell.omega / (ModuleBase::PI), 1.0 / 3.0);
			// To cal Cs, we still cal all Vs(R) in r space
			// moment_abfs->cal_VR(ucell,
			//					 this->abfs,
			//					 list_As_Vs,
			//					 orb_cutoff_,
			//					 hf_Rcut,
			//					 this->exx_objs[settings_list.first].cv,
			//					 Vs_full);
			delete moment_abfs;
			moment_abfs = nullptr;
			malloc_trim(0);
		}

		flag = 0;
		for (const auto& IJRc: this->exx_objs[settings_list.first].cv.Vws)
		{
			const auto& JRc = IJRc.second;
			for (const auto& JRc_tensor: JRc)
			{
				const auto Rc = JRc_tensor.second;
				for (const auto& Rc_tensor: Rc)
				{
					flag += 1;
				}
			}
		}
		std::cout << "Coulomb: number of all atom-pairs is " << flag << ". " << std::endl;
		// ======rotate ABFs end======

		if (settings_list.first == Conv_Coulomb_Pot_K::Coulomb_Method::Ewald)
		{
			this->exx_objs[settings_list.first].evq.init_ions(ucell, period_Vs);
			const auto& coulomb_param = settings_list.second.second;
			std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_ewald;
			for (const auto& param_list: coulomb_param)
			{
				std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_ewald_temp;
				switch (param_list.first)
				{
				case Conv_Coulomb_Pot_K::Coulomb_Type::Fock: {
					double chi
						= this->exx_objs[settings_list.first].evq.get_singular_chi(ucell, param_list.second, 2.0);
					Vs_ewald_temp = this->exx_objs[settings_list.first].evq.cal_Vs(ucell, chi, Vs_temp);
					break;
				}
				default: {
					throw std::invalid_argument(std::string(__FILE__) + " line " + std::to_string(__LINE__));
				}
				}
				// Vs_temp cannot be covered here
				Vs_ewald = Vs_ewald.empty() ? Vs_ewald_temp : LRI_CV_Tools::add(Vs_ewald, Vs_ewald_temp);
			}
			Vs_temp = Vs_ewald;
		}

		Vs_full = Vs.empty() ? Vs_temp : LRI_CV_Tools::add(Vs_full, Vs_temp);
	}

	if (write_cv && GlobalV::MY_RANK == 0)
	{
		LRI_CV_Tools::write_Vs_abf(Vs_full, PARAM.globalv.global_out_dir + "Vs_full");
	}
	// this->exx_lri.set_Vs(std::move(Vs_full), this->info.V_threshold);

	// const std::array<Tcell,Ndim> period_Cs = LRI_CV_Tools::cal_latvec_range<Tcell>(2, ucell,orb_cutoff_);
	// const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>>
	// 	list_As_Cs = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Cs, 2, false);
	// std::pair<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,
	// 		  std::map<TA,std::map<TAC,std::array<RI::Tensor<Tdata>,3>>>>
	// 	Cs_dCs = this->exx_objs[coulomb_method].cv.cal_Cs_dCs(ucell,
	// 		list_As_Cs.first, list_As_Cs.second[0],
	// 		{{"cal_dC",PARAM.inp.cal_force||PARAM.inp.cal_stress},
	// 		 {"writable_Cws",true},
	// 		 {"writable_dCws",true},
	// 		 {"writable_Vws",false},
	// 		 {"writable_dVws",false}});
	// Cs = std::get<0>(Cs_dCs);
	// this->exx_objs[coulomb_method].cv.Cws = LRI_CV_Tools::get_CVws(ucell,Cs);
	// if(write_cv && GlobalV::MY_RANK==0)
	// {
	// 	LRI_CV_Tools::write_Cs_ao(Cs, PARAM.globalv.global_out_dir + "Cs");
	// }
	// this->exx_lri.set_Cs(Cs, this->info.C_threshold);
	ModuleBase::timer::end("Exx_LRI", "cal_ewald_coulomb");
}
	#endif

template <typename Tdata>
void Exx_LRI<Tdata>::cal_cut_coulomb_cs(
	std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_cut,
	std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs,
	const UnitCell& ucell,
	const bool write_cv)
{
	ModuleBase::TITLE("Exx_LRI", "cal_cut_coulomb_cs");
	ModuleBase::timer::start("Exx_LRI", "cal_cut_coulomb_cs");

	std::vector<TA> atoms(ucell.nat);
	for (int iat = 0; iat < ucell.nat; ++iat)
	{
		atoms[iat] = iat;
	}
	std::map<TA, TatomR> atoms_pos;
	for (int iat = 0; iat < ucell.nat; ++iat)
	{
		atoms_pos[iat] = RI_Util::Vector3_to_array3(ucell.atoms[ucell.iat2it[iat]].tau[ucell.iat2ia[iat]]);
	}
	const std::array<TatomR, Ndim> latvec = {
		RI_Util::Vector3_to_array3(ucell.a1),
		RI_Util::Vector3_to_array3(ucell.a2),
		RI_Util::Vector3_to_array3(ucell.a3)};
	const std::array<Tcell, Ndim> period = {this->p_kv->nmp[0], this->p_kv->nmp[1], this->p_kv->nmp[2]};

	this->exx_lri.set_parallel(this->mpi_comm, atoms_pos, latvec, period);

	const auto center2_method = Conv_Coulomb_Pot_K::Coulomb_Method::Center2;
	auto center2_obj_it = this->exx_objs.find(center2_method);
	if (center2_obj_it == this->exx_objs.end())
	{
		throw std::invalid_argument(std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}

	const std::array<Tcell, Ndim> period_Vs
		= LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->info.ccp_rmesh_times, ucell, orb_cutoff_);
	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>> list_As_Vs
		= RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Vs, 2, false);

	Vs_cut = center2_obj_it->second.cv.cal_Vs(
		ucell,
		list_As_Vs.first,
		list_As_Vs.second[0],
		{{"writable_Vws", true}});
	center2_obj_it->second.cv.Vws = LRI_CV_Tools::get_CVws(ucell, Vs_cut);
	if (write_cv && GlobalV::MY_RANK == 0)
	{
		LRI_CV_Tools::write_Vs_abf(Vs_cut, PARAM.globalv.global_out_dir + "Vs_cut");
	}
	this->exx_lri.set_Vs(Vs_cut, this->info.V_threshold);

	const std::array<Tcell, Ndim> period_Cs = LRI_CV_Tools::cal_latvec_range<Tcell>(2, ucell, orb_cutoff_);
	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>> list_As_Cs
		= RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Cs, 2, false);
	std::pair<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>,
			  std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, 3>>>>
		Cs_dCs = center2_obj_it->second.cv.cal_Cs_dCs(
			ucell,
			list_As_Cs.first,
			list_As_Cs.second[0],
			{{"cal_dC", false},
			 {"writable_Cws", true},
			 {"writable_dCws", true},
			 {"writable_Vws", false},
			 {"writable_dVws", false}});
	Cs = std::get<0>(Cs_dCs);
	center2_obj_it->second.cv.Cws = LRI_CV_Tools::get_CVws(ucell, Cs);
	if (write_cv && GlobalV::MY_RANK == 0)
	{
		LRI_CV_Tools::write_Cs_ao(Cs, PARAM.globalv.global_out_dir + "Cs");
	}
	this->exx_lri.set_Cs(Cs, this->info.C_threshold);

	ModuleBase::timer::end("Exx_LRI", "cal_cut_coulomb_cs");
}

template <typename Tdata>
void Exx_LRI<Tdata>::cal_ewald_coulomb(
	std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs_full,
	std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs,
	const UnitCell& ucell,
	const bool write_cv)
{
	ModuleBase::TITLE("Exx_LRI", "cal_ewald_coulomb");
	ModuleBase::timer::start("Exx_LRI", "cal_ewald_coulomb");

	std::vector<TA> atoms(ucell.nat);
	for (int iat = 0; iat < ucell.nat; ++iat)
	{
		atoms[iat] = iat;
	}
	std::map<TA, TatomR> atoms_pos;
	for (int iat = 0; iat < ucell.nat; ++iat)
	{
		atoms_pos[iat] = RI_Util::Vector3_to_array3(ucell.atoms[ucell.iat2it[iat]].tau[ucell.iat2ia[iat]]);
	}
	const std::array<TatomR, Ndim> latvec = {
		RI_Util::Vector3_to_array3(ucell.a1),
		RI_Util::Vector3_to_array3(ucell.a2),
		RI_Util::Vector3_to_array3(ucell.a3)};
	const std::array<Tcell, Ndim> period = {this->p_kv->nmp[0], this->p_kv->nmp[1], this->p_kv->nmp[2]};

	this->exx_lri.set_parallel(this->mpi_comm, atoms_pos, latvec, period);

	const std::array<Tcell, Ndim> period_Vs
		= LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->info.ccp_rmesh_times, ucell, orb_cutoff_);
	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>> list_As_Vs
		= RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Vs, 2, false);

	for (const auto& settings_list : this->coulomb_settings)
	{
		std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_temp
			= this->exx_objs[settings_list.first].cv.cal_Vs(
				ucell,
				list_As_Vs.first,
				list_As_Vs.second[0],
				{{"writable_Vws", true}});
		this->exx_objs[settings_list.first].cv.Vws = LRI_CV_Tools::get_CVws(ucell, Vs_temp);

		if (settings_list.first == Conv_Coulomb_Pot_K::Coulomb_Method::Ewald)
		{
            ExxLriDetail::overwrite_ewald_far_field_with_moment(
                this->info,
                ucell,
                this->abfs,
                list_As_Vs,
                this->exx_objs[settings_list.first].cv,
                Vs_temp);
			this->exx_objs[settings_list.first].evq.init_ions(ucell, period_Vs);
			std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_ewald;
			for (const auto& param_list : settings_list.second.second)
			{
				std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_ewald_temp;
				switch (param_list.first)
				{
				case Conv_Coulomb_Pot_K::Coulomb_Type::Fock:
				{
					double chi = this->exx_objs[settings_list.first].evq.get_singular_chi(ucell, param_list.second, 2.0);
					Vs_ewald_temp = this->exx_objs[settings_list.first].evq.cal_Vs(ucell, chi, Vs_temp);
					break;
				}
				default:
				{
					throw std::invalid_argument(std::string(__FILE__) + " line " + std::to_string(__LINE__));
				}
				}
				Vs_ewald = Vs_ewald.empty() ? Vs_ewald_temp : LRI_CV_Tools::add(Vs_ewald, Vs_ewald_temp);
			}
			Vs_temp = Vs_ewald;
		}

		Vs_full = Vs_full.empty() ? Vs_temp : LRI_CV_Tools::add(Vs_full, Vs_temp);
	}

	if (write_cv && GlobalV::MY_RANK == 0)
	{
		LRI_CV_Tools::write_Vs_abf(Vs_full, PARAM.globalv.global_out_dir + "Vs_full");
	}

	Cs.clear();
	ModuleBase::timer::end("Exx_LRI", "cal_ewald_coulomb");
}

template<typename Tdata>
void Exx_LRI<Tdata>::cal_exx_elec(const std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>>& Ds,
	const UnitCell& ucell,
	const Parallel_Orbitals& pv,
	const ModuleSymmetry::Symmetry_rotation* p_symrot)
{
	ModuleBase::TITLE("Exx_LRI","cal_exx_elec");
	ModuleBase::timer::start("Exx_LRI", "cal_exx_elec");

	const std::vector<std::tuple<std::set<TA>, std::set<TA>>> judge = RI_2D_Comm::get_2D_judge(ucell,pv);
    std::set<TA> all_atoms;
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        all_atoms.insert(iat);
    }
    const bool debug_parallel_exx = ExxLriDetail::debug_parallel_exx_enabled();

	if(p_symrot)
		{
            this->exx_lri.set_symmetry(true, p_symrot->get_irreducible_sector());
        }
	else
		{
            this->exx_lri.set_symmetry(false, {});
        }

    double full_cal_hs_time = 0.0;
    double short_cal_hs_time = 0.0;
    double long_cal_hs_time = 0.0;
    const int cal_hs_benchmark_repeat = ExxLriDetail::get_cal_hs_benchmark_repeat();
    const bool h_only_rt_mode = (PARAM.inp.esolver_type == "tddft");

	    auto run_exx_channel =
	        [&](RI::Exx<TA, Tcell, Ndim, Tdata>& exx_channel,
	            const std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& D_in,
	            const int spin_index,
	            const std::string& ds_suffix,
	            const std::string& cv_suffix,
	            double& cal_hs_time_acc) -> std::pair<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>, double>
	    {
	        if (debug_parallel_exx)
	        {
	            const auto d_local_stats = ExxLriDetail::tensor_map_stats(D_in);
	            const auto d_full_stats = ExxLriDetail::tensor_map_stats(
	                RI_2D_Comm::comm_map2_first(this->mpi_comm, D_in, all_atoms, all_atoms));
	            if (GlobalV::MY_RANK == 0)
	            {
	                std::cout << "EXX debug D[" << ds_suffix << "] local: "
	                          << ExxLriDetail::format_tensor_map_stats(d_local_stats) << std::endl;
	                std::cout << "EXX debug D[" << ds_suffix << "] full : "
	                          << ExxLriDetail::format_tensor_map_stats(d_full_stats) << std::endl;
	            }
	        }

	        if (h_only_rt_mode)
	        {
	            exx_channel.set_Ds_no_post_2d(D_in, this->info.dm_threshold, ds_suffix);
	        }
	        else
	        {
	            exx_channel.set_Ds(D_in, this->info.dm_threshold, ds_suffix);
	        }
	        const auto cal_hs_t0 = std::chrono::steady_clock::now();
	        for (int irepeat = 0; irepeat < cal_hs_benchmark_repeat; ++irepeat)
	        {
	            if (h_only_rt_mode)
	            {
	                exx_channel.cal_Hs_only({cv_suffix, cv_suffix, ds_suffix});
	            }
	            else
	            {
	                exx_channel.cal_Hs({cv_suffix, cv_suffix, ds_suffix});
	            }
	        }
	        const auto cal_hs_t1 = std::chrono::steady_clock::now();
        cal_hs_time_acc
            += std::chrono::duration<double>(cal_hs_t1 - cal_hs_t0).count() / static_cast<double>(cal_hs_benchmark_repeat);

        if (debug_parallel_exx)
        {
            const auto h_local_stats = ExxLriDetail::tensor_map_stats(exx_channel.Hs);
	            const auto h_full_stats = ExxLriDetail::tensor_map_stats(
	                RI_2D_Comm::comm_map2_first(this->mpi_comm, exx_channel.Hs, all_atoms, all_atoms));
	            if (GlobalV::MY_RANK == 0)
	            {
	                std::cout << "EXX debug H[" << ds_suffix << "] local: "
	                          << ExxLriDetail::format_tensor_map_stats(h_local_stats) << std::endl;
	                std::cout << "EXX debug H[" << ds_suffix << "] full : "
	                          << ExxLriDetail::format_tensor_map_stats(h_full_stats)
	                          << ", energy=" << std::real(exx_channel.energy) << std::endl;
	            }
        }

        if (!p_symrot)
        {
            return std::make_pair(
                RI::Communicate_Tensors_Map_Judge::comm_map2_first(
                    this->mpi_comm,
                    std::move(exx_channel.Hs),
                    std::get<0>(judge[spin_index]),
                    std::get<1>(judge[spin_index])),
                h_only_rt_mode ? 0.0 : std::real(exx_channel.energy));
        }

	        auto Hs_a2D = exx_channel.post_2D.set_tensors_map2(exx_channel.Hs);
	        Hs_a2D = p_symrot->restore_HR(ucell.symm, ucell.atoms, ucell.st, 'H', Hs_a2D);
	        if (!h_only_rt_mode)
	        {
	            exx_channel.energy = exx_channel.post_2D.cal_energy(exx_channel.post_2D.saves["Ds_" + ds_suffix],
	                                                                exx_channel.post_2D.set_tensors_map2(Hs_a2D));
	        }
	        return std::make_pair(
	            RI::Communicate_Tensors_Map_Judge::comm_map2_first(
                this->mpi_comm,
                std::move(Hs_a2D),
                std::get<0>(judge[spin_index]),
                std::get<1>(judge[spin_index])),
            h_only_rt_mode ? 0.0 : std::real(exx_channel.energy));
    };

	// (nspin=4) the spinor H(R) rotation mixes the 4 spin channels, so they cannot be rotated
	// independently in the per-spin loop below; hand off to the SOC implementation.
	if (p_symrot && PARAM.inp.nspin == 4)
	{
		this->cal_exx_elec_soc(Ds, ucell, judge, p_symrot);
		this->exx_lri.set_symmetry(false, {});
		ModuleBase::timer::end("Exx_LRI", "cal_exx_elec");
		return;
	}

	this->Hexxs.resize(PARAM.inp.nspin);
	this->Eexx = 0;
	for(int is=0; is<PARAM.inp.nspin; ++is)
	{
		const std::string suffix = ((PARAM.inp.cal_force || PARAM.inp.cal_stress) ? std::to_string(is) : "");

	        auto short_channel = run_exx_channel(this->exx_lri,
	                                             Ds[is],
	                                             is,
	                                             suffix,
	                                             this->use_rotated_n0_long_range ? "short" : "",
	                                             this->use_rotated_n0_long_range ? short_cal_hs_time : full_cal_hs_time);
	        this->Hexxs[is] = std::move(short_channel.first);
	        this->Eexx += short_channel.second;
	        if (this->use_rotated_n0_long_range)
	        {
	            auto long_channel = run_exx_channel(this->exx_lri,
	                                                Ds[is],
	                                                is,
	                                                suffix + "_lr",
	                                                "long",
	                                                long_cal_hs_time);
            this->Hexxs[is] = LRI_CV_Tools::add(this->Hexxs[is], long_channel.first);
            this->Eexx += long_channel.second;
        }
		post_process_Hexx(this->Hexxs[is]);
	}
	this->Eexx = h_only_rt_mode ? 0.0 : post_process_Eexx(this->Eexx);
	this->exx_lri.set_symmetry(false, {});
    if (GlobalV::MY_RANK == 0)
    {
        if (this->use_rotated_n0_long_range)
        {
            std::cout << "EXX cal_Hs timing summary: short = " << short_cal_hs_time
                      << " s, long = " << long_cal_hs_time
                      << " s, repeat = " << cal_hs_benchmark_repeat << std::endl;
        }
        else
        {
            std::cout << "EXX cal_Hs timing summary: full = " << full_cal_hs_time
                      << " s, repeat = " << cal_hs_benchmark_repeat << std::endl;
        }
    }
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_elec");
}

template<typename Tdata>
void Exx_LRI<Tdata>::post_process_Hexx( std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> &Hexxs_io ) const
{
	ModuleBase::TITLE("Exx_LRI","post_process_Hexx");
	constexpr Tdata frac = -1 * 2;								// why?	Hartree to Ry?
	const std::function<void(RI::Tensor<Tdata>&)>
		multiply_frac = [&frac](RI::Tensor<Tdata> &t)
		{ t = t*frac; };
	RI::Map_Operator::for_each( Hexxs_io, multiply_frac );
}

template<typename Tdata>
double Exx_LRI<Tdata>::post_process_Eexx(const double& Eexx_in) const
{
	ModuleBase::TITLE("Exx_LRI","post_process_Eexx");
	const double SPIN_multiple = std::map<int, double>{ {1,2}, {2,1}, {4,1} }.at(PARAM.inp.nspin);				// why?
	const double frac = -SPIN_multiple;
	return frac * Eexx_in;
}

/*
post_process_old
{
	// D
	const std::map<int,double> SPIN_multiple = {{1,0.5}, {2,1}, {4,1}};							// ???
	DR *= SPIN_multiple.at(NSPIN);

	// H
	HR *= -2;

	// E
	const std::map<int,double> SPIN_multiple = {{1,2}, {2,1}, {4,1}};							// ???
	energy *= SPIN_multiple.at(PARAM.inp.nspin);			// ?
	energy /= 2;					// /2 for Ry
}
*/

template<typename Tdata>
void Exx_LRI<Tdata>::cal_exx_force(const int& nat)
{
	ModuleBase::TITLE("Exx_LRI","cal_exx_force");
	ModuleBase::timer::start("Exx_LRI", "cal_exx_force");

	this->force_exx.create(nat, Ndim);
	for(int is=0; is<PARAM.inp.nspin; ++is)
	{
		this->exx_lri.cal_force({"","",std::to_string(is),"",""});
		for(std::size_t idim=0; idim<Ndim; ++idim) {
			for(const auto &force_item : this->exx_lri.force[idim]) {
				this->force_exx(force_item.first, idim) += std::real(force_item.second);
					} 		}
	}

	const double SPIN_multiple = std::map<int,double>{{1,2}, {2,1}, {4,1}}.at(PARAM.inp.nspin);				// why?
	const double frac = -2 * SPIN_multiple;		// why?
	this->force_exx *= frac;
	ModuleBase::timer::end("Exx_LRI", "cal_exx_force");
}


template<typename Tdata>
void Exx_LRI<Tdata>::cal_exx_stress(const double& omega, const double& lat0)
{
	ModuleBase::TITLE("Exx_LRI","cal_exx_stress");
	ModuleBase::timer::start("Exx_LRI", "cal_exx_stress");

	this->stress_exx.create(Ndim, Ndim);
	for(int is=0; is<PARAM.inp.nspin; ++is)
	{
		this->exx_lri.cal_stress({"","",std::to_string(is),"",""});
		for(std::size_t idim0=0; idim0<Ndim; ++idim0) {
			for(std::size_t idim1=0; idim1<Ndim; ++idim1) {
				this->stress_exx(idim0,idim1) += std::real(this->exx_lri.stress(idim0,idim1));
				} 	}
	}

	const double SPIN_multiple = std::map<int,double>{{1,2}, {2,1}, {4,1}}.at(PARAM.inp.nspin);				// why?
	const double frac = 2 * SPIN_multiple / omega * lat0;		// why?
	this->stress_exx *= frac;

	ModuleBase::timer::end("Exx_LRI", "cal_exx_stress");
}

template<typename Tdata>
std::vector<std::vector<int>> Exx_LRI<Tdata>::get_abfs_nchis() const
{
	std::vector<std::vector<int>> abfs_nchis;
	for (const auto& abfs_T : this->abfs)
	{
		std::vector<int> abfs_nchi_T;
		for (const auto& abfs_L : abfs_T)
			{ abfs_nchi_T.push_back(abfs_L.size()); }
		abfs_nchis.push_back(abfs_nchi_T);
	}
	return abfs_nchis;
}

#endif
