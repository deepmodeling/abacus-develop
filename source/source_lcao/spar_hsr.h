#ifndef SPARSE_FORMAT_HSR_H
#define SPARSE_FORMAT_HSR_H

#include "source_lcao/LCAO_HS_arrays.hpp"
#include "source_lcao/hamilt_lcao.h"
#include "source_lcao/module_dftu/dftu.h" // mohan add 20251107

namespace sparse_format
{
#ifdef __MPI
// Synchronize all processes' R coordinates,
// otherwise HS_Arrays.all_R_coor would have different sizes on different processes,
// causing MPI communication errors.
void sync_all_R_coor(std::set<Abfs::Vector3_Order<int>>& all_R_coor, MPI_Comm comm);
#endif

template <typename T>
std::set<Abfs::Vector3_Order<int>> get_R_range(const hamilt::HContainer<T>& hR)
{
    std::set<Abfs::Vector3_Order<int>> all_R_coor;
    for (int iap = 0; iap < hR.size_atom_pairs(); ++iap)
    {
        const hamilt::AtomPair<T>& atom_pair = hR.get_atom_pair(iap);
        for (int iR = 0; iR < atom_pair.get_R_size(); ++iR)
        {
            const auto& r_index = atom_pair.get_R_index(iR);
            Abfs::Vector3_Order<int> dR(r_index.x, r_index.y, r_index.z);
            all_R_coor.insert(dR);
        }
    }

#ifdef __MPI
    // Fix: Sync all_R_coor across processes
    sparse_format::sync_all_R_coor(all_R_coor, MPI_COMM_WORLD);
#endif

    return all_R_coor;
};

using TAC = std::pair<int, std::array<int, 3>>;
template <typename TK>
void cal_HSR(const UnitCell& ucell,
		Plus_U &dftu, // mohan add 2025-11-07
		const Parallel_Orbitals& pv,
		LCAO_HS_Arrays& HS_Arrays,
		const Grid_Driver& grid,
		const int& current_spin,
		const double& sparse_thr,
		const int (&nmp)[3],
		hamilt::Hamilt<TK>* p_ham
#ifdef __EXX
		,
		const std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd = nullptr,
		const std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc = nullptr
#endif
		);

template <typename TI, typename TO = TI>
void cal_HContainer(const Parallel_Orbitals& pv,
                    const double& sparse_thr,
                    const hamilt::HContainer<TI>& hR,
                    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, TO>>>& target);

void clear_zero_elements(LCAO_HS_Arrays& HS_Arrays, const int& current_spin, const double& sparse_thr);

void destroy_HS_R_sparse(LCAO_HS_Arrays& HS_Arrays);

} // namespace sparse_format

#endif
