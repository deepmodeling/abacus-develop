#ifndef SPARSE_FORMAT_HSR_H
#define SPARSE_FORMAT_HSR_H

#include "source_lcao/hamilt_lcao.h"

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
void cal_HSR(const UnitCell& ucell,
             const Parallel_Orbitals& pv,
             LCAO_HS_Arrays& HS_Arrays,
             const Grid_Driver& grid,
             const int& current_spin,
             const double& sparse_thr,
             const int (&nmp)[3],
             hamilt::Hamilt<std::complex<double>>* p_ham
#ifdef __EXX
             ,
             const std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd = nullptr,
             const std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc = nullptr
#endif
);

void cal_HContainer_d(const Parallel_Orbitals& pv,
                      const int& current_spin,
                      const double& sparse_threshold,
                      const hamilt::HContainer<double>& hR,
                      std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>& target);

void cal_HContainer_cd(
    const Parallel_Orbitals& pv,
    const int& current_spin,
    const double& sparse_threshold,
    const hamilt::HContainer<std::complex<double>>& hR,
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>& target);

void cal_HContainer_td(
    const Parallel_Orbitals& pv,
    const int& current_spin,
    const double& sparse_threshold,
    const hamilt::HContainer<double>& hR,
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>& target);

void clear_zero_elements(LCAO_HS_Arrays& HS_Arrays, const int& current_spin, const double& sparse_thr);

void destroy_HS_R_sparse(LCAO_HS_Arrays& HS_Arrays);

} // namespace sparse_format

#endif