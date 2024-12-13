#ifndef SPARSE_FORMAT_HSR_H
#define SPARSE_FORMAT_HSR_H

#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"

namespace sparse_format {
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

void cal_HContainer_d(
    const Parallel_Orbitals& pv,
    const int& current_spin,
    const double& sparse_threshold,
    const hamilt::HContainer<double>& hR,
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, double>>>& target);

void cal_HContainer_cd(
    const Parallel_Orbitals& pv,
    const int& current_spin,
    const double& sparse_threshold,
    const hamilt::HContainer<std::complex<double>>& hR,
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, std::complex<double>>>>& target);

void cal_HContainer_td(
    const Parallel_Orbitals& pv,
    const int& current_spin,
    const double& sparse_threshold,
    const hamilt::HContainer<double>& hR,
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, std::complex<double>>>>& target);

void clear_zero_elements(LCAO_HS_Arrays& HS_Arrays,
                         const int& current_spin,
                         const double& sparse_thr);

void destroy_HS_R_sparse(LCAO_HS_Arrays& HS_Arrays);

} // namespace sparse_format

#endif