#include "spar_hsr.h"

#include "module_parameter/parameter.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/module_tddft/td_velocity.h"
#include "spar_dh.h"
#include "spar_exx.h"
#include "spar_u.h"

void sparse_format::cal_HSR(const UnitCell& ucell,
                            const Parallel_Orbitals& pv,
                            LCAO_HS_Arrays& HS_Arrays,
                            const Grid_Driver& grid,
                            const int& current_spin,
                            const double& sparse_thr,
                            const int (&nmp)[3],
                            hamilt::Hamilt<std::complex<double>>* p_ham
#ifdef __EXX
                            ,
                            const std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd,
                            const std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc
#endif
) {
    ModuleBase::TITLE("sparse_format", "cal_HSR");

    sparse_format::set_R_range(HS_Arrays.all_R_coor, grid);

    const int nspin = PARAM.inp.nspin;

    // cal_STN_R_sparse(current_spin, sparse_thr);
    if (nspin == 1 || nspin == 2) {
        hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao
            = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(
                p_ham);

        if (TD_Velocity::tddft_velocity) {
            sparse_format::cal_HContainer_td(
                pv,
                current_spin,
                sparse_thr,
                *(p_ham_lcao->getHR()),
                TD_Velocity::td_vel_op->HR_sparse_td_vel[current_spin]);
        } else {
            sparse_format::cal_HContainer_d(pv,
                                            current_spin,
                                            sparse_thr,
                                            *(p_ham_lcao->getHR()),
                                            HS_Arrays.HR_sparse[current_spin]);
        }

        sparse_format::cal_HContainer_d(pv,
                                        current_spin,
                                        sparse_thr,
                                        *(p_ham_lcao->getSR()),
                                        HS_Arrays.SR_sparse);
    } else if (nspin == 4) {
        hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*
            p_ham_lcao
            = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>,
                                              std::complex<double>>*>(p_ham);

        sparse_format::cal_HContainer_cd(pv,
                                         current_spin,
                                         sparse_thr,
                                         *(p_ham_lcao->getHR()),
                                         HS_Arrays.HR_soc_sparse);

        sparse_format::cal_HContainer_cd(pv,
                                         current_spin,
                                         sparse_thr,
                                         *(p_ham_lcao->getSR()),
                                         HS_Arrays.SR_soc_sparse);
    } else {
        ModuleBase::WARNING_QUIT("cal_HSR", "check the value of nspin.");
    }

    // only old DFT+U method need to cal extra contribution to HR
    if (PARAM.inp.dft_plus_u == 2) {
        if (nspin == 1 || nspin == 2) {
            cal_HR_dftu(pv,
                        HS_Arrays.all_R_coor,
                        HS_Arrays.SR_sparse,
                        HS_Arrays.HR_sparse,
                        current_spin,
                        sparse_thr);
        } else if (nspin == 4) {
            cal_HR_dftu_soc(pv,
                            HS_Arrays.all_R_coor,
                            HS_Arrays.SR_soc_sparse,
                            HS_Arrays.HR_soc_sparse,
                            current_spin,
                            sparse_thr);
        } else {
            ModuleBase::WARNING_QUIT("cal_HSR", "check the value of nspin.");
        }
    }

#ifdef __EXX
#ifdef __MPI
    // if EXX is considered
    if (GlobalC::exx_info.info_global.cal_exx) {
        if (Hexxd && GlobalC::exx_info.info_ri.real_number)
        {
            sparse_format::cal_HR_exx(ucell,
                pv,
                HS_Arrays,
                current_spin,
                sparse_thr,
                nmp,
                *Hexxd);
        }
        else if (Hexxc && !GlobalC::exx_info.info_ri.real_number)
        {
            sparse_format::cal_HR_exx(ucell,
                pv,
                HS_Arrays,
                current_spin,
                sparse_thr,
                nmp,
                *Hexxc);
        }
    }
#endif // __MPI
#endif // __EXX

    sparse_format::clear_zero_elements(HS_Arrays, current_spin, sparse_thr);

    return;
}

void sparse_format::cal_HContainer_d(
    const Parallel_Orbitals& pv,
    const int& current_spin,
    const double& sparse_thr,
    const hamilt::HContainer<double>& hR,
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, double>>>& target) {
    ModuleBase::TITLE("sparse_format", "cal_HContainer_d");

    auto row_indexes = pv.get_indexes_row();
    auto col_indexes = pv.get_indexes_col();
    for (int iap = 0; iap < hR.size_atom_pairs(); ++iap) {
        int atom_i = hR.get_atom_pair(iap).get_atom_i();
        int atom_j = hR.get_atom_pair(iap).get_atom_j();
        int start_i = pv.atom_begin_row[atom_i];
        int start_j = pv.atom_begin_col[atom_j];
        int row_size = pv.get_row_size(atom_i);
        int col_size = pv.get_col_size(atom_j);
        for (int iR = 0; iR < hR.get_atom_pair(iap).get_R_size(); ++iR) {
            auto& matrix = hR.get_atom_pair(iap).get_HR_values(iR);
            const ModuleBase::Vector3<int> r_index
                = hR.get_atom_pair(iap).get_R_index(iR);
            Abfs::Vector3_Order<int> dR(r_index.x, r_index.y, r_index.z);
            for (int i = 0; i < row_size; ++i) {
                int mu = row_indexes[start_i + i];
                for (int j = 0; j < col_size; ++j) {
                    int nu = col_indexes[start_j + j];
                    const auto& value_tmp = matrix.get_value(i, j);
                    if (std::abs(value_tmp) > sparse_thr) {
                        target[dR][mu][nu] = value_tmp;
                    }
                }
            }
        }
    }

    return;
}

void sparse_format::cal_HContainer_cd(
    const Parallel_Orbitals& pv,
    const int& current_spin,
    const double& sparse_thr,
    const hamilt::HContainer<std::complex<double>>& hR,
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, std::complex<double>>>>&
        target) {
    ModuleBase::TITLE("sparse_format", "cal_HContainer_cd");

    auto row_indexes = pv.get_indexes_row();
    auto col_indexes = pv.get_indexes_col();
    for (int iap = 0; iap < hR.size_atom_pairs(); ++iap) {
        int atom_i = hR.get_atom_pair(iap).get_atom_i();
        int atom_j = hR.get_atom_pair(iap).get_atom_j();
        int start_i = pv.atom_begin_row[atom_i];
        int start_j = pv.atom_begin_col[atom_j];
        int row_size = pv.get_row_size(atom_i);
        int col_size = pv.get_col_size(atom_j);
        for (int iR = 0; iR < hR.get_atom_pair(iap).get_R_size(); ++iR) {
            auto& matrix = hR.get_atom_pair(iap).get_HR_values(iR);
            const ModuleBase::Vector3<int> r_index
                = hR.get_atom_pair(iap).get_R_index(iR);
            Abfs::Vector3_Order<int> dR(r_index.x, r_index.y, r_index.z);
            for (int i = 0; i < row_size; ++i) {
                int mu = row_indexes[start_i + i];
                for (int j = 0; j < col_size; ++j) {
                    int nu = col_indexes[start_j + j];
                    const auto& value_tmp = matrix.get_value(i, j);
                    if (std::abs(value_tmp) > sparse_thr) {
                        target[dR][mu][nu] = value_tmp;
                    }
                }
            }
        }
    }

    return;
}

void sparse_format::cal_HContainer_td(
    const Parallel_Orbitals& pv,
    const int& current_spin,
    const double& sparse_thr,
    const hamilt::HContainer<double>& hR,
    std::map<Abfs::Vector3_Order<int>,
             std::map<size_t, std::map<size_t, std::complex<double>>>>&
        target) {
    ModuleBase::TITLE("sparse_format", "cal_HContainer_td");

    auto row_indexes = pv.get_indexes_row();
    auto col_indexes = pv.get_indexes_col();
    for (int iap = 0; iap < hR.size_atom_pairs(); ++iap) {
        int atom_i = hR.get_atom_pair(iap).get_atom_i();
        int atom_j = hR.get_atom_pair(iap).get_atom_j();
        int start_i = pv.atom_begin_row[atom_i];
        int start_j = pv.atom_begin_col[atom_j];
        int row_size = pv.get_row_size(atom_i);
        int col_size = pv.get_col_size(atom_j);
        for (int iR = 0; iR < hR.get_atom_pair(iap).get_R_size(); ++iR) {
            auto& matrix = hR.get_atom_pair(iap).get_HR_values(iR);
            const ModuleBase::Vector3<int> r_index
                = hR.get_atom_pair(iap).get_R_index(iR);
            Abfs::Vector3_Order<int> dR(r_index.x, r_index.y, r_index.z);
            for (int i = 0; i < row_size; ++i) {
                int mu = row_indexes[start_i + i];
                for (int j = 0; j < col_size; ++j) {
                    int nu = col_indexes[start_j + j];
                    const auto& value_tmp
                        = std::complex<double>(matrix.get_value(i, j), 0.0);
                    if (std::abs(value_tmp) > sparse_thr) {
                        target[dR][mu][nu] += value_tmp;
                    }
                }
            }
        }
    }

    return;
}

// in case there are elements smaller than the threshold
void sparse_format::clear_zero_elements(LCAO_HS_Arrays& HS_Arrays,
                                        const int& current_spin,
                                        const double& sparse_thr) {
    ModuleBase::TITLE("sparse_format", "clear_zero_elements");

    if (PARAM.inp.nspin != 4) {
        for (auto& R_loop: HS_Arrays.HR_sparse[current_spin]) {
            for (auto& row_loop: R_loop.second) {
                auto& col_map = row_loop.second;
                auto iter = col_map.begin();
                while (iter != col_map.end()) {
                    if (std::abs(iter->second) <= sparse_thr) {
                        col_map.erase(iter++);
                    } else {
                        iter++;
                    }
                }
            }
        }
        if (TD_Velocity::tddft_velocity) {
            for (auto& R_loop:
                 TD_Velocity::td_vel_op->HR_sparse_td_vel[current_spin]) {
                for (auto& row_loop: R_loop.second) {
                    auto& col_map = row_loop.second;
                    auto iter = col_map.begin();
                    while (iter != col_map.end()) {
                        if (std::abs(iter->second) <= sparse_thr) {
                            col_map.erase(iter++);
                        } else {
                            iter++;
                        }
                    }
                }
            }
        }

        for (auto& R_loop: HS_Arrays.SR_sparse) {
            for (auto& row_loop: R_loop.second) {
                auto& col_map = row_loop.second;
                auto iter = col_map.begin();
                while (iter != col_map.end()) {
                    if (std::abs(iter->second) <= sparse_thr) {
                        col_map.erase(iter++);
                    } else {
                        iter++;
                    }
                }
            }
        }
    } else {
        for (auto& R_loop: HS_Arrays.HR_soc_sparse) {
            for (auto& row_loop: R_loop.second) {
                auto& col_map = row_loop.second;
                auto iter = col_map.begin();
                while (iter != col_map.end()) {
                    if (std::abs(iter->second) <= sparse_thr) {
                        col_map.erase(iter++);
                    } else {
                        iter++;
                    }
                } // end while iter
            }     // end row loop
        }         // end R loop

        for (auto& R_loop: HS_Arrays.SR_soc_sparse) {
            for (auto& row_loop: R_loop.second) {
                auto& col_map = row_loop.second;
                auto iter = col_map.begin();
                while (iter != col_map.end()) {
                    if (std::abs(iter->second) <= sparse_thr) {
                        col_map.erase(iter++);
                    } else {
                        iter++;
                    }
                } // end while iter
            }     // end row_loop
        }         // end R_loop
    }

    return;
}

void sparse_format::destroy_HS_R_sparse(LCAO_HS_Arrays& HS_Arrays) {
    ModuleBase::TITLE("sparse_format", "destroy_HS_R_sparse");

    if (PARAM.inp.nspin != 4) {
        std::map<Abfs::Vector3_Order<int>,
                 std::map<size_t, std::map<size_t, double>>>
            empty_HR_sparse_up;
        std::map<Abfs::Vector3_Order<int>,
                 std::map<size_t, std::map<size_t, double>>>
            empty_HR_sparse_down;
        std::map<Abfs::Vector3_Order<int>,
                 std::map<size_t, std::map<size_t, double>>>
            empty_SR_sparse;
        HS_Arrays.HR_sparse[0].swap(empty_HR_sparse_up);
        HS_Arrays.HR_sparse[1].swap(empty_HR_sparse_down);
        HS_Arrays.SR_sparse.swap(empty_SR_sparse);
    } else {
        std::map<Abfs::Vector3_Order<int>,
                 std::map<size_t, std::map<size_t, std::complex<double>>>>
            empty_HR_soc_sparse;
        std::map<Abfs::Vector3_Order<int>,
                 std::map<size_t, std::map<size_t, std::complex<double>>>>
            empty_SR_soc_sparse;
        HS_Arrays.HR_soc_sparse.swap(empty_HR_soc_sparse);
        HS_Arrays.SR_soc_sparse.swap(empty_SR_soc_sparse);
    }

    // 'all_R_coor' has a small memory requirement and does not need to be
    // deleted. std::set<Abfs::Vector3_Order<int>> empty_all_R_coor;
    // all_R_coor.swap(empty_all_R_coor);

    return;
}