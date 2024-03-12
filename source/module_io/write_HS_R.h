#ifndef WRITE_HS_R_H
#define WRITE_HS_R_H

#include "module_base/matrix.h"
#include "module_cell/klist.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_hamilt_general/hamilt.h"

namespace ModuleIO
{
        void output_HS_R(
            const int &istep,
            const ModuleBase::matrix& v_eff,
            LCAO_Hamilt &UHM,
            const K_Vectors& kv,
            hamilt::Hamilt<std::complex<double>>* p_ham,
            const std::string& SR_filename = "data-SR-sparse_SPIN0.csr",
            const std::string& HR_filename_up = "data-HR-sparse_SPIN0.csr",
            const std::string HR_filename_down = "data-HR-sparse_SPIN1.csr",
            const bool& binary = false,
            const double& sparse_threshold = 1e-10); //LiuXh add 2019-07-15, modify in 2021-12-3

        void output_dH_R(
            const int &istep,
            const ModuleBase::matrix& v_eff,
            LCAO_Hamilt &UHM,
            const K_Vectors& kv,
            const bool& binary = false,
            const double& sparse_threshold = 1e-10);

        void output_T_R(
            const int istep,
            LCAO_Hamilt &UHM,
            const std::string& TR_filename = "data-TR-sparse_SPIN0.csr",
            const bool& binary = false,
            const double& sparse_threshold = 1e-10);

        void output_S_R(
            LCAO_Hamilt &UHM,
            hamilt::Hamilt<std::complex<double>>* p_ham,
            const std::string& SR_filename = "data-SR-sparse_SPIN0.csr",
            const bool& binary = false,
            const double& sparse_threshold = 1e-10);
}

#endif
