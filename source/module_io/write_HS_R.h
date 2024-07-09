#ifndef WRITE_HS_R_H
#define WRITE_HS_R_H

#include "module_base/matrix.h"
#include "module_basis/module_nao/two_center_bundle.h"
#include "module_cell/klist.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace ModuleIO
{
void output_HSR(const int& istep,
                const ModuleBase::matrix& v_eff,
                const Parallel_Orbitals& pv,
                LCAO_Matrix& lm,
                LCAO_HS_Arrays& HS_Arrays,
                Grid_Driver& grid, // mohan add 2024-04-06
                const K_Vectors& kv,
                hamilt::Hamilt<std::complex<double>>* p_ham,
                const std::string& SR_filename = "data-SR-sparse_SPIN0.csr",
                const std::string& HR_filename_up = "data-HR-sparse_SPIN0.csr",
                const std::string HR_filename_down = "data-HR-sparse_SPIN1.csr",
                const bool& binary = false,
                const double& sparse_threshold = 1e-10); // LiuXh add 2019-07-15, modify in 2021-12-3

void output_dHR(const int& istep,
                const ModuleBase::matrix& v_eff,
                Gint_k& gint_k,    // mohan add 2024-04-01
                LCAO_Matrix& lm,   // mohan add 2024-04-01
                LCAO_HS_Arrays& HS_Arrays,
                Grid_Driver& grid, // mohan add 2024-04-06
                const TwoCenterBundle& two_center_bundle,
                const K_Vectors& kv,
                const bool& binary = false,
                const double& sparse_threshold = 1e-10);

void output_TR(const int istep,
               const UnitCell& ucell,
               const Parallel_Orbitals& pv,
               LCAO_Matrix& lm,
               LCAO_HS_Arrays& HS_Arrays,
               Grid_Driver& grid,
               const TwoCenterBundle& two_center_bundle,
               const std::string& TR_filename = "data-TR-sparse_SPIN0.csr",
               const bool& binary = false,
               const double& sparse_threshold = 1e-10);

void output_SR(Parallel_Orbitals& pv,
               LCAO_Matrix& lm,
               Grid_Driver& grid,
               hamilt::Hamilt<std::complex<double>>* p_ham,
               const std::string& SR_filename = "data-SR-sparse_SPIN0.csr",
               const bool& binary = false,
               const double& sparse_threshold = 1e-10);
} // namespace ModuleIO

#endif
