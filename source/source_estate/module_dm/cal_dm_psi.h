#ifndef CAL_DM_PSI_H
#define CAL_DM_PSI_H

#include "source_base/matrix.h"
#include "source_psi/psi.h"
#include "density_matrix.h"
#include "source_lcao/module_gint/gint_helper.h"

namespace elecstate
{
    // for Gamma-Only case where DMK is double
    void cal_dm_psi(const Parallel_Orbitals* ParaV, const ModuleBase::matrix& wg, const psi::Psi<double>& wfc, elecstate::DensityMatrix<double, double>& DM);

    // for Multi-k case where DMK is std::complex<double>
    template <typename TR>
    void cal_dm_psi(const Parallel_Orbitals* ParaV, const ModuleBase::matrix& wg, const psi::Psi<std::complex<double>>& wfc, elecstate::DensityMatrix<std::complex<double>, TR>& DM);

    // for Multi-k case with mixed precision (fp32 GEMM internally)
    void cal_dm_psi_mixed(const Parallel_Orbitals* ParaV, const ModuleBase::matrix& wg, const psi::Psi<std::complex<double>>& wfc, elecstate::DensityMatrix<std::complex<double>, double>& DM);

    // for Gamma-Only case with MPI
    void psiMulPsiMpi(const psi::Psi<double>& psi1,
                            const psi::Psi<double>& psi2,
                            double* dm_out,
                            const int* desc_psi,
                            const int* desc_dm);
    
    // for multi-k case with MPI
    void psiMulPsiMpi(const psi::Psi<std::complex<double>>& psi1,
                            const psi::Psi<std::complex<double>>& psi2,
                            std::complex<double>* dm_out,
                            const int* desc_psi,
                            const int* desc_dm);

    // for multi-k case with MPI and fp32 GEMM
    void psiMulPsiMpiMixed(const psi::Psi<std::complex<double>>& psi1,
                                  const psi::Psi<std::complex<double>>& psi2,
                                  std::complex<double>* dm_out,
                                  const int* desc_psi,
                                  const int* desc_dm,
                                  int nloc_dm);

    // for Gamma-Only case without MPI
    void psiMulPsi(const psi::Psi<double>& psi1, const psi::Psi<double>& psi2, double* dm_out);

    // for multi-k case without MPI
    void psiMulPsi(const psi::Psi<std::complex<double>>& psi1,
                        const psi::Psi<std::complex<double>>& psi2,
                        std::complex<double>* dm_out);

    // for multi-k case without MPI with fp32 GEMM
    void psiMulPsiMixed(const psi::Psi<std::complex<double>>& psi1,
                              const psi::Psi<std::complex<double>>& psi2,
                              std::complex<double>* dm_out);
};
#endif
