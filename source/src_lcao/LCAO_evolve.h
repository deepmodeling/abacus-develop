#ifndef EVOLVE_LCAO_MATRIX_H
#define EVOLVE_LCAO_MATRIX_H

#include "../module_base/complexmatrix.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "module_hamilt/hamilt_lcao.h"
#include "module_psi/psi.h"
#include "src_lcao/LCAO_matrix.h"
#include "src_lcao/local_orbital_wfc.h"

class Evolve_LCAO_Matrix
{
  public:
    Evolve_LCAO_Matrix(const Parallel_Orbitals* pv)
    {
        this->ParaV = pv;
    }
    ~Evolve_LCAO_Matrix();

    void evolve_complex_matrix(const int& ik,
                               hamilt::Hamilt* phami,
                               Local_Orbital_wfc& lowf,
                               psi::Psi<std::complex<double>>* psi,
                               psi::Psi<std::complex<double>>* psi_laststep,
                               double* ekb) const;

  private:
    // LCAO_Matrix* LM;
    const Parallel_Orbitals* ParaV;

    void using_LAPACK_complex(const int& ik,
                              hamilt::Hamilt* phami,
                              std::complex<double>** wfc_k_grid,
                              std::complex<double>* wfc_k,
                              std::complex<double>* wfc_k_laststep,
                              Local_Orbital_wfc& lowf,
                              double* ekb) const;
#ifdef __MPI
    void using_ScaLAPACK_complex(const int& ik,
                                 hamilt::Hamilt* phami,
                                 std::complex<double>** wfc_k_grid,
                                 std::complex<double>* wfc_k,
                                 std::complex<double>* wfc_k_laststep,
                                 Local_Orbital_wfc& lowf,
                                 double* ekb) const;
#endif
};
#endif
