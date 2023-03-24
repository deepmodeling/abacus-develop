#ifndef EVOLVE_LCAO_MATRIX_H
#define EVOLVE_LCAO_MATRIX_H

#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_psi/psi.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"

class Evolve_LCAO_Matrix
{
  public:
    Evolve_LCAO_Matrix(const Parallel_Orbitals* pv)
    {
        this->ParaV = pv;
    }
    ~Evolve_LCAO_Matrix();

    void evolve_complex_matrix(const int& ik,
                               hamilt::Hamilt<double>* p_hamilt,
                               psi::Psi<std::complex<double>>* psi_k,
                               psi::Psi<std::complex<double>>* psi_k_laststep,
                               double* ekb) const;

  private:
    // LCAO_Matrix* LM;
    const Parallel_Orbitals* ParaV;

    void using_LAPACK_complex(const int& ik,
                              hamilt::Hamilt<double>* p_hamilt,
                              std::complex<double>* psi_k,
                              std::complex<double>* psi_k_laststep,
                              double* ekb) const;
#ifdef __MPI
    void using_ScaLAPACK_complex(
                  const int nband,
                  const int nlocal,         
                  const std::complex<double>* psi_k_laststep,
                  hamilt::Hamilt<double>* p_hamilt,
                  std::complex<double>* psi_k,
                  double* ekb) const;

    void compute_U_operator(
                  const int nband,
                  const int nlocal,     
                  const std::complex<double>* Stmp,
                  const std::complex<double>* Htmp,
                  std::complex<double>* U_operator,
                  const int print_matrix) const;

    void U_to_wfc(
                  const int nband,
                  const int nlocal,     
                  const std::complex<double>* U_operator,
                  const std::complex<double>* psi_k_laststep,
                  std::complex<double>* psi_k) const;

    void norm_wfc(
                  const int nband,
                  const int nlocal,     
                  const std::complex<double>* Stmp,
                  std::complex<double>* psi_k,
                  const int print_matrix) const;
          

    void compute_ekb(
                  const int nband,
                  const int nlocal,     
                  const std::complex<double>* Htmp,
                  const std::complex<double>* psi_k,
                  double* ekb) const;

#endif
};
#endif
