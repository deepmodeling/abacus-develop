#ifndef OUTPUT_MAT_SPARSE_H
#define OUTPUT_MAT_SPARSE_H

#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hsolver/hsolver_lcao.h"
#include "output_interface.h"

namespace ModuleIO
{
/// @brief the output interface to write the sparse matrix of H, S, T, and r
class Output_Mat_Sparse : public Output_Interface
{
  public:
    Output_Mat_Sparse(int out_mat_hsR,
                      int out_mat_dh,
                      int out_mat_t,
                      int out_mat_r,
                      int istep,
                      const ModuleBase::matrix& v_eff,
                      const Parallel_Orbitals& pv,
                      LCAO_Hamilt& UHM,
                      LCAO_Matrix& LM,
                      const K_Vectors& kv,
                      hamilt::Hamilt<std::complex<double>>* p_ham);
    void write() override;

  private:
    int _out_mat_hsR;
    int _out_mat_dh;
    int _out_mat_t;
    int _out_mat_r;
    int _istep;
    const ModuleBase::matrix& _v_eff;
    const Parallel_Orbitals& _pv;
    LCAO_Hamilt& _UHM;
    LCAO_Matrix& _LM;
    const K_Vectors& _kv;
    hamilt::Hamilt<std::complex<double>>* _p_ham;
};
} // namespace ModuleIO

#endif // OUTPUT_MAT_SPARSE_H