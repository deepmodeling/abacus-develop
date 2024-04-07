#ifndef OUTPUT_MAT_SPARSE_H
#define OUTPUT_MAT_SPARSE_H

#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_gen_fixedH.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "module_hsolver/hsolver_lcao.h"
#include "output_interface.h"

namespace ModuleIO
{
    /// @brief the output interface to write the sparse matrix of H, S, T, and r
    template<typename T>
    class Output_Mat_Sparse : public Output_Interface
{
  public:
    Output_Mat_Sparse(
        int out_mat_hsR,
        int out_mat_dh,
        int out_mat_t,
        int out_mat_r,
        int istep,
        const ModuleBase::matrix &v_eff,
        const Parallel_Orbitals &pv,
        LCAO_gen_fixedH &gen_h, // mohan add 2024-04-02
        Gint_k &gint_k, // mohan add 2024-04-01
        LCAO_Matrix &lm,
        Grid_Driver &grid, // mohan add 2024-04-06
        const K_Vectors &kv,
        hamilt::Hamilt<T> *p_ham);

    void write() override;

  private:

    //! generate a file containing the Hamiltonian and S(overlap) matrices 
    int _out_mat_hsR;

    //! generate a file containing the derivatives of the Hamiltonian matrix (in Ry/Bohr)
    int _out_mat_dh;

    //! generate a file containing the kinetic energy matrix
    int _out_mat_t;

    //! generate a file containing the matrix representation of the position matrix (in Bohr)
    int _out_mat_r;

    int _istep;

    const ModuleBase::matrix& _v_eff;

    const Parallel_Orbitals& _pv;

    LCAO_gen_fixedH& _gen_h; // mohan add 2024-04-02

    Gint_k& _gint_k; // mohan add 2024-04-01

    LCAO_Matrix& _lm;

    // mohan fix bug 2024-04-07, a typical bug!!!
    Grid_Driver& _grid; // mohan add 2024-04-06

    const K_Vectors& _kv;

    hamilt::Hamilt<T>* _p_ham;
};

} // namespace ModuleIO

#endif // OUTPUT_MAT_SPARSE_H
