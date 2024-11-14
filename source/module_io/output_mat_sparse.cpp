#include "output_mat_sparse.h"

#include "module_io/cal_r_overlap_R.h"
#include "module_io/write_HS_R.h"

namespace ModuleIO
{

template <>
void output_mat_sparse(const bool& out_mat_hsR,
                       const bool& out_mat_dh,
                       const bool& out_mat_t,
                       const bool& out_mat_r,
                       const int& istep,
                       const ModuleBase::matrix& v_eff,
                       const Parallel_Orbitals& pv,
                       Gint_k& gint_k,
                       const TwoCenterBundle& two_center_bundle,
                       const LCAO_Orbitals& orb,
                       UnitCell& ucell,
                       Grid_Driver& grid,
                       const K_Vectors& kv,
                       hamilt::Hamilt<double>* p_ham)
{
}

template <>
void output_mat_sparse(const bool& out_mat_hsR,
                       const bool& out_mat_dh,
                       const bool& out_mat_t,
                       const bool& out_mat_r,
                       const int& istep,
                       const ModuleBase::matrix& v_eff,
                       const Parallel_Orbitals& pv,
                       Gint_k& gint_k,
                       const TwoCenterBundle& two_center_bundle,
                       const LCAO_Orbitals& orb,
                       UnitCell& ucell,
                       Grid_Driver& grid,
                       const K_Vectors& kv,
                       hamilt::Hamilt<std::complex<double>>* p_ham)
{
    LCAO_HS_Arrays HS_Arrays; // store sparse arrays

    //! generate a file containing the Hamiltonian and S(overlap) matrices
    if (out_mat_hsR)
    {
        output_HSR(istep, v_eff, pv, HS_Arrays, grid, kv, p_ham);
    }

    //! generate a file containing the kinetic energy matrix
    if (out_mat_t)
    {
        output_TR(istep, ucell, pv, HS_Arrays, grid, two_center_bundle, orb);
    }

    //! generate a file containing the derivatives of the Hamiltonian matrix (in Ry/Bohr)
    if (out_mat_dh)
    {
        output_dHR(istep,
                   v_eff,
                   gint_k, // mohan add 2024-04-01
                   pv,
                   HS_Arrays,
                   grid, // mohan add 2024-04-06
                   two_center_bundle,
                   orb,
                   kv); // LiuXh add 2019-07-15
    }

    // add by jingan for out r_R matrix 2019.8.14
    if (out_mat_r)
    {
        cal_r_overlap_R r_matrix;
        r_matrix.init(pv, orb);
        if (out_mat_hsR)
        {
            r_matrix.out_rR_other(istep, HS_Arrays.output_R_coor);
        }
        else
        {
            r_matrix.out_rR(istep);
        }
    }

    return;
}

} // namespace ModuleIO
