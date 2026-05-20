#include "output_mat_sparse.h"

#include "cal_r_overlap_R.h"
#include "source_io/module_hs/write_HS_R.h"
#include "source_io/module_parameter/parameter.h"

namespace ModuleIO
{
template <typename T>
void output_mat_sparse(const bool& out_mat_dh,
                       const bool& out_mat_ds,
                       const bool& out_mat_t,
                       const bool& out_mat_r,
                       const int& istep,
                       const ModuleBase::matrix& v_eff,
                       const Parallel_Orbitals& pv,
                       const TwoCenterBundle& two_center_bundle,
                       const LCAO_Orbitals& orb,
                       UnitCell& ucell,
                       const Grid_Driver& grid,
                       const K_Vectors& kv,
                       hamilt::Hamilt<T>* p_ham,
                       Plus_U* p_dftu)
{
    LCAO_HS_Arrays HS_Arrays; // store sparse arrays

    const int out_mat_hsR = PARAM.inp.out_mat_hs2[0];
    const bool reduce = (out_mat_hsR != 3 && out_mat_hsR != 4);
    const bool binary = (out_mat_hsR == 2 || out_mat_hsR == 4);
    //! generate a file containing the kinetic energy matrix
    if (out_mat_t)
    {
        output_TR(istep, ucell, pv, HS_Arrays, grid, two_center_bundle, orb, "trs1_nao.csr", binary, 1e-10, reduce);
    }

    //! generate a file containing the derivatives of the Hamiltonian matrix (in Ry/Bohr)
    if (out_mat_dh)
    {
        output_dHR(istep,
                   v_eff,
                   ucell,
                   pv,
                   HS_Arrays,
                   grid,
                   two_center_bundle,
                   orb,
                   kv,
                   binary,
                   1e-10,
                   reduce);
    }
    //! generate a file containing the derivatives of the overlap matrix (in Ry/Bohr)
    if (out_mat_ds)
    {
        output_dSR(istep,
                   ucell,
                   pv,
                   HS_Arrays,
                   grid,
                   two_center_bundle,
                   orb,
                   kv,
                   binary,
                   1e-10,
                   reduce);
    }

    // add by jingan for out r_R matrix 2019.8.14
    if (out_mat_r)
    {
        cal_r_overlap_R r_matrix;
        r_matrix.init(ucell, pv, orb);
        r_matrix.binary = binary;
        if (out_mat_hsR != 0)
        {
            r_matrix.out_rR_other(ucell, istep, HS_Arrays.output_R_coor, reduce);
        }
        else
        {
            r_matrix.out_rR(ucell, grid, istep, reduce);
        }
    }

    return;
}

template void output_mat_sparse<double>(const bool& out_mat_dh,
                                        const bool& out_mat_ds,
                                        const bool& out_mat_t,
                                        const bool& out_mat_r,
                                        const int& istep,
                                        const ModuleBase::matrix& v_eff,
                                        const Parallel_Orbitals& pv,
                                        const TwoCenterBundle& two_center_bundle,
                                        const LCAO_Orbitals& orb,
                                        UnitCell& ucell,
                                        const Grid_Driver& grid,
                                        const K_Vectors& kv,
										hamilt::Hamilt<double>* p_ham,
										Plus_U* p_dftu);

template void output_mat_sparse<std::complex<double>>(const bool& out_mat_dh,
                                                      const bool& out_mat_ds,
                                                      const bool& out_mat_t,
                                                      const bool& out_mat_r,
                                                      const int& istep,
                                                      const ModuleBase::matrix& v_eff,
                                                      const Parallel_Orbitals& pv,
                                                      const TwoCenterBundle& two_center_bundle,
                                                      const LCAO_Orbitals& orb,
                                                      UnitCell& ucell,
                                                      const Grid_Driver& grid,
                                                      const K_Vectors& kv,
													  hamilt::Hamilt<std::complex<double>>* p_ham,
													  Plus_U* p_dftu);

} // namespace ModuleIO
