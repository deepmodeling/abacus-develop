#include "output_mat_sparse.h"

#include "cal_r_overlap_R.h"
#include "source_io/module_hs/write_HS_R.h"

namespace ModuleIO
{
template <typename T>
void output_mat_sparse(const MatSparseOutputOptions& options,
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

    //! generate a file containing the kinetic energy matrix
    if (options.out_mat_t)
    {
        output_TR(istep,
                  ucell,
                  pv,
                  HS_Arrays,
                  grid,
                  two_center_bundle,
                  orb,
                  "trs1_nao.csr",
                  options.binary,
                  options.sparse_threshold,
                  options.t_precision);
    }

    //! generate a file containing the derivatives of the Hamiltonian matrix (in Ry/Bohr)
    if (options.out_mat_dh)
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
                   options.binary,
                   options.sparse_threshold,
                   options.dh_precision);
    }
    //! generate a file containing the derivatives of the overlap matrix (in Ry/Bohr)
    if (options.out_mat_ds)
    {
        output_dSR(istep,
                   ucell,
                   pv,
                   HS_Arrays,
                   grid,
                   two_center_bundle,
                   orb,
                   kv,
                   options.binary,
                   options.sparse_threshold,
                   options.ds_precision);
    }

    // add by jingan for out r_R matrix 2019.8.14
    if (options.out_mat_r)
    {
        cal_r_overlap_R r_matrix;
        r_matrix.binary = options.binary;
        r_matrix.sparse_threshold = options.sparse_threshold;
        r_matrix.init(ucell, pv, orb);
        r_matrix.out_rR(ucell, grid, istep, options.r_precision);
    }

    return;
}

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
    MatSparseOutputOptions options;
    options.out_mat_dh = out_mat_dh;
    options.out_mat_ds = out_mat_ds;
    options.out_mat_t = out_mat_t;
    options.out_mat_r = out_mat_r;
    output_mat_sparse(options,
                      istep,
                      v_eff,
                      pv,
                      two_center_bundle,
                      orb,
                      ucell,
                      grid,
                      kv,
                      p_ham,
                      p_dftu);
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

template void output_mat_sparse<double>(const MatSparseOutputOptions& options,
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

template void output_mat_sparse<std::complex<double>>(const MatSparseOutputOptions& options,
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
