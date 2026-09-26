#include "output_mat_sparse.h"

#include "cal_r_overlap_r.h"
#include "source_io/module_hs/write_hs_r.h"
#include "source_io/module_parameter/parameter.h"

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
                       Plus_U_Base* p_dftu)
{
    LCAO_HS_Arrays HS_Arrays; // store sparse arrays

    const std::string& global_out_dir = PARAM.globalv.global_out_dir;
    const std::string& global_matrix_dir = PARAM.globalv.global_matrix_dir;
    const std::string& calculation = PARAM.inp.calculation;
    const bool out_app_flag = PARAM.inp.out_app_flag;
    const int nspin = PARAM.inp.nspin;
    const bool gamma_only_local = PARAM.globalv.gamma_only_local;
    const int npol = PARAM.globalv.npol;
    const int nlocal = PARAM.globalv.nlocal;

    //! generate a file containing the kinetic energy matrix
    if (options.out_mat_t)
    {
        const std::string tr_filename = "trs1_nao.csr";
        output_TR(istep,
                  ucell,
                  pv,
                  HS_Arrays,
                  grid,
                  two_center_bundle,
                  orb,
                  tr_filename,
                  options.binary,
                  options.sparse_threshold,
                  options.t_precision,
                  global_out_dir,
                  global_matrix_dir,
                  calculation,
                  out_app_flag,
                  nspin);
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
                   options.dh_precision,
                   global_out_dir,
                   global_matrix_dir,
                   calculation,
                   out_app_flag,
                   nspin,
                   gamma_only_local,
                   npol,
                   nlocal);
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
                   options.ds_precision,
                   global_out_dir,
                   global_matrix_dir,
                   calculation,
                   out_app_flag,
                   nspin,
                   gamma_only_local,
                   npol,
                   nlocal);
    }

    // add by jingan for out r_R matrix 2019.8.14
    if (options.out_mat_r)
    {
        cal_r_overlap_R r_matrix;
        r_matrix.binary = options.binary;
        r_matrix.sparse_threshold = options.sparse_threshold;
        const bool cal_force = PARAM.inp.cal_force;
        r_matrix.init(ucell, pv, orb, cal_force, nlocal);
        r_matrix.out_rR(ucell,
                        grid,
                        istep,
                        options.r_precision,
                        global_out_dir,
                        global_matrix_dir,
                        calculation,
                        out_app_flag,
                        nlocal,
                        npol);
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
                       Plus_U_Base* p_dftu)
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
                                        Plus_U_Base* p_dftu);

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
                                                      Plus_U_Base* p_dftu);

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
                                        Plus_U_Base* p_dftu);

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
                                                      Plus_U_Base* p_dftu);

} // namespace ModuleIO
