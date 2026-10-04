#ifndef CHG_PARALLEL_H
#define CHG_PARALLEL_H

// MPI reductions of charge-density arrays across k-point pools and band
// groups. Stateless free functions extracted from Charge member functions;
// the charge buffers and the parallel grid are supplied by the Charge
// argument. The pool/band parallelization settings are passed explicitly
// by the callers instead of reading GlobalV/PARAM.

#ifdef __MPI

class Charge;

namespace module_charge
{

/**
 * @brief Reduce a real-space array across k-point pools and band groups.
 *
 * @param array_rho real-space array [chr.nrxx], reduced in place
 * @param chr charge object supplying the parallel grid and the local grid size
 * @param kpar number of k-point pools (GlobalV::KPAR)
 * @param all_ks_run whether all processes run KS calculations (PARAM.globalv.all_ks_run)
 * @param bndpar number of band groups (PARAM.inp.bndpar)
 */
void reduce_diff_pools(double* array_rho, const Charge& chr, const int kpar,
                       const bool all_ks_run, const int bndpar);

/**
 * @brief Reduce rho across pools; also reduce kin_r when its buffer is
 *        allocated (meta-GGA functionals, or ELF output requested).
 *
 * @param chr charge object supplying rho/kin_r buffers
 * @param kpar number of k-point pools (GlobalV::KPAR)
 * @param all_ks_run whether all processes run KS calculations (PARAM.globalv.all_ks_run)
 * @param bndpar number of band groups (PARAM.inp.bndpar)
 * @param nspin number of spin channels (PARAM.inp.nspin)
 */
void rho_mpi(Charge& chr, const int kpar, const bool all_ks_run,
             const int bndpar, const int nspin);

/**
 * @brief Reduce kin_r across pools when its buffer is allocated
 *        (meta-GGA functionals, or ELF output requested).
 *
 * @param chr charge object supplying kin_r buffers
 * @param kpar number of k-point pools (GlobalV::KPAR)
 * @param all_ks_run whether all processes run KS calculations (PARAM.globalv.all_ks_run)
 * @param bndpar number of band groups (PARAM.inp.bndpar)
 * @param nspin number of spin channels (PARAM.inp.nspin)
 */
void kin_r_mpi(Charge& chr, const int kpar, const bool all_ks_run,
               const int bndpar, const int nspin);

} // namespace module_charge

#endif

#endif // CHG_PARALLEL_H
