#include "dm_from_psi.h"

#include <cassert>
#include <complex>
#include <vector>

#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"
#include "source_psi/psi.h"

namespace module_dm
{
namespace
{
/**
 * @brief Conjugation wrapper with an exact TK return type
 *
 * std::conj(double) is an additional overload whose availability/return type
 * differs across standard libraries; these helpers guarantee that the Gamma-only
 * (TK = double) instantiation stays real instead of promoting to complex.
 */
inline double conj_value(const double x)
{
    return x;
}

inline std::complex<double> conj_value(const std::complex<double> x)
{
    return std::conj(x);
}

/**
 * @brief Build the weighted (and conjugated) left factor of the DM GEMM for one k-point
 *
 * wg_wfc(ib, iw) = factor_ib * conj(wfc(ib, iw)); for TK = double std::conj is
 * the identity, so the same template covers the Gamma-only case.
 *
 * The local-to-global band mapping is taken verbatim from the historical
 * cal_dm implementation: ib_global advances monotonically while scanning
 * ParaV->global2local_col(). A local band whose global index does not fall
 * into the columns of wg keeps the legacy factor 1.0 instead of 0.0.
 *
 * @param ParaV orbital distribution, provides global2local_col()
 * @param wg band weights for the current k-point
 * @param ik k-point index; wfc must already be fixed to it
 * @param wfc wavefunction block of the current k-point
 * @param wg_wfc preallocated buffer of the same (nbands_local, nbasis_local) shape
 */
template <typename TK>
void fill_weighted_wfc(const Parallel_Orbitals* ParaV,
                       const ModuleBase::matrix& wg,
                       const int ik,
                       const psi::Psi<TK>& wfc,
                       psi::Psi<TK>& wg_wfc)
{
    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();

    // Resolve every per-band factor first: the global-band scan is serial.
    std::vector<double> factor(nbands_local, 1.0);
    int ib_global = 0;
    for (int ib_local = 0; ib_local < nbands_local; ++ib_local)
    {
        while (ib_local != ParaV->global2local_col(ib_global))
        {
            ++ib_global;
            if (ib_global >= wg.nc)
            {
                break;
            }
        }
        if (ib_global < wg.nc)
        {
            factor[ib_local] = wg(ik, ib_global);
        }
    }

    // Fuse the copy, the conjugation (complex case) and the weighting into one pass.
    // get_pointer() addresses the block selected by the caller's fix_k(); using the
    // three-argument operator() here would always address k-point block 0.
    const TK* wfc_block = wfc.get_pointer();
    TK* wg_wfc_block = wg_wfc.get_pointer();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ib_local = 0; ib_local < nbands_local; ++ib_local)
    {
        const int offset = ib_local * nbasis_local;
        const double factor_ib = factor[ib_local];
        for (int iw = 0; iw < nbasis_local; ++iw)
        {
            wg_wfc_block[offset + iw] = factor_ib * conj_value(wfc_block[offset + iw]);
        }
    }
}

#ifdef __MPI
/**
 * @brief Distributed GEMM: dmk = wg_wfc * wfc^T (column-major perspective)
 *
 * The row-major wavefunction wfc(ib, iw) is seen by the column-major BLAS as
 * its transpose. Using 'N' on the pre-conjugated wg_wfc and 'T' on wfc yields
 *     dmk(iw1, iw2) = sum_ib wg_wfc(ib, iw1) * wfc(ib, iw2),
 * i.e. the conjugation lives on the first index. 'C' must not be substituted
 * for 'T': 'T' and 'C' give identical operand dimensions, but 'C' additionally
 * conjugates the transposed operand, moving the conjugation to the second
 * index and producing the transpose of the stored conj-first DM block.
 */
void gemm_dm(const psi::Psi<double>& psi1,
             const psi::Psi<double>& psi2,
             double* dm_out,
             const int* desc_psi,
             const int* desc_dm)
{
    ModuleBase::timer::start("dmk_from_psi", "pdgemm");
    const double one_float = 1.0;
    const double zero_float = 0.0;
    const int one_int = 1;
    const char n_char = 'N';
    const char t_char = 'T';
    const int nlocal = desc_dm[2];
    const int nbands = desc_psi[3];
    ScalapackConnector::gemm(n_char,
                             t_char,
                             nlocal,
                             nlocal,
                             nbands,
                             one_float,
                             psi1.get_pointer(),
                             one_int,
                             one_int,
                             desc_psi,
                             psi2.get_pointer(),
                             one_int,
                             one_int,
                             desc_psi,
                             zero_float,
                             dm_out,
                             one_int,
                             one_int,
                             desc_dm);
    ModuleBase::timer::end("dmk_from_psi", "pdgemm");
}

void gemm_dm(const psi::Psi<std::complex<double>>& psi1,
             const psi::Psi<std::complex<double>>& psi2,
             std::complex<double>* dm_out,
             const int* desc_psi,
             const int* desc_dm)
{
    ModuleBase::timer::start("dmk_from_psi", "pzgemm");
    const std::complex<double> one_complex = {1.0, 0.0};
    const std::complex<double> zero_complex = {0.0, 0.0};
    const int one_int = 1;
    const char n_char = 'N';
    const char t_char = 'T';
    const int nlocal = desc_dm[2];
    const int nbands = desc_psi[3];
    ScalapackConnector::gemm(n_char,
                             t_char,
                             nlocal,
                             nlocal,
                             nbands,
                             one_complex,
                             psi1.get_pointer(),
                             one_int,
                             one_int,
                             desc_psi,
                             psi2.get_pointer(),
                             one_int,
                             one_int,
                             desc_psi,
                             zero_complex,
                             dm_out,
                             one_int,
                             one_int,
                             desc_dm);
    ModuleBase::timer::end("dmk_from_psi", "pzgemm");
}
#else
/**
 * @brief Serial GEMM: dmk = wg_wfc * wfc^T, see the MPI overload for the 'T' rationale
 */
void gemm_dm(const psi::Psi<double>& psi1, const psi::Psi<double>& psi2, double* dm_out)
{
    const double one_float = 1.0;
    const double zero_float = 0.0;
    const int one_int = 1;
    const char n_char = 'N';
    const char t_char = 'T';
    const int nlocal = psi1.get_nbasis();
    const int nbands = psi1.get_nbands();
    BlasConnector::gemm_cm(n_char,
                           t_char,
                           nlocal,
                           nlocal,
                           nbands,
                           one_float,
                           psi1.get_pointer(),
                           nlocal,
                           psi2.get_pointer(),
                           nlocal,
                           zero_float,
                           dm_out,
                           nlocal);
}

void gemm_dm(const psi::Psi<std::complex<double>>& psi1,
             const psi::Psi<std::complex<double>>& psi2,
             std::complex<double>* dm_out)
{
    const std::complex<double> one_complex = {1.0, 0.0};
    const std::complex<double> zero_complex = {0.0, 0.0};
    const int one_int = 1;
    const char n_char = 'N';
    const char t_char = 'T';
    const int nlocal = psi1.get_nbasis();
    const int nbands = psi1.get_nbands();
    BlasConnector::gemm_cm(n_char,
                           t_char,
                           nlocal,
                           nlocal,
                           nbands,
                           one_complex,
                           psi1.get_pointer(),
                           nlocal,
                           psi2.get_pointer(),
                           nlocal,
                           zero_complex,
                           dm_out,
                           nlocal);
}
#endif

/**
 * @brief Fill the weighted wavefunction buffer and run the DM GEMM for one k-point
 *
 * @param wg_wfc reusable scratch buffer, shape (1, nbands_local, nbasis_local)
 */
template <typename TK>
void dmk_from_psi_impl(const Parallel_Orbitals* ParaV,
                      const ModuleBase::matrix& wg,
                      const int ik,
                      const psi::Psi<TK>& wfc,
                      TK* dmk_out,
                      psi::Psi<TK>& wg_wfc)
{
    wfc.fix_k(ik);
    fill_weighted_wfc(ParaV, wg, ik, wfc, wg_wfc);
#ifdef __MPI
    gemm_dm(wg_wfc, wfc, dmk_out, ParaV->desc_wfc, ParaV->desc);
#else
    gemm_dm(wg_wfc, wfc, dmk_out);
#endif
}
} // namespace

// for Gamma-Only case where DMK is double
void dm_from_psi(const Parallel_Orbitals* ParaV,
                const ModuleBase::matrix& wg,
                const psi::Psi<double>& wfc,
                module_dm::DensityMatrix<double, double>& DM)
{
    assert(ParaV != nullptr);
    ModuleBase::TITLE("elecstate", "dm_from_psi");
    ModuleBase::timer::start("elecstate", "dm_from_psi");

    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();

    // Allocate the weighted-wavefunction scratch once and reuse it for every k-point.
    psi::Psi<double> wg_wfc(1, nbands_local, nbasis_local, nbasis_local, true);

    for (int ik = 0; ik < wfc.get_nk(); ++ik)
    {
        double* dmk_pointer = DM.get_dmk_ptr(ik);
        dmk_from_psi_impl(ParaV, wg, ik, wfc, dmk_pointer, wg_wfc);
    }
    ModuleBase::timer::end("elecstate", "dm_from_psi");
}

template <typename TR>
void dm_from_psi(const Parallel_Orbitals* ParaV,
                const ModuleBase::matrix& wg,
                const psi::Psi<std::complex<double>>& wfc,
                module_dm::DensityMatrix<std::complex<double>, TR>& DM)
{
    assert(ParaV != nullptr);
    ModuleBase::TITLE("elecstate", "dm_from_psi");
    ModuleBase::timer::start("elecstate", "dm_from_psi");

    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();

    // Allocate the weighted-wavefunction scratch once and reuse it for every k-point.
    psi::Psi<std::complex<double>> wg_wfc(1, nbands_local, nbasis_local, nbasis_local, true);

    for (int ik = 0; ik < wfc.get_nk(); ++ik)
    {
        std::complex<double>* dmk_pointer = DM.get_dmk_ptr(ik);
        dmk_from_psi_impl(ParaV, wg, ik, wfc, dmk_pointer, wg_wfc);
    }

    ModuleBase::timer::end("elecstate", "dm_from_psi");
}

void dmk_from_psi(const Parallel_Orbitals* ParaV,
                 const ModuleBase::matrix& wg,
                 const int ik,
                 const psi::Psi<double>& wfc,
                 double* dmk_out)
{
    assert(ParaV != nullptr);
    assert(dmk_out != nullptr);
    assert(ik >= 0 && ik < wfc.get_nk());

    psi::Psi<double> wg_wfc(1, wfc.get_nbands(), wfc.get_nbasis(), wfc.get_nbasis(), true);
    dmk_from_psi_impl(ParaV, wg, ik, wfc, dmk_out, wg_wfc);
}

void dmk_from_psi(const Parallel_Orbitals* ParaV,
                 const ModuleBase::matrix& wg,
                 const int ik,
                 const psi::Psi<std::complex<double>>& wfc,
                 std::complex<double>* dmk_out)
{
    assert(ParaV != nullptr);
    assert(dmk_out != nullptr);
    assert(ik >= 0 && ik < wfc.get_nk());

    psi::Psi<std::complex<double>> wg_wfc(1,
                                          wfc.get_nbands(),
                                          wfc.get_nbasis(),
                                          wfc.get_nbasis(),
                                          true);
    dmk_from_psi_impl(ParaV, wg, ik, wfc, dmk_out, wg_wfc);
}

template void dm_from_psi(const Parallel_Orbitals* ParaV,
                         const ModuleBase::matrix& wg,
                         const psi::Psi<std::complex<double>>& wfc,
                         module_dm::DensityMatrix<std::complex<double>, std::complex<double>>& DM);
template void dm_from_psi(const Parallel_Orbitals* ParaV,
                         const ModuleBase::matrix& wg,
                         const psi::Psi<std::complex<double>>& wfc,
                         module_dm::DensityMatrix<std::complex<double>, double>& DM);
} // namespace module_dm
