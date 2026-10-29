#include <algorithm>
#include "diago_dav_subspace.h"

#include "diago_iter_assist.h"

#include "source_base/module_device/device.h"
#include "source_base/parallel_cell.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/parallel_reduce.h"
#include "source_base/module_container/base/third_party/lapack.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/kernels/dsp/dsp_connector.h"
// #include "source_base/module_container/ATen/kernels/lapack.h"

#include <ATen/kernels/lapack.h>

#include "source_hsolver/kernels/hegvd_op.h"
#include "source_hsolver/diag_hs_para.h"
#include "source_hsolver/kernels/bpcg_kernel_op.h" // normalize_op, precondition_op, apply_eigenvalues_op

#include <vector>
#include <limits>
#include <sstream>

#ifdef __MPI
#include <mpi.h>
#include "source_base/parallel_comm.h"
#include "source_base/parallel_device.h"
#endif

using namespace hsolver;

namespace
{
template <typename T>
T dav_conjugate(const T value)
{
    return value;
}

template <typename Real>
std::complex<Real> dav_conjugate(const std::complex<Real> value)
{
    return std::conj(value);
}

// Keep INFO and M by reference: the general connector currently passes M by value.
void dav_hegvx(const int n, const int nband, std::complex<double>* h, std::complex<double>* s,
               double* eigenvalues, std::complex<double>* vectors, std::complex<double>* work,
               const int lwork, double* rwork, int* iwork, int* ifail, int& found, int& info)
{
    const int itype = 1;
    const char job = 'V';
    const char range = 'I';
    const char uplo = 'U';
    const int first = 1;
    const double zero = 0;
    zhegvx_(&itype, &job, &range, &uplo, &n, h, &n, s, &n, &zero, &zero, &first, &nband,
            &zero, &found, eigenvalues, vectors, &n, work, &lwork, rwork, iwork, ifail, &info);
}

void dav_hegvx(const int n, const int nband, std::complex<float>* h, std::complex<float>* s,
               float* eigenvalues, std::complex<float>* vectors, std::complex<float>* work,
               const int lwork, float* rwork, int* iwork, int* ifail, int& found, int& info)
{
    const int itype = 1;
    const char job = 'V';
    const char range = 'I';
    const char uplo = 'U';
    const int first = 1;
    const float zero = 0;
    chegvx_(&itype, &job, &range, &uplo, &n, h, &n, s, &n, &zero, &zero, &first, &nband,
            &zero, &found, eigenvalues, vectors, &n, work, &lwork, rwork, iwork, ifail, &info);
}

void dav_hegvx(const int n, const int nband, double* h, double* s, double* eigenvalues,
               double* vectors, double* work, const int lwork, double* rwork, int* iwork,
               int* ifail, int& found, int& info)
{
    const int itype = 1;
    const char job = 'V';
    const char range = 'I';
    const char uplo = 'U';
    const int first = 1;
    const double zero = 0;
    dsygvx_(&itype, &job, &range, &uplo, &n, h, &n, s, &n, &zero, &zero, &first, &nband,
            &zero, &found, eigenvalues, vectors, &n, work, &lwork, iwork, ifail, &info);
}
} // namespace

template <typename T, typename Device>
Diago_DavSubspace<T, Device>::Diago_DavSubspace(const std::vector<Real>& precondition_in,
                                                const int& nband_in,
                                                const int& nbasis_in,
                                                const int& david_ndim_in,
                                                const double& diag_thr_in,
                                                const int& diag_nmax_in,
                                                const diag_comm_info& diag_comm_in,
                                                const int diag_subspace_in,
                                                const int diago_subspace_bs_in)
    : precondition(precondition_in), n_band(nband_in), dim(nbasis_in), nbase_x(nband_in * david_ndim_in),
      diag_thr(diag_thr_in), iter_nmax(diag_nmax_in), diag_comm(diag_comm_in),
        diag_subspace(diag_subspace_in), diago_subspace_bs(diago_subspace_bs_in)
{
    this->device = base_device::get_device_type(this->ctx);

    this->one = &one_;
    this->zero = &zero_;
    this->neg_one = &neg_one_;

    assert(david_ndim_in > 1);
    // Block ownership can be empty or uneven; use the actual communicator sum.
    ModuleBase::CommunicationDomain domain;
#ifdef __MPI
    domain.initialize(this->diag_comm.comm);
#endif
    const int global_dimension = domain.sum(this->dim);
    assert(this->nbase_x < global_dimension);
    assert(diag_subspace >= 0 && diag_subspace < 3);

    // TODO: Added memory usage statistics

    const int local_storage = std::max(this->nbase_x * this->dim, 1);
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    resmem_complex_op()(this->psi_in_iter, local_storage, "DAV::psi_in_iter");
    setmem_complex_op()(this->psi_in_iter, 0, this->nbase_x * this->dim);

    // the product of H and psi in the reduced psi set
    resmem_complex_op()(this->hpsi, local_storage, "DAV::hpsi");
    setmem_complex_op()(this->hpsi, 0, this->nbase_x * this->dim);

    // the product of S and psi in the reduced psi set
    resmem_complex_op()(this->spsi, local_storage, "DAV::spsi");
    setmem_complex_op()(this->spsi, 0, this->nbase_x * this->dim);

    // Hamiltonian on the reduced psi set
    resmem_complex_op()(this->hcc, this->nbase_x * this->nbase_x, "DAV::hcc");
    setmem_complex_op()(this->hcc, 0, this->nbase_x * this->nbase_x);

    // Overlap on the reduced psi set
    resmem_complex_op()(this->scc, this->nbase_x * this->nbase_x, "DAV::scc");
    setmem_complex_op()(this->scc, 0, this->nbase_x * this->nbase_x);

    // Eigenvectors
    resmem_complex_op()(this->vcc, this->nbase_x * this->nbase_x, "DAV::vcc");
    setmem_complex_op()(this->vcc, 0, this->nbase_x * this->nbase_x);
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        resmem_real_op()(this->d_precondition, nbasis_in);
        syncmem_var_h2d_op()(this->d_precondition, this->precondition.data(), nbasis_in);
        resmem_complex_op()(this->d_scc, this->nbase_x * this->nbase_x);
        resmem_real_op()(this->d_eigenvalue, this->nbase_x);
    }
#endif
}

template <typename T, typename Device>
Diago_DavSubspace<T, Device>::~Diago_DavSubspace()
{
    delmem_complex_op()(this->psi_in_iter);

    delmem_complex_op()(this->hpsi);
    delmem_complex_op()(this->spsi);
    delmem_complex_op()(this->hcc);
    delmem_complex_op()(this->scc);
    delmem_complex_op()(this->vcc);

#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        delmem_real_op()(this->d_precondition);
        delmem_complex_op()(this->d_scc);
        delmem_real_op()(this->d_eigenvalue);
    }
#endif
}

template <typename T, typename Device>
int Diago_DavSubspace<T, Device>::orthonormalize_cpu(const HSOperator<T, Device>& op,
                                                  const int first,
                                                  const int count)
{
    // Scale before testing rank, so small residuals are not discarded just for
    // their magnitude. Test the loss of S norm relative to the scaled input.
    const Real rank_tolerance = 64 * std::numeric_limits<Real>::epsilon();
    const int matrix_ld = std::max(this->dim, 1);
    // Keep valid storage on empty ranks while passing the physical size to H/S.
    std::vector<T> candidate(matrix_ld);
    std::vector<T> s_candidate(matrix_ld);
    std::vector<T> overlaps(this->nbase_x);
    int retained = 0;
    for (int column = 0; column < count; ++column)
    {
        const int source_offset = (first + column) * this->dim;
        Real scale = 0;
        int invalid = 0;
        for (int row = 0; row < this->dim; ++row)
        {
            candidate[row] = this->psi_in_iter[source_offset + row];
            const Real magnitude = std::abs(candidate[row]);
            invalid += !std::isfinite(magnitude);
            scale = std::max(scale, magnitude);
        }
        if (this->diag_comm.nproc > 1)
        {
            Parallel_Reduce::reduce_pool(&invalid, 1);
            double global_scale = scale;
            Parallel_Reduce::reduce_max_pool(this->diag_comm.nproc, global_scale);
            scale = static_cast<Real>(global_scale);
        }
        if (invalid != 0)
        {
            ModuleBase::WARNING_QUIT("Diago_DavSubspace", "Non-finite CPU basis direction.");
        }
        if (scale == 0)
        {
            continue;
        }
        for (int row = 0; row < this->dim; ++row)
        {
            candidate[row] /= scale;
        }
        op.spsi(candidate.data(), s_candidate.data(), this->dim, 1);
        T initial_norm = 0;
        for (int row = 0; row < this->dim; ++row)
        {
            initial_norm += dav_conjugate(candidate[row]) * s_candidate[row];
        }
        if (this->diag_comm.nproc > 1)
        {
            Parallel_Reduce::reduce_pool(&initial_norm, 1);
        }
        const Real initial_norm_squared = std::real(initial_norm);
        if (!std::isfinite(initial_norm_squared) || initial_norm_squared <= 0)
        {
            ModuleBase::WARNING_QUIT("Diago_DavSubspace", "The CPU basis has a non-positive or non-finite S norm.");
        }
        const int accepted = first + retained;
        for (int pass = 0; pass < 2; ++pass)
        {
            if (accepted > 0)
            {
                if (this->dim == 0)
                {
                    // The preceding collective stores global coefficients here;
                    // an empty rank must contribute zero again on the next pass.
                    const T zero_value = 0;
                    std::fill_n(overlaps.data(), accepted, zero_value);
                }
                ModuleBase::gemv_op<T, base_device::DEVICE_CPU>()
                    ('C', this->dim, accepted, this->one, this->psi_in_iter, matrix_ld,
                     s_candidate.data(), 1, this->zero, overlaps.data(), 1);
            }
            if (this->diag_comm.nproc > 1 && accepted > 0)
            {
                Parallel_Reduce::reduce_pool(overlaps.data(), accepted);
            }
            if (accepted > 0)
            {
                ModuleBase::gemv_op<T, base_device::DEVICE_CPU>()
                    ('N', this->dim, accepted, this->neg_one, this->psi_in_iter, matrix_ld,
                     overlaps.data(), 1, this->one, candidate.data(), 1);
            }
            // Reapply S to avoid estimating a tiny norm by cancellation of
            // large overlap terms, including when S is not the identity.
            op.spsi(candidate.data(), s_candidate.data(), this->dim, 1);
        }
        T final_norm = 0;
        for (int row = 0; row < this->dim; ++row)
        {
            final_norm += dav_conjugate(candidate[row]) * s_candidate[row];
        }
        if (this->diag_comm.nproc > 1)
        {
            Parallel_Reduce::reduce_pool(&final_norm, 1);
        }
        const Real final_norm_squared = std::real(final_norm);
        if (!std::isfinite(final_norm_squared))
        {
            ModuleBase::WARNING_QUIT("Diago_DavSubspace", "Non-finite projected CPU S norm.");
        }
        if (final_norm_squared <= rank_tolerance * initial_norm_squared)
        {
            continue;
        }
        const Real norm = std::sqrt(final_norm_squared);
        const int target_offset = accepted * this->dim;
        for (int row = 0; row < this->dim; ++row)
        {
            this->psi_in_iter[target_offset + row] = candidate[row] / norm;
            this->spsi[target_offset + row] = s_candidate[row] / norm;
        }
        ++retained;
    }
    return retained;
}

template <typename T, typename Device>
bool Diago_DavSubspace<T, Device>::restarted_basis_is_orthonormal_cpu() const
{
    // Ritz coefficients are S-orthonormal, so rotation normally preserves the
    // basis metric. Check the actual rotated vectors before using that identity.
    const int gram_size = this->n_band * this->n_band;
    const int matrix_ld = std::max(this->dim, 1);
    std::vector<T> gram(gram_size);
    ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()
        ('C', 'N', this->n_band, this->n_band, this->dim, this->one,
         this->psi_in_iter, matrix_ld, this->spsi, matrix_ld, this->zero,
         gram.data(), this->n_band);
    if (this->diag_comm.nproc > 1)
    {
        Parallel_Reduce::reduce_pool(gram.data(), gram_size);
    }
    const Real tolerance = 64 * std::numeric_limits<Real>::epsilon() * this->n_band;
    for (int column = 0; column < this->n_band; ++column)
    {
        for (int row = 0; row < this->n_band; ++row)
        {
            T error = gram[column * this->n_band + row];
            if (row == column)
            {
                error -= 1;
            }
            const Real magnitude = std::abs(error);
            if (!std::isfinite(magnitude) || magnitude > tolerance)
            {
                return false;
            }
        }
    }
    return true;
}

template <typename T, typename Device>
bool Diago_DavSubspace<T, Device>::residuals_converged_cpu(const int nbase,
                                                        const std::vector<Real>& eigenvalues,
                                                        const std::vector<double>& ethr_band) const
{
    // This guard is used only when rank filtering leaves no expansion. A
    // repeated Ritz value alone is not convergence in that situation.
    for (int band = 0; band < this->n_band; ++band)
    {
        Real norms[3] = {};
        for (int row = 0; row < this->dim; ++row)
        {
            T hx = 0;
            T sx = 0;
            for (int basis = 0; basis < nbase; ++basis)
            {
                const T coefficient = this->vcc[band * this->nbase_x + basis];
                hx += this->hpsi[basis * this->dim + row] * coefficient;
                sx += this->spsi[basis * this->dim + row] * coefficient;
            }
            const T residual = hx - eigenvalues[band] * sx;
            const T residual_squared = dav_conjugate(residual) * residual;
            const T h_squared = dav_conjugate(hx) * hx;
            const T s_squared = dav_conjugate(sx) * sx;
            norms[0] += std::real(residual_squared);
            norms[1] += std::real(h_squared);
            norms[2] += std::real(s_squared);
        }
        if (this->diag_comm.nproc > 1)
        {
            Parallel_Reduce::reduce_pool(norms, 3);
        }
        const Real residual_norm = std::sqrt(norms[0]);
        const Real operator_scale = std::sqrt(norms[1]) + std::abs(eigenvalues[band]) * std::sqrt(norms[2]);
        const Real roundoff = 32 * std::numeric_limits<Real>::epsilon() * nbase * operator_scale;
        const Real target = static_cast<Real>(ethr_band[band]);
        const Real tolerance = std::max(target, roundoff);
        if (!std::isfinite(residual_norm) || residual_norm > tolerance)
        {
            return false;
        }
    }
    return true;
}

template <typename T, typename Device>
int Diago_DavSubspace<T, Device>::diag_once(const HSOperator<T, Device>& op,
                                            T* psi_in,
                                            const int psi_in_dmax,
                                            Real* eigenvalue_in_hsolver,
                                            const std::vector<double>& ethr_band)
{
    ModuleBase::timer::start("Diago_DavSubspace", "diag_once");
    const int matrix_ld = std::max(this->dim, 1);
    const int psi_ld = std::max(psi_in_dmax, 1);

    // the eigenvalues in dav iter
    std::vector<Real> eigenvalue_iter(this->nbase_x, 0.0);

    // convflag[m] = true if the m th band is convergent
    std::vector<bool> convflag(this->n_band, false);

    // unconv[m] store the number of the m th unconvergent band
    std::vector<int> unconv(this->n_band);

    // the dimension of the reduced psi set
    int nbase = 0;

    // the number of the unconvergent bands
    this->notconv = this->n_band;

    ModuleBase::timer::start("Diago_DavSubspace", "first");

    syncmem_complex_2d_op()(this->psi_in_iter, this->dim, psi_in, psi_in_dmax, this->dim, this->n_band);
    for (int m = 0; m < this->n_band; m++)
    {
        unconv[m] = m;
    }

    if (this->device == base_device::CpuDevice)
    {
        const int retained = this->orthonormalize_cpu(op, 0, this->n_band);
        if (retained != this->n_band)
        {
            ModuleBase::WARNING_QUIT("Diago_DavSubspace", "Initial CPU trial vectors are linearly dependent in S.");
        }
    }

    // compute h*psi_in_iter
    // NOTE: bands after the first n_band should yield zero
    // hpsi[:, 0:nbase_x] = H * psi_in_iter[:, 0:nbase_x]
    op.hpsi(this->psi_in_iter, this->hpsi, this->dim, this->notconv);

    // compute s*psi_in_iter
    // spsi[:, 0:nbase_x] = S * psi_in_iter[:, 0:nbase_x]
    op.spsi(this->psi_in_iter, this->spsi, this->dim, this->notconv);

    this->cal_elem(this->dim, nbase, this->notconv, this->psi_in_iter, this->spsi, this->hpsi, this->hcc, this->scc);

    this->diag_zhegvx(nbase, this->notconv, this->hcc, this->scc, this->nbase_x, &eigenvalue_iter, this->vcc);

    for (size_t m = 0; m < this->n_band; m++)
    {
        eigenvalue_in_hsolver[m] = eigenvalue_iter[m];
    }

    ModuleBase::timer::end("Diago_DavSubspace", "first");

    int dav_iter = 0;

    do
    {
        dav_iter++;

        const int added = this->cal_grad(op,
                       this->dim,
                       nbase,
                       this->notconv,
                       this->psi_in_iter,
                       this->hpsi,
                       this->spsi,
                       this->vcc,
                       unconv.data(),
                       &eigenvalue_iter);

        if (added == 0)
        {
            const bool converged = this->residuals_converged_cpu(nbase, eigenvalue_iter, ethr_band);
            if (!converged)
            {
                ModuleBase::WARNING_QUIT("Diago_DavSubspace",
                                         "No independent CPU expansion direction remains, but the true residual is not converged.");
            }
            this->notconv = 0;
        }
        else
        {
            this->cal_elem(this->dim,
                           nbase,
                           added,
                           this->psi_in_iter,
                           this->spsi,
                           this->hpsi,
                           this->hcc,
                           this->scc);

            this->diag_zhegvx(nbase, this->n_band, this->hcc, this->scc, this->nbase_x, &eigenvalue_iter, this->vcc);

            // check convergence and update eigenvalues
            ModuleBase::timer::start("Diago_DavSubspace", "check_update");

            this->notconv = 0;
            for (int m = 0; m < this->n_band; m++)
            {
                const Real eigenvalue_change = eigenvalue_iter[m] - eigenvalue_in_hsolver[m];
                convflag[m] = (std::abs(eigenvalue_change) < ethr_band[m]);

                if (!convflag[m])
                {
                    unconv[this->notconv] = m;
                    this->notconv++;
                }

                eigenvalue_in_hsolver[m] = eigenvalue_iter[m];
            }

            ModuleBase::timer::end("Diago_DavSubspace", "check_update");
        }

        if ((this->notconv == 0) || (nbase + this->notconv + 1 > this->nbase_x) || (dav_iter == this->iter_nmax))
        {
            ModuleBase::timer::start("Diago_DavSubspace", "last");

            // updata eigenvectors of Hamiltonian
            setmem_complex_op()(psi_in, 0, n_band * psi_in_dmax);

#ifdef __DSP
            ModuleBase::gemm_op_mt<T, Device>() // In order to not coding another whole template, using this method to
                                                // minimize the code change.
#else
            ModuleBase::gemm_op<T, Device>()
#endif
                ('N',
                 'N',
                 this->dim,
                 this->n_band,
                 nbase,
                 this->one,
                 this->psi_in_iter,
                 matrix_ld,
                 this->vcc,
                 this->nbase_x,
                 this->zero,
                 psi_in,
                 psi_ld);

            if (!this->notconv || (dav_iter == this->iter_nmax))
            {
                // overall convergence or last iteration: exit the iteration

                ModuleBase::timer::end("Diago_DavSubspace", "last");
                break;
            }
            else
            {
                // if the dimension of the reduced basis set is becoming too large,
                // then replace the first N (=nband) basis vectors with the current
                // estimate of the eigenvectors and set the basis dimension to N;

                // update this->psi_in_iter according to psi_in
                syncmem_complex_2d_op()(this->psi_in_iter, this->dim, psi_in, psi_in_dmax, this->dim, this->n_band);

                // Rotate cached H and S images with the same Ritz coefficients.
                // This avoids reapplying H and solving an unchanged Ritz problem
                // at each restart. Repair the CPU basis only if its metric drifts.
                this->refresh(this->dim,
                              this->n_band,
                              nbase,
                              eigenvalue_in_hsolver,
                              this->psi_in_iter,
                              this->hpsi,
                              this->spsi,
                              this->hcc,
                              this->scc,
                              this->vcc);

                if (this->device == base_device::CpuDevice && !this->restarted_basis_is_orthonormal_cpu())
                {
                    const int retained = this->orthonormalize_cpu(op, 0, this->n_band);
                    if (retained != this->n_band)
                    {
                        ModuleBase::WARNING_QUIT("Diago_DavSubspace", "Restarted CPU Ritz vectors have lost S rank.");
                    }
                    op.hpsi(this->psi_in_iter, this->hpsi, this->dim, retained);
                    nbase = 0;
                    this->cal_elem(this->dim, nbase, retained, this->psi_in_iter,
                                   this->spsi, this->hpsi, this->hcc, this->scc);
                    this->diag_zhegvx(nbase, this->n_band, this->hcc, this->scc,
                                      this->nbase_x, &eigenvalue_iter, this->vcc);
                    for (int band = 0; band < this->n_band; ++band)
                    {
                        eigenvalue_in_hsolver[band] = eigenvalue_iter[band];
                    }
                }

                ModuleBase::timer::end("Diago_DavSubspace", "last");
            }
        }

    } while (true);

    ModuleBase::timer::end("Diago_DavSubspace", "diag_once");

    return dav_iter;
}

template <typename T, typename Device>
int Diago_DavSubspace<T, Device>::cal_grad(const HSOperator<T, Device>& op,
                                            const int& dim,
                                            const int& nbase,
                                            const int& notconv,
                                            T* psi_iter,
                                            T* hpsi,
                                            T* spsi,
                                            T* vcc,
                                            const int* unconv,
                                            std::vector<Real>* eigenvalue_iter)
{
    ModuleBase::timer::start("Diago_DavSubspace", "cal_grad");
    const int matrix_ld = std::max(this->dim, 1);

    if (this->device == base_device::CpuDevice)
    {
        std::vector<Real> expansion_eigenvalues(notconv);
        const int coefficient_count = nbase * notconv;
        std::vector<T> coefficients(coefficient_count);
        for (int column = 0; column < notconv; ++column)
        {
            const int band = unconv[column];
            const Real eigenvalue = (*eigenvalue_iter)[band];
            expansion_eigenvalues[column] = eigenvalue;
            for (int basis = 0; basis < nbase; ++basis)
            {
                coefficients[column * nbase + basis] = vcc[band * this->nbase_x + basis];
            }
        }
        T* residuals = psi_iter + nbase * dim;
        ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()
            ('N', 'N', dim, notconv, nbase, this->one, hpsi, matrix_ld,
             coefficients.data(), nbase, this->zero, residuals, matrix_ld);
        for (int column = 0; column < notconv; ++column)
        {
            const Real negative_eigenvalue = -expansion_eigenvalues[column];
            for (int basis = 0; basis < nbase; ++basis)
            {
                coefficients[column * nbase + basis] *= negative_eigenvalue;
            }
        }
        ModuleBase::gemm_op<T, base_device::DEVICE_CPU>()
            ('N', 'N', dim, notconv, nbase, this->one, spsi, matrix_ld,
             coefficients.data(), nbase, this->one, residuals, matrix_ld);
        precondition_op<T, base_device::DEVICE_CPU>()(dim, psi_iter, nbase, notconv,
                                                       this->precondition.data(), expansion_eigenvalues.data());
        const int added = this->orthonormalize_cpu(op, nbase, notconv);
        if (added > 0)
        {
            T* new_basis = psi_iter + nbase * dim;
            T* new_h_basis = hpsi + nbase * dim;
            op.hpsi(new_basis, new_h_basis, dim, added);
        }
        ModuleBase::timer::end("Diago_DavSubspace", "cal_grad");
        return added;
    }

    for (size_t i = 0; i < notconv; i++)
    {
        if (unconv[i] != i)
        {
            syncmem_complex_op()(vcc + i * this->nbase_x, vcc + unconv[i] * this->nbase_x, nbase);
            (*eigenvalue_iter)[i] = (*eigenvalue_iter)[unconv[i]];
        }
    }

    if (notconv > 1){

#ifdef __DSP
    ModuleBase::gemm_op_mt<T, Device>()
#else
    ModuleBase::gemm_op<T, Device>()
#endif
                        ('N',
                         'N',
                         this->dim,
                         notconv,
                         nbase,
                         this->one,
                         hpsi,
                         matrix_ld,
                         vcc,
                         this->nbase_x,
                         this->zero,
                         psi_iter + (nbase) * this->dim,
                         matrix_ld);
    } else
    {

#ifdef __DSP
    ModuleBase::gemv_op_mt<T, Device>()
#else
    ModuleBase::gemv_op<T, Device>()
#endif
                        ('N',
                         this->dim,     // m: row of A
                         nbase,         // n: col of A
                         this->one,     // alpha
                         hpsi,          // A dim * nbase
                         matrix_ld,     // LDA: if(N) max(1,m)
                         vcc,           // X nbase
                         1,             // incx
                         this->zero,    // beta
                         psi_iter + (nbase) * this->dim,  // Y dim
                         1              // incy
                        );
    }


    // Eigenvalues operation section
    Real* e_temp_hd = eigenvalue_iter->data();
    if (this->device == base_device::GpuDevice)
    {
        syncmem_var_h2d_op()(this->d_eigenvalue, eigenvalue_iter->data(), this->nbase_x);
        e_temp_hd = this->d_eigenvalue;
    }

    // vcc = - vcc * eigenvalue
    ModuleBase::matrix_mul_vector_op<T, Device>()(nbase, notconv, vcc, this->nbase_x, e_temp_hd, -1.0, vcc, this->nbase_x);

    if (notconv > 1){

#ifdef __DSP
    ModuleBase::gemm_op_mt<T, Device>()
#else
    ModuleBase::gemm_op<T, Device>()
#endif
        ('N',
         'N',
         this->dim,
         notconv,
         nbase,
         this->one,
         spsi,
         matrix_ld,
         vcc,
         this->nbase_x,
         this->one,
         psi_iter + nbase * this->dim,
         matrix_ld);
    } else
    {
#ifdef __DSP
    ModuleBase::gemv_op_mt<T, Device>()
#else
    ModuleBase::gemv_op<T, Device>()
#endif
        ('N',
         this->dim,     // m: row of A
         nbase,         // n: col of A
         this->one,     // alpha
         spsi,          // A dim * nbase
         matrix_ld,     // LDA: if(N) max(1,m)
         vcc,           // X nbase
         1,             // incx
         this->one,    // beta
         psi_iter + nbase * this->dim,  // Y dim
         1              // incy
        );
    }

    // Precondition section
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        precondition_op<T, Device>()(this->dim,
                                    psi_iter,
                                    nbase,
                                    notconv,
                                    d_precondition,
                                    this->d_eigenvalue);
    }
    else
#endif
    {
        precondition_op<T, Device>()(this->dim,
                                    psi_iter,
                                    nbase,
                                    notconv,
                                    this->precondition.data(),
                                    (*eigenvalue_iter).data());
    }

    // Normalize section
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        Real* psi_norm = nullptr;
        resmem_real_op()(psi_norm, notconv);
        setmem_real_op()(psi_norm, 0.0, notconv);

        normalize_op<T, Device>()(this->dim,
                                psi_iter,
                                nbase,
                                notconv,
                                psi_norm);

        // Check for zero norms (GPU path: copy norms from device to host)
        Real* psi_norm_host = nullptr;
        resmem_real_h_op()(psi_norm_host, notconv);
        syncmem_var_d2h_op()(psi_norm_host, psi_norm, notconv);
        for (int i = 0; i < notconv; i++)
        {
            if (psi_norm_host[i] <= 1.0e-12)
            {
                std::cout << "Diago_DavSubspace::cal_grad: psi_norm <= 0 for band " << i << std::endl;
                std::cout << "This may be due to npwx < nbands: the number of plane waves is less than" << std::endl;
                std::cout << "the number of bands, leading to a rank-deficient problem." << std::endl;
                std::cout << "Please increase ecutwfc or reduce nbands." << std::endl;
                delmem_real_h_op()(psi_norm_host);
                delmem_real_op()(psi_norm);
                ModuleBase::WARNING_QUIT("cal_grad", "psi_norm <= 0");
            }
        }
        delmem_real_h_op()(psi_norm_host);
        delmem_real_op()(psi_norm);
    }
    else
#endif
    {
        Real* psi_norm = nullptr;
        resmem_real_h_op()(psi_norm, notconv);
        setmem_real_h_op()(psi_norm, 0.0, notconv);

        normalize_op<T, Device>()(this->dim,
                                psi_iter,
                                nbase,
                                notconv,
                                psi_norm);

        // Check for zero norms (CPU path)
        for (int i = 0; i < notconv; i++)
        {
            if (psi_norm[i] <= 1.0e-12)
            {
                std::cout << "Diago_DavSubspace::cal_grad: psi_norm <= 0 for band " << i << std::endl;
                std::cout << "This may be due to npwx < nbands: the number of plane waves is less than" << std::endl;
                std::cout << "the number of bands, leading to a rank-deficient problem." << std::endl;
                std::cout << "Please increase ecutwfc or reduce nbands." << std::endl;
                delmem_real_h_op()(psi_norm);
                ModuleBase::WARNING_QUIT("cal_grad", "psi_norm <= 0");
            }
        }
        delmem_real_h_op()(psi_norm);
    }

    // update hpsi[:, nbase:nbase+notconv]
    // hpsi[:, nbase:nbase+notconv] = H * psi_iter[:, nbase:nbase+notconv]
    op.hpsi(psi_iter + nbase * dim, hpsi + nbase * this->dim, this->dim, notconv);
    op.spsi(psi_iter + nbase * dim, spsi + nbase * this->dim, this->dim, notconv);

    ModuleBase::timer::end("Diago_DavSubspace", "cal_grad");
    return notconv;
}

template <typename T, typename Device>
void Diago_DavSubspace<T, Device>::cal_elem(const int& dim,
                                            int& nbase,
                                            const int& notconv,
                                            const T* psi_iter,
                                            const T* spsi,
                                            const T* hpsi,
                                            T* hcc,
                                            T* scc)
{
    ModuleBase::timer::start("Diago_DavSubspace", "cal_elem");
    const int matrix_ld = std::max(this->dim, 1);

    if (this->dim == 0)
    {
        const int contribution_size = this->nbase_x * notconv;
        T* const h_contribution = hcc + nbase * this->nbase_x;
        T* const s_contribution = scc + nbase * this->nbase_x;
        // GEMV does not overwrite its output when the local row count is zero.
        setmem_complex_op()(h_contribution, 0, contribution_size);
        setmem_complex_op()(s_contribution, 0, contribution_size);
    }

    if (notconv > 1){
#ifdef __DSP
    ModuleBase::gemm_op_mt<T, Device>()
#else
    ModuleBase::gemm_op<T, Device>()
#endif
        ('C',
         'N',
         nbase + notconv,
         notconv,
         this->dim,
         this->one,
         psi_iter,
         matrix_ld,
         &hpsi[nbase * this->dim],
         matrix_ld,
         this->zero,
         &hcc[nbase * this->nbase_x],
         this->nbase_x);

#ifdef __DSP
    ModuleBase::gemm_op_mt<T, Device>()
#else
    ModuleBase::gemm_op<T, Device>()
#endif
        ('C',
         'N',
         nbase + notconv,
         notconv,
         this->dim,
         this->one,
         psi_iter,
         matrix_ld,
         spsi + nbase * this->dim,
         matrix_ld,
         this->zero,
         &scc[nbase * this->nbase_x],
         this->nbase_x);

        } else {

#ifdef __DSP
    ModuleBase::gemv_op_mt<T, Device>()
#else
    ModuleBase::gemv_op<T, Device>()
#endif
        ('C',
         this->dim,                 // m: row of A
         nbase + notconv,          // n: col of A
         this->one,                // alpha
         psi_iter,                 // A dim * nbase
         matrix_ld,                // LDA: if(N) max(1,m)
         &hpsi[nbase * this->dim], // X nbase
         1,                        // incx
         this->zero,               // beta
         &hcc[nbase * this->nbase_x], // Y dim
         1                         // incy
        );
#ifdef __DSP
    ModuleBase::gemv_op_mt<T, Device>()
#else
    ModuleBase::gemv_op<T, Device>()
#endif
        ('C',
         this->dim,                 // m: row of A
         nbase + notconv,          // n: col of A
         this->one,                // alpha
         psi_iter,                 // A dim * nbase
         matrix_ld,                // LDA: if(N) max(1,m)
         spsi + nbase * this->dim, // X nbase
         1,                        // incx
         this->zero,               // beta
         &scc[nbase * this->nbase_x], // Y dim
         1                         // incy
        );

        }


#ifdef __MPI
    if (this->diag_comm.nproc > 1)
    {
#ifdef __DSP
        // Only on dsp hardware need an extra space to reduce data
        mtfunc::dsp_dav_subspace_reduce(hcc, scc, nbase, this->nbase_x, notconv, this->diag_comm.comm);
#else
        assert(this->diag_comm.comm == POOL_WORLD);
        Parallel_Reduce::reduce_pool(hcc + nbase * this->nbase_x, notconv * this->nbase_x);
        Parallel_Reduce::reduce_pool(scc + nbase * this->nbase_x, notconv * this->nbase_x);
#endif
    }
#endif

    const size_t last_nbase = nbase; // init: last_nbase = 0
    nbase = nbase + notconv;

    ModuleBase::timer::end("Diago_DavSubspace", "cal_elem");
    return;
}

template <typename T, typename Device>
void Diago_DavSubspace<T, Device>::diag_zhegvx(const int& nbase,
                                               const int& nband,
                                               T* hcc,
                                               T* scc,
                                               const int& nbase_x,
                                               std::vector<Real>* eigenvalue_iter,
                                               T* vcc)
{
    ModuleBase::timer::start("Diago_DavSubspace", "diag_zhegvx");
    assert(nbase_x >= std::max(1, nbase));

    if (this->device == base_device::GpuDevice)
    {
#if defined(__CUDA) || defined(__ROCM)
        if (this->diag_comm.rank == 0)
        {
            syncmem_complex_op()(this->d_scc, scc, nbase * this->nbase_x);
            ct::kernels::lapack_hegvd<T, ct_Device>()(nbase, this->nbase_x, this->hcc, this->d_scc, this->d_eigenvalue, this->vcc);
            syncmem_var_d2h_op()((*eigenvalue_iter).data(), this->d_eigenvalue, this->nbase_x);
        }
#endif
    }
    else
    {
        if (this->diag_subspace == 0)
        {
            const int matrix_size = nbase * nbase;
            const int vector_size = nbase * nband;
            const T scalar_zero = 0;
            const Real real_zero = 0;
            std::vector<T> vectors(vector_size, scalar_zero);
            std::vector<Real> eigenvalues(nbase, real_zero);
            // Query INFO, solve INFO, and the actual returned eigenpair count.
            int status[3] = {};
            if (this->diag_comm.rank == 0)
            {
                std::vector<T> h_diag(matrix_size);
                std::vector<T> s_diag(matrix_size);
                for (int column = 0; column < nbase; ++column)
                {
                    for (int row = 0; row <= column; ++row)
                    {
                        const int source = column * nbase_x + row;
                        const int target = column * nbase + row;
                        h_diag[target] = hcc[source];
                        s_diag[target] = scc[source];
                        const Real h_magnitude = std::abs(h_diag[target]);
                        const Real s_magnitude = std::abs(s_diag[target]);
                        if (!std::isfinite(h_magnitude) || !std::isfinite(s_magnitude))
                        {
                            status[0] = -1000;
                        }
                    }
                }
                const int real_work_size = 7 * nbase;
                const int integer_work_size = 5 * nbase;
                std::vector<Real> rwork(real_work_size);
                std::vector<int> iwork(integer_work_size);
                std::vector<int> ifail(nbase);
                T work_query = 0;
                if (status[0] == 0)
                {
                    dav_hegvx(nbase, nband, h_diag.data(), s_diag.data(), eigenvalues.data(), vectors.data(),
                               &work_query, -1, rwork.data(), iwork.data(), ifail.data(), status[2], status[0]);
                }
                const Real requested_work = std::real(work_query);
                const Real maximum_work = static_cast<Real>(std::numeric_limits<int>::max());
                if (status[0] == 0
                    && (!std::isfinite(requested_work) || requested_work < 1 || requested_work >= maximum_work))
                {
                    status[0] = -1001;
                }
                if (status[0] == 0)
                {
                    const int queried_work_size = static_cast<int>(requested_work);
                    const int minimum_work_size = 8 * nbase;
                    const int work_size = std::max(queried_work_size, minimum_work_size);
                    std::vector<T> work(work_size);
                    status[2] = 0;
                    dav_hegvx(nbase, nband, h_diag.data(), s_diag.data(), eigenvalues.data(), vectors.data(),
                               work.data(), work_size, rwork.data(), iwork.data(), ifail.data(), status[2], status[1]);
                }
                if (status[0] == 0 && status[1] == 0 && status[2] == nband)
                {
                    for (int band = 0; band < nband; ++band)
                    {
                        if (!std::isfinite(eigenvalues[band]))
                        {
                            status[1] = -1002;
                        }
                    }
                    for (int index = 0; index < vector_size; ++index)
                    {
                        const Real magnitude = std::abs(vectors[index]);
                        if (!std::isfinite(magnitude))
                        {
                            status[1] = -1002;
                        }
                    }
                }
            }
            if (this->diag_comm.nproc > 1)
            {
                // Only the pool root contributes; publish status before any
                // rank decides whether the temporary Ritz result is usable.
                Parallel_Reduce::reduce_pool(status, 3);
            }
            if (status[0] != 0 || status[1] != 0 || status[2] != nband)
            {
                std::ostringstream message;
                message << "CPU reduced solve failed: query INFO=" << status[0]
                        << ", solve INFO=" << status[1] << ", M=" << status[2]
                        << ", requested=" << nband << ", dimension=" << nbase;
                const std::string diagnostic = message.str();
                ModuleBase::WARNING_QUIT("Diago_DavSubspace", diagnostic);
            }
            if (this->diag_comm.nproc > 1)
            {
                Parallel_Reduce::reduce_pool(eigenvalues.data(), nband);
                Parallel_Reduce::reduce_pool(vectors.data(), vector_size);
            }
            for (int band = 0; band < nband; ++band)
            {
                (*eigenvalue_iter)[band] = eigenvalues[band];
                for (int row = 0; row < nbase; ++row)
                {
                    vcc[band * nbase_x + row] = vectors[band * nbase + row];
                }
            }
            ModuleBase::timer::end("Diago_DavSubspace", "diag_zhegvx");
            return;
        }
        else
        {
#ifdef __MPI
            std::vector<T> h_diag;
            std::vector<T> s_diag;
            std::vector<T> vcc_tmp;
            if (this->diag_comm.rank == 0)
            {
                h_diag.resize(nbase * nbase, *this->zero);
                s_diag.resize(nbase * nbase, *this->zero);
                vcc_tmp.resize(nbase * nbase, *this->zero);
                for (size_t i = 0; i < nbase; i++)
                {
                    for (size_t j = 0; j < nbase; j++)
                    {
                        h_diag[i * nbase + j] = hcc[i * this->nbase_x + j];
                        s_diag[i * nbase + j] = scc[i * this->nbase_x + j];
                    }
                }
            }
            diago_hs_para(h_diag.data(),
                          s_diag.data(),
                          nbase,
                          nband,
                          (*eigenvalue_iter).data(),
                          vcc_tmp.data(),
                          this->diag_comm.comm,
                          this->diag_subspace,
                          this->diago_subspace_bs);
            if (this->diag_comm.rank == 0)
            {
                for (size_t i = 0; i < nband; i++)
                {
                    for (size_t j = 0; j < nbase; j++)
                    {
                        vcc[i * this->nbase_x + j] = vcc_tmp[i * nbase + j];
                    }
                }
            }
#else
            std::cout << "Error: parallel diagonalization is not supported in serial mode." << std::endl;
            exit(1);
#endif
        }
    }

#ifdef __MPI
    if (this->diag_comm.nproc > 1)
    {
        // vcc: nbase * nband
        for (int i = 0; i < nband; i++)
        {
            // fix a bug here when vcc is real: bcast according to the true type
            Parallel_Common::bcast_data(&vcc[i * this->nbase_x], nbase, this->diag_comm.comm, 0);
        }
        Parallel_Common::bcast_data((*eigenvalue_iter).data(), nband, this->diag_comm.comm, 0);
    }
#endif

    ModuleBase::timer::end("Diago_DavSubspace", "diag_zhegvx");
    return;
}

template <typename T, typename Device>
void Diago_DavSubspace<T, Device>::refresh(const int& dim,
                                           const int& nband,
                                           int& nbase,
                                           const Real* eigenvalue_in_hsolver,
                                           //    const psi::Psi<T, Device>& psi,
                                           T* psi_iter,
                                           T* hpsi,
                                           T* spsi,
                                           T* hcc,
                                           T* scc,
                                           T* vcc)
{
    ModuleBase::timer::start("Diago_DavSubspace", "refresh");
    const int matrix_ld = std::max(this->dim, 1);

#ifdef __DSP
    ModuleBase::gemm_op_mt<T, Device>()
#else
    ModuleBase::gemm_op<T, Device>()
#endif
                        ('N',
                         'N',
                         this->dim,
                         nband,
                         nbase,
                         this->one,
                         this->hpsi,
                         matrix_ld,
                         this->vcc,
                         this->nbase_x,
                         this->zero,
                         psi_iter + nband * this->dim,
                         matrix_ld);

    // update hpsi
    syncmem_complex_op()(hpsi, psi_iter + nband * this->dim, this->dim * nband);

#ifdef __DSP
    ModuleBase::gemm_op_mt<T, Device>()
#else
    ModuleBase::gemm_op<T, Device>()
#endif
        ('N',
         'N',
         this->dim,
         nband,
         nbase,
         this->one,
         this->spsi,
         matrix_ld,
         this->vcc,
         this->nbase_x,
         this->zero,
         psi_iter + nband * this->dim,
         matrix_ld);

    // update spsi
    syncmem_complex_op()(spsi, psi_iter + nband * this->dim, this->dim * nband);

    nbase = nband;

    // set hcc/scc/vcc to 0
    setmem_complex_2d_op()(hcc, this->nbase_x, 0, nbase, nbase);
    setmem_complex_2d_op()(scc, this->nbase_x, 0, nbase, nbase);
    setmem_complex_2d_op()(vcc, this->nbase_x, 0, nbase, nbase);

    if (this->device == base_device::GpuDevice)
    {
#if defined(__CUDA) || defined(__ROCM)
        syncmem_var_h2d_op()(this->d_eigenvalue, eigenvalue_in_hsolver, nbase);
#endif
        refresh_hcc_scc_vcc_op<T, Device>()(nbase, hcc, scc, vcc, this->nbase_x, this->d_eigenvalue, this->one_);
    }
    else
    {
        for (int i = 0; i < nbase; i++)
        {
            hcc[i * this->nbase_x + i] = eigenvalue_in_hsolver[i];
            scc[i * this->nbase_x + i] = this->one[0];
            vcc[i * this->nbase_x + i] = this->one[0];
        }
    }
    ModuleBase::timer::end("Diago_DavSubspace", "refresh");

    return;
}

template <typename T, typename Device>
int Diago_DavSubspace<T, Device>::diag(const HSOperator<T, Device>& op,
                                       T* psi_in,
                                       const int psi_in_dmax,
                                       Real* eigenvalue_in_hsolver,
                                       const std::vector<double>& ethr_band,
                                       const bool& scf_type)
{
    /// record the times of trying iterative diagonalization
    this->notconv = 0;

    int sum_iter = 0;
    int ntry = 0;

    do
    {

        sum_iter += this->diag_once(op, psi_in, psi_in_dmax, eigenvalue_in_hsolver, ethr_band);

        ++ntry;

    } while (this->test_exit_cond(ntry, this->notconv, scf_type));

    if (notconv > std::max(5, this->n_band / 4))
    {
        std::cout << "\n notconv = " << this->notconv;
        std::cout << "\n Diago_DavSubspace::diag', too many bands are not converged! \n";
    }

    return sum_iter;
}

template <typename T, typename Device>
bool Diago_DavSubspace<T, Device>::test_exit_cond(const int& ntry, const int& notconv, const bool& scf)
{
    // scf = true; // scf
    // scf = false; // nscf

    // If ntry <=5, try to do it better, if ntry > 5, exit.
    const bool f1 = (ntry <= 5);

    // In non-self consistent calculation, do until totally converged.
    const bool f2 = ((!scf && (notconv > 0)));

    // if self consistent calculation, if not converged > 5
    const bool f3 = ((scf && (notconv > 5)));

    return (f1 && (f2 || f3));
}

namespace hsolver
{

template class Diago_DavSubspace<std::complex<float>, base_device::DEVICE_CPU>;
template class Diago_DavSubspace<std::complex<double>, base_device::DEVICE_CPU>;

#if ((defined __CUDA) || (defined __ROCM))
template class Diago_DavSubspace<std::complex<float>, base_device::DEVICE_GPU>;
template class Diago_DavSubspace<std::complex<double>, base_device::DEVICE_GPU>;
#endif

#ifdef __LCAO
template class Diago_DavSubspace<double, base_device::DEVICE_CPU>;

#if ((defined __CUDA) || (defined __ROCM))
template class Diago_DavSubspace<double, base_device::DEVICE_GPU>;
#endif

#endif
} // namespace hsolver
