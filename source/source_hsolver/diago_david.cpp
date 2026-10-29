#include <algorithm>
#include "diago_david.h"

#include "source_base/parallel_cell.h"
#include "source_base/timer.h"                  // ModuleBase::timer::tick
#include "source_base/tool_title.h"             // ModuleBase::TITLE
#include "source_base/module_device/device.h"

#include "source_hsolver/kernels/hegvd_op.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_comm.h"

#include <cmath>
#include <limits>
#include <stdexcept>


using namespace hsolver;


template <typename T, typename Device>
DiagoDavid<T, Device>::DiagoDavid(const Real* precondition_in,
                                  const int nband_in,
                                  const int dim_in,
                                  const int david_ndim_in,
                                  const diag_comm_info& diag_comm_in)
    : nband(nband_in), dim(dim_in), nbase_x(david_ndim_in * nband_in), david_ndim(david_ndim_in), diag_comm(diag_comm_in)
{
    this->device = base_device::get_device_type(this->ctx);
    this->precondition = precondition_in;

    this->one = &one_;
    this->zero = &zero_;
    this->neg_one = &neg_one_;

    test_david = 2;
    // 1: check which function is called and which step is executed
    // 2: check the eigenvalues of the result of each iteration
    // 3: check the eigenvalues and errors of the last result
    // default: no check


    // set auxiliary memory
    // !!!
    // allocate operaions MUST be paired with deallocate operations
    // in the destructor
    assert(this->david_ndim > 1);
    // Block ownership can be empty or uneven; use the actual communicator sum.
    ModuleBase::CommunicationDomain domain;
#ifdef __MPI
    domain.initialize(this->diag_comm.comm);
#endif
    const int global_dimension = domain.sum(this->dim);
    assert(this->nbase_x < global_dimension);

    /// initialize variables
    /// k_first = 0 means that nks is more like a dimension of "basis" to be contracted in "HPsi".In LR-TDDFT the formula writes :
    /// $$\sum_{ jb\mathbf{k}'}A^I_{ia\mathbf{k}, jb\mathbf{k}' }X ^ I_{ jb\mathbf{k}'}$$
    /// In the code :
    /// - "H" means A
    /// - "Psi" means X
    /// - "band" means the superscript I : the number of excited states to be solved
    /// - k : k-points, the same meaning as the ground state
    /// - "basis" : number of occupied ks-orbitals(subscripts i,j) * number of unoccupied ks-orbitals(subscripts a,b), corresponding to "bands" of the ground state

    // the lowest N eigenvalues
    base_device::memory::resize_memory_op<Real, base_device::DEVICE_CPU>()(this->eigenvalue, nbase_x, "DAV::eig");
    base_device::memory::set_memory_op<Real, base_device::DEVICE_CPU>()(this->eigenvalue, 0, nbase_x);

    const int local_storage = std::max(this->nbase_x * this->dim, 1);
    // basis(dim, nbase_x), leading dimension = dim
    resmem_complex_op()(basis, local_storage, "DAV::basis");
    setmem_complex_op()(basis, 0, nbase_x * dim);

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // hpsi(nbase_x, dim); // the product of H and psi in the reduced basis set
    resmem_complex_op()(this->hpsi, local_storage, "DAV::hpsi");
    setmem_complex_op()(this->hpsi, 0, nbase_x * dim);

    // spsi(nbase_x, dim); // the Product of S and psi in the reduced basis set
    resmem_complex_op()(this->spsi, local_storage, "DAV::spsi");
    setmem_complex_op()(this->spsi, 0, nbase_x * dim);

    // hcc(nbase_x, nbase_x); // Hamiltonian on the reduced basis
    resmem_complex_op()(this->hcc, nbase_x * nbase_x, "DAV::hcc");
    setmem_complex_op()(this->hcc, 0, nbase_x * nbase_x);

    // vcc(nbase_x, nbase_x); // Eigenvectors of hcc
    resmem_complex_op()(this->vcc, nbase_x * nbase_x, "DAV::vcc");
    setmem_complex_op()(this->vcc, 0, nbase_x * nbase_x);
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // lagrange_matrix(nband, nband); // for orthogonalization
    resmem_complex_op()(this->lagrange_matrix, nband * nband);
    setmem_complex_op()(this->lagrange_matrix, 0, nband * nband);

#if defined(__CUDA) || defined(__ROCM)
    // device precondition array
    if (this->device == base_device::GpuDevice)
    {
        resmem_var_op()(this->d_precondition, dim);
        syncmem_var_h2d_op()(this->d_precondition, this->precondition, dim);
    }
#endif
}

template <typename T, typename Device>
DiagoDavid<T, Device>::~DiagoDavid()
{
    delmem_complex_op()(this->basis);
    delmem_complex_op()(this->hpsi);
    delmem_complex_op()(this->spsi);
    delmem_complex_op()(this->hcc);
    delmem_complex_op()(this->vcc);
    delmem_complex_op()(this->lagrange_matrix);
    base_device::memory::delete_memory_op<Real, base_device::DEVICE_CPU>()(this->eigenvalue);
    // If the device is a GPU device, free the d_precondition array.
#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        delmem_var_op()(this->d_precondition);
    }
#endif
}

template <typename T, typename Device>
int DiagoDavid<T, Device>::diag_once(const HSOperator<T, Device>& op,
                                     const int dim,
                                     const int nband,
                                     const int ld_psi,
                                     T *psi_in,
                                     Real* eigenvalue_in,
                                     const std::vector<double>& ethr_band,
                                     const int david_maxiter)
{
    if (test_david == 1)
    {
        ModuleBase::TITLE("DiagoDavid", "diag_once");
    }
    ModuleBase::timer::start("DiagoDavid", "diag_once");
    const int matrix_ld = std::max(dim, 1);
    const int psi_ld = std::max(ld_psi, 1);

    // convflag[m] = true if the m th band is converged
    std::vector<bool> convflag(nband, false);
    // unconv[m] store the number of the m th unconverged band
    std::vector<int> unconv(nband);

    int nbase = 0; // the dimension of the reduced basis set

    this->notconv = nband; // the number of unconverged eigenvalues

    for (int m = 0; m < nband; m++) {
        unconv[m] = m;
    }

    ModuleBase::timer::start("DiagoDavid", "first");

    // orthogonalise the initial trial psi(0~nband-1)

    // begin SchmidtOrth
    for (int m = 0; m < nband; m++)
    {
        syncmem_complex_op()(basis + dim*m, psi_in + m*ld_psi, dim);

        T* coefficients = this->lagrange_matrix + m * nband;
        const bool independent = this->SchmidtOrth(op, dim, m, this->spsi, coefficients);
        if (!independent)
        {
            throw std::runtime_error("DiagoDavid: initial wavefunctions are linearly dependent in the S metric");
        }
    }

    // end of SchmidtOrth and calculate H|psi>
    // hpsi_info dav_hpsi_in(&basis, psi::Range(true, 0, 0, nband - 1), this->hpsi);
    // phm_in->ops->hPsi(dav_hpsi_in);
    // hpsi[:, 0:nband] = H basis[:, 0:nband]
    // slice index in this piece of code is in C manner. i.e. 0:id stands for [0,id)
    op.hpsi(basis, hpsi, dim, nband);

    this->cal_elem(dim, nbase, nbase_x, this->notconv, this->hpsi, this->spsi, this->hcc);

    this->diag_zhegvx(nbase, nband, this->hcc, nbase_x, this->eigenvalue, this->vcc);

    for (int m = 0; m < nband; m++)
    {
        eigenvalue_in[m] = this->eigenvalue[m];
    }

    ModuleBase::timer::end("DiagoDavid", "first");

    int dav_iter = 0;
    do
    {
        dav_iter++;

        std::vector<bool> check_residual(nband, false);
        const int added = this->cal_grad(op,
                       dim,
                       nbase,
                       nbase_x,
                       this->notconv,
                       this->hpsi,
                       this->spsi,
                       this->vcc,
                       unconv.data(),
                       this->eigenvalue,
                       ethr_band,
                       check_residual);

        if (added > 0)
        {
            this->cal_elem(dim, nbase, nbase_x, added, this->hpsi, this->spsi, this->hcc);
            this->diag_zhegvx(nbase, nband, this->hcc, nbase_x, this->eigenvalue, this->vcc);
        }

        // check convergence and update eigenvalues
        ModuleBase::timer::start("DiagoDavid", "check_update");

        this->notconv = 0;
        for (int m = 0; m < nband; m++)
        {
            convflag[m] = (std::abs(this->eigenvalue[m] - eigenvalue_in[m]) < ethr_band[m]);
            if (check_residual[m])
            {
                // A skipped or dependent direction is not evidence of convergence.
                // Recheck the current Ritz state after the small solve: it may
                // have changed even if its previous residual was already small.
                const double residual_norm = this->ritz_residual_norm(nbase, m);
                convflag[m] = std::isfinite(residual_norm) && residual_norm <= ethr_band[m];
            }
            if (!convflag[m])
            {
                unconv[this->notconv] = m;
                this->notconv++;
            }
            eigenvalue_in[m] = this->eigenvalue[m];
        }

        ModuleBase::timer::end("DiagoDavid", "check_update");
        const bool expansion_stalled = added == 0 && this->notconv > 0;
        if (expansion_stalled && diag_comm.rank == 0)
        {
            std::cout << "DiagoDavid: no independent expansion direction remains for "
                      << this->notconv << " unconverged bands; retaining the residual failure for retry."
                      << std::endl;
        }
        if (!this->notconv || expansion_stalled || (nbase + this->notconv > nbase_x)
            || (dav_iter == david_maxiter))
        {
            ModuleBase::timer::start("DiagoDavid", "last");

            // update eigenvectors of Hamiltonian

            setmem_complex_op()(psi_in, 0, nband * ld_psi);
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            ModuleBase::gemm_op<T, Device>()('N',
                                             'N',
                                             dim,   // m: row of A,C
                                             nband, // n: col of B,C
                                             nbase, // k: col of A, row of B
                                             this->one,
                                             basis, // A dim * nbase
                                             matrix_ld,
                                             this->vcc, // B nbase * nband
                                             nbase_x,
                                             this->zero,
                                             psi_in, // C dim * nband
                                             psi_ld);

            if (!this->notconv || expansion_stalled || (dav_iter == david_maxiter))
            {
                // overall convergence or last iteration: exit the iteration
                ModuleBase::timer::end("DiagoDavid", "last");
                break;
            }
            else
            {
                // if the dimension of the reduced basis set is becoming too large,
                // then replace the first N (=nband) basis vectors with the current
                // estimate of the eigenvectors and set the basis dimension to N;
                this->refresh(dim,
                              nband,
                              nbase,
                              nbase_x,
                              eigenvalue_in,
                              psi_in,
                              ld_psi,
                              this->hpsi,
                              this->spsi,
                              this->hcc,
                              this->vcc);
                ModuleBase::timer::end("DiagoDavid", "last");
            }

        } // end of if

    } while (true);

    ModuleBase::timer::end("DiagoDavid", "diag_once");

    return dav_iter;
}


template <typename T, typename Device>
int DiagoDavid<T, Device>::cal_grad(const HSOperator<T, Device>& op,
                                        const int& dim,
                                        const int& nbase,   // current dimension of the reduced basis
                                        const int nbase_x,  // maximum dimension of the reduced basis set
                                        const int& notconv,
                                        T* hpsi,
                                        T* spsi,
                                        const T* vcc,
                                        const int* unconv,
                                        const Real* eigenvalue,
                                        const std::vector<double>& ethr_band,
                                        std::vector<bool>& check_residual)
{
    if (test_david == 1) {
        ModuleBase::TITLE("DiagoDavid", "cal_grad");
    }
    // if all eigenpairs have converged, return
    if (notconv == 0) {
        return 0;
    }
    ModuleBase::timer::start("DiagoDavid", "cal_grad");
    const int matrix_ld = std::max(dim, 1);

    // use template pointer for accelerate
    // std::complex<double>* spsi;
    // std::complex<double>* ppsi;

    // expand the reduced basis set with the new basis vectors P|Real(psi)>...
    // in which psi are the last eigenvectors
    // we define |Real(psi)> as (H-ES)*|Psi>, E = <psi|H|psi>/<psi|S|psi>

    // vc_ev_vector(notconv, nbase);
    // eigenvectors of unconverged index extracted from vcc
    T* vc_ev_vector = nullptr;
    resmem_complex_op()(vc_ev_vector, notconv * nbase);
    setmem_complex_op()(vc_ev_vector, 0, notconv * nbase);

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // for (int m = 0; m < notconv; m++)
    // {
    //     for (int i = 0; i < nbase; i++)
    //     {
    //         // vc_ev_vector(m, i) = vc(i, unconv[m]);
    //         vc_ev_vector[m * nbase + i] = vcc[i * nbase_x + unconv[m]];
    //     }
    // }
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    // set vc_ev_vector = vcc (unconverged index, leading dimension = nbase_x)
    // vc_ev_vector(m, i) = vc(i, unconv[m]);
    //         vc_ev_vector[m * nbase + i] = vcc[i * nbase_x + unconv[m]];
    for (int m = 0; m < notconv; m++)
    {
        syncmem_complex_op()(vc_ev_vector + m * nbase,
                             vcc + unconv[m] * nbase_x,
                             nbase);
    }

    // basis[nbase] = hpsi * vc_ev_vector = hpsi*vcc
    // basis'        =   vc_ev_vector' * hpsi'
    // (dim, notconv)  (dim, nbase) (nbase, notconv)
    if (notconv == 1){
        //Reuse gemv for vector case to avoid potential bug using gemm call with n=1
        ModuleBase::gemv_op<T, Device>()('N',
                                     dim,                 // m: row of A
                                     nbase,               // n: col of A
                                     this->one,           // alpha
                                     hpsi,                // A dim * nbase
                                     matrix_ld,                 // LDA: if(N) max(1,m)
                                     vc_ev_vector,        // X nbase
                                     1,                   // incx
                                     this->zero,          // beta
                                     basis + dim * nbase, // Y dim
                                     1                    // incy
        );

    }else
    {
        ModuleBase::gemm_op<T, Device>()('N',
                                     'N',
                                     dim,                 // m: row of A,C
                                     notconv,             // n: col of B,C
                                     nbase,               // k: col of A, row of B
                                     this->one,           // alpha
                                     hpsi,                // A dim * nbase
                                     matrix_ld,                 // LDA: if(N) max(1,m) if(T) max(1,k)
                                     vc_ev_vector,        // B nbase * notconv
                                     nbase,               // LDB: if(N) max(1,k) if(T) max(1,n)
                                     this->zero,          // belta
                                     basis + dim * nbase, // C dim * notconv
                                     matrix_ld                  // LDC: if(N) max(1, m)
        );
    }

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // for (int m = 0; m < notconv; m++)
    // {
    //     for (int i = 0; i < nbase; i++)
    //     {
    //         vc_ev_vector[m * nbase + i] *= -1 * this->eigenvalue[unconv[m]];
    //     }
    // }
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // e_temp_cpu = {-lambda}
    // vc_ev_vector[nbase] = vc_ev_vector[nbase] * e_temp_cpu
    // now vc_ev_vector[nbase] = - lambda * ev = -lambda * vcc
    for (int m = 0; m < notconv; m++)
    {
        std::vector<Real> e_temp_cpu(nbase, (-1.0 * this->eigenvalue[unconv[m]]));

        if (this->device == base_device::GpuDevice)
        {
#if defined(__CUDA) || defined(__ROCM)
            Real* e_temp_gpu = nullptr;
            resmem_var_op()(e_temp_gpu, nbase);
            syncmem_var_h2d_op()(e_temp_gpu, e_temp_cpu.data(), nbase);
            ModuleBase::vector_mul_vector_op<T, Device, Real>()(nbase,
                                                                vc_ev_vector + m * nbase,
                                                                vc_ev_vector + m * nbase,
                                                                e_temp_gpu,
                                                                false);
            delmem_var_op()(e_temp_gpu);
#endif
        }
        else
        {
            ModuleBase::vector_mul_vector_op<T, Device, Real>()(nbase,
                                                                vc_ev_vector + m * nbase,
                                                                vc_ev_vector + m * nbase,
                                                                e_temp_cpu.data(),
                                                                false);
        }
    }
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // basis[nbase] = basis[nbase] - spsi * vc_ev_vector
    //              = hpsi - spsi * lambda * vcc
    //              = (H - lambda * S) * psi * vcc
    //              = (H - lambda * S) * psi_new
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (notconv == 1){
        //Use gemv for vector case to avoid potential bug using gemm call with n=1
        ModuleBase::gemv_op<T, Device>()('N',
                                        dim,                  // m: row of A
                                        nbase,                // n: col of A
                                        this->one,            // alpha
                                        spsi,                 // A dim * nbase
                                        matrix_ld,                  // LDA: if(N) max(1,m)
                                        vc_ev_vector,         // X nbase
                                        1,                    // incx
                                        this->one,            // beta
                                        basis + dim * nbase,  // Y dim
                                        1                     //incy
        );
    } else 
    {
        ModuleBase::gemm_op<T, Device>()('N',
                                        'N',
                                        dim,                 // m: row of A,C
                                        notconv,             // n: col of B,C
                                        nbase,               // k: col of A, row of B
                                        this->one,           // alpha
                                        spsi,                // A
                                        matrix_ld,                 // LDA: if(N) max(1,m) if(T) max(1,k)
                                        vc_ev_vector,        // B
                                        nbase,               // LDB: if(N) max(1,k) if(T) max(1,n)
                                        this->one,           // beta
                                        basis + dim * nbase, // C dim * notconv
                                        matrix_ld                  // LDC: if(N) max(1, m)
        );
    }
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // Keep the number of accepted directions distinct from the number of
    // unconverged eigenpairs. Zero residuals and redundant block directions
    // must not introduce zero columns into the reduced problem.
    T* lagrange = nullptr;
    const int coefficient_count = nbase + notconv;
    resmem_complex_op()(lagrange, coefficient_count);
    int added = 0;
    for (int m = 0; m < notconv; m++)
    {
        const int band = unconv[m];
        T* candidate = basis + dim * (nbase + m);
        const double residual2 = ModuleBase::dot_real_op<T, Device>()(dim, candidate, candidate);
        if (!std::isfinite(residual2) || residual2 < 0.0)
        {
            throw std::runtime_error("DiagoDavid: non-finite residual before preconditioning");
        }
        const double residual_norm = std::sqrt(residual2);
        if (residual_norm <= ethr_band[band])
        {
            check_residual[band] = true;
            continue;
        }

        const int target = nbase + added;
        T* packed_candidate = basis + dim * target;
        if (added != m)
        {
            syncmem_complex_op()(packed_candidate, candidate, dim);
        }
        if (this->device == base_device::GpuDevice)
        {
#if defined(__CUDA) || defined(__ROCM)
            ModuleBase::vector_div_vector_op<T, Device, Real>()(dim,
                                                               packed_candidate,
                                                               packed_candidate,
                                                               this->d_precondition);
#endif
        }
        else
        {
            ModuleBase::vector_div_vector_op<T, Device, Real>()(dim,
                                                               packed_candidate,
                                                               packed_candidate,
                                                               this->precondition);
        }
        const bool independent = this->SchmidtOrth(op, dim, target, spsi, lagrange);
        if (independent)
        {
            ++added;
        }
        else
        {
            check_residual[band] = true;
        }
    }

    if (added > 0)
    {
        T* new_basis = basis + nbase * dim;
        T* new_hpsi = hpsi + nbase * dim;
        op.hpsi(new_basis, new_hpsi, dim, added);
    }
    delmem_complex_op()(lagrange);
    delmem_complex_op()(vc_ev_vector);

    ModuleBase::timer::end("DiagoDavid", "cal_grad");
    return added;
}


template <typename T, typename Device>
double DiagoDavid<T, Device>::ritz_residual_norm(const int nbase, const int band) const
{
    ModuleBase::timer::start("DiagoDavid", "ritz_residual_norm");
    const int matrix_ld = std::max(dim, 1);
    T* residual = nullptr;
    T* overlap = nullptr;
    resmem_complex_op()(residual, matrix_ld);
    resmem_complex_op()(overlap, matrix_ld);
    const T* coefficients = vcc + band * nbase_x;
    ModuleBase::gemv_op<T, Device>()('N', dim, nbase, this->one, hpsi, matrix_ld,
                                    coefficients, 1, this->zero, residual, 1);
    ModuleBase::gemv_op<T, Device>()('N', dim, nbase, this->one, spsi, matrix_ld,
                                    coefficients, 1, this->zero, overlap, 1);
    const Real unit = 1.0;
    const Real negative_eigenvalue = -eigenvalue[band];
    ModuleBase::vector_add_vector_op<T, Device, Real>()(dim, residual, residual, unit,
                                                        overlap, negative_eigenvalue);
    const double residual2 = ModuleBase::dot_real_op<T, Device>()(dim, residual, residual);
    delmem_complex_op()(residual);
    delmem_complex_op()(overlap);
    const double residual_norm = std::sqrt(residual2);
    ModuleBase::timer::end("DiagoDavid", "ritz_residual_norm");
    return residual_norm;
}


template <typename T, typename Device>
void DiagoDavid<T, Device>::cal_elem(const int& dim,
                                          int& nbase,           // current dimension of the reduced basis
                                          const int nbase_x,    // maximum dimension of the reduced basis set
                                          const int& notconv,   // number of newly added basis vectors
                                          const T* hpsi,
                                          const T* spsi,
                                          T* hcc)
{
    if (test_david == 1) {
        ModuleBase::TITLE("DiagoDavid", "cal_elem");
    }

    if (notconv == 0) {
        return;
    }
    ModuleBase::timer::start("DiagoDavid", "cal_elem");
    const int matrix_ld = std::max(dim, 1);

    // hcc[nbase](notconv, nbase + notconv)= basis[nbase]' * hpsi
    ModuleBase::gemm_op<T, Device>()('C',
                                     'N',
                                     notconv,
                                     nbase + notconv,
                                     dim,
                                     this->one,
                                     basis + dim * nbase, // basis(:,nbase:)  dim * notconv
                                     matrix_ld,
                                     hpsi, // dim * (nbase + notconv)
                                     matrix_ld,
                                     this->zero,
                                     hcc + nbase, // notconv * (nbase + notconv)
                                     nbase_x);


#ifdef __MPI
    if (diag_comm.nproc > 1)
    {
        ModuleBase::matrixTranspose_op<T, Device>()(nbase_x, nbase_x, hcc, hcc);

        assert(diag_comm.comm == POOL_WORLD);
        Parallel_Reduce::reduce_pool(hcc + nbase * nbase_x, notconv * nbase_x);

        ModuleBase::matrixTranspose_op<T, Device>()(nbase_x, nbase_x, hcc, hcc);
    }
#endif

    nbase += notconv;
    ModuleBase::timer::end("DiagoDavid", "cal_elem");
    return;
}

//==============================================================================
// optimize diag_zhegvx().
// 09-05-09 wangjp
// fixed a bug in diag_zhegvx().
// modify the dimension of h and s as (n,n) and copy the leading N*N
// part of hc & sc into h & s
// 09-05-10 wangjp
// As the complexmatrixs will be copied again in the subroutine ZHEGVX(...  ),
// i.e ZHEGVX(...) will not destroy the input complexmatrixs,
// we needn't creat another two complexmatrixs in diag_zhegvx().
//==============================================================================
template <typename T, typename Device>
void DiagoDavid<T, Device>::diag_zhegvx(const int& nbase,
                                             const int& nband,
                                             const T* hcc,
                                             const int& nbase_x,
                                             Real* eigenvalue, // in CPU
                                             T* vcc)
{
    ModuleBase::timer::start("DiagoDavid", "diag_zhegvx");
    if (diag_comm.rank == 0)
    {
        assert(nbase_x >= std::max(1, nbase));

        if (this->device == base_device::GpuDevice)
        {
#if defined(__CUDA) || defined(__ROCM)
            Real* eigenvalue_gpu = nullptr;
            resmem_var_op()(eigenvalue_gpu, nbase_x);
            syncmem_var_h2d_op()(eigenvalue_gpu, this->eigenvalue, nbase_x);

            // heevx_op<T, Device>()(this->ctx, nbase, nbase_x, hcc, nband, eigenvalue_gpu, vcc);
            ct::kernels::lapack_heevx<T, ct_Device>()(nbase, nbase_x, hcc, nband, eigenvalue_gpu, vcc);

            syncmem_var_d2h_op()(this->eigenvalue, eigenvalue_gpu, nbase_x);
            delmem_var_op()(eigenvalue_gpu);
#endif
        }
        else
        {
            //heevx_op<T, Device>()(this->ctx, nbase, nbase_x, hcc, nband, this->eigenvalue, vcc);
            ct::kernels::lapack_heevx<T, ct_Device>()(nbase, nbase_x, hcc, nband, this->eigenvalue, vcc);
        }
    }

#ifdef __MPI
    if (diag_comm.nproc > 1)
    {
        // vcc: nbase * nband
        for (int i = 0; i < nband; i++)
        {
            MPI_Bcast(&vcc[i * nbase_x], nbase, MPI_DOUBLE_COMPLEX, 0, diag_comm.comm);
        }
        MPI_Bcast(this->eigenvalue, nband, MPI_DOUBLE, 0, diag_comm.comm);
    }
#endif

    ModuleBase::timer::end("DiagoDavid", "diag_zhegvx");
    return;
}


template <typename T, typename Device>
void DiagoDavid<T, Device>::refresh(const int& dim,
                                         const int& nband,
                                         int& nbase,
                                         const int nbase_x, // maximum dimension of the reduced basis set
                                         const Real* eigenvalue_in,
                                         const T *psi_in,
                                         const int ld_psi,
                                         T* hpsi,
                                         T* spsi,
                                         T* hcc,
                                         T* vcc)
{
    if (test_david == 1) {
        ModuleBase::TITLE("DiagoDavid", "refresh");
    }
    ModuleBase::timer::start("DiagoDavid", "refresh");
    const int matrix_ld = std::max(dim, 1);

    // update hp,sp
    setmem_complex_op()(basis , 0, nbase_x * dim);

    // basis(dim, nband) = hpsi(dim, nbase) * vcc(nbase, nband)
    ModuleBase::gemm_op<T, Device>()('N',
                                     'N',
                                     dim,   // m: row of A,C
                                     nband, // n: col of B,C
                                     nbase, // k: col of A, row of B
                                     this->one,
                                     hpsi, // A dim * nbase
                                     matrix_ld,
                                     vcc, // B nbase * nband
                                     nbase_x,
                                     zero,
                                     basis, // C dim * nband
                                     matrix_ld);

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // basis[nband] = spsi * vcc
    ModuleBase::gemm_op<T, Device>()('N',
                                     'N',
                                     dim,   // m: row of A,C
                                     nband, // n: col of B,C
                                     nbase, // k: col of A, row of B
                                     this->one,
                                     spsi, // A dim * nbase
                                     matrix_ld,
                                     vcc, // B nbase * nband
                                     nbase_x,
                                     this->zero,
                                     basis + dim * nband, // C dim * nband
                                     matrix_ld);

    // hpsi = basis, spsi = basis[nband]
    syncmem_complex_op()(hpsi, basis, dim * nband);
    syncmem_complex_op()(spsi, basis + dim*nband, dim * nband);
    /*for (int m = 0; m < nband; m++) {
        for (int ig = 0; ig < dim; ig++)
        {
            hp(m, ig) = basis(m, ig);
            sp(m, ig) = basis(m + nband, ig);
        }
    }*/

    // update basis
    setmem_complex_op()(basis , 0, nbase_x * dim);

    for (int m = 0; m < nband; m++)
    {
        syncmem_complex_op()(basis + dim*m,psi_in + m*ld_psi, dim);
        /*for (int ig = 0; ig < npw; ig++)
            basis(m, ig) = psi(m, ig);*/
    }

    // update the reduced Hamiltonian
    // basis set size reset to nband
    nbase = nband;

    setmem_complex_op()(hcc, 0, nbase_x * nbase_x);
    // The restarted basis already consists of Ritz vectors. Its coefficient
    // matrix is the identity, including zero off-diagonal entries.
    const int coefficient_matrix_size = nbase_x * nbase_x;
    setmem_complex_op()(vcc, 0, coefficient_matrix_size);

    if (this->device == base_device::GpuDevice)
    {
#if defined(__CUDA) || defined(__ROCM)
        T* hcc_cpu = nullptr;
        T* vcc_cpu = nullptr;
        base_device::memory::resize_memory_op<T, base_device::DEVICE_CPU>()(hcc_cpu,
                                                                            nbase_x * nbase_x,
                                                                            "DAV::hcc");
        base_device::memory::resize_memory_op<T, base_device::DEVICE_CPU>()(vcc_cpu,
                                                                            nbase_x * nbase_x,
                                                                            "DAV::vcc");

        syncmem_d2h_op()(hcc_cpu, hcc, nbase_x * nbase_x);
        syncmem_d2h_op()(vcc_cpu, vcc, nbase_x * nbase_x);

        for (int i = 0; i < nbase; i++)
        {
            hcc_cpu[i * nbase_x + i] = eigenvalue_in[i];
            vcc_cpu[i * nbase_x + i] = this->one[0];
        }

        syncmem_h2d_op()(hcc, hcc_cpu, nbase_x * nbase_x);
        syncmem_h2d_op()(vcc, vcc_cpu, nbase_x * nbase_x);

        base_device::memory::delete_memory_op<T, base_device::DEVICE_CPU>()(hcc_cpu);
        base_device::memory::delete_memory_op<T, base_device::DEVICE_CPU>()(vcc_cpu);
#endif
    }
    else
    {
        for (int i = 0; i < nbase; i++)
        {
            hcc[i * nbase_x + i] = eigenvalue_in[i];
            // sc(i, i) = this->one;
            // vc(i, i) = this->one;
            vcc[i * nbase_x + i] = this->one[0];
        }
    }
    ModuleBase::timer::end("DiagoDavid", "refresh");
    return;
}


template <typename T, typename Device>
bool DiagoDavid<T, Device>::SchmidtOrth(const HSOperator<T, Device>& op,
                                       const int& dim,
                                       const int m,
                                       T* spsi,
                                       T* lagrange_m)
{
    ModuleBase::timer::start("DiagoDavid", "SchmidtOrth");
    const int matrix_ld = std::max(dim, 1);
    assert(m >= 0);
    T* candidate = basis + dim * m;
    T* s_candidate = spsi + dim * m;
    op.spsi(candidate, s_candidate, dim, 1);
    const Real initial_squared = ModuleBase::dot_real_op<T, Device>()(dim, candidate, s_candidate);
    if (!std::isfinite(initial_squared) || initial_squared < 0)
    {
        throw std::runtime_error("DiagoDavid: non-positive or non-finite overlap norm");
    }
    if (initial_squared == 0)
    {
        ModuleBase::timer::end("DiagoDavid", "SchmidtOrth");
        return false;
    }

    // Normalize the input before testing its independence: the absolute size
    // of a preconditioned residual does not measure the rank of its direction.
    const Real initial_scale = 1 / std::sqrt(initial_squared);
    ModuleBase::vector_mul_real_op<T, Device>()(dim, candidate, candidate, initial_scale);
    ModuleBase::vector_mul_real_op<T, Device>()(dim, s_candidate, s_candidate, initial_scale);
    for (int pass = 0; pass < 2; ++pass)
    {
        if (m > 0)
        {
            if (dim == 0)
            {
                // GEMV leaves its output untouched on empty local rows.
                // Clear each pass before all ranks reduce the projection.
                setmem_complex_op()(lagrange_m, 0, m);
            }
            ModuleBase::gemv_op<T, Device>()('C', dim, m, this->one, basis, matrix_ld,
                                            s_candidate, 1, this->zero, lagrange_m, 1);
            // Keep the collective host-side even in a device instantiation.
            std::vector<T> coefficients(m);
            syncmem_d2h_op()(coefficients.data(), lagrange_m, m);
            Parallel_Reduce::reduce_pool(coefficients.data(), m);
            syncmem_h2d_op()(lagrange_m, coefficients.data(), m);
            ModuleBase::gemv_op<T, Device>()('N', dim, m, this->neg_one, basis, matrix_ld,
                                            lagrange_m, 1, this->one, candidate, 1);
        }
        op.spsi(candidate, s_candidate, dim, 1);
    }

    // Measure the projected vector itself; subtracting nearly equal squared
    // norms loses the information needed to normalize a small direction.
    const Real projected_squared = ModuleBase::dot_real_op<T, Device>()(dim, candidate, s_candidate);
    if (!std::isfinite(projected_squared))
    {
        throw std::runtime_error("DiagoDavid: non-finite projected overlap norm");
    }
    const Real rank_tolerance = 64 * std::numeric_limits<Real>::epsilon();
    const Real rank_tolerance_squared = rank_tolerance * rank_tolerance;
    if (projected_squared <= rank_tolerance_squared)
    {
        ModuleBase::timer::end("DiagoDavid", "SchmidtOrth");
        return false;
    }
    const Real projected_scale = 1 / std::sqrt(projected_squared);
    ModuleBase::vector_mul_real_op<T, Device>()(dim, candidate, candidate, projected_scale);
    ModuleBase::vector_mul_real_op<T, Device>()(dim, s_candidate, s_candidate, projected_scale);
    ModuleBase::timer::end("DiagoDavid", "SchmidtOrth");
    return true;
}


template <typename T, typename Device>
void DiagoDavid<T, Device>::planSchmidtOrth(const int nband, std::vector<int>& pre_matrix_mm_m, std::vector<int>& pre_matrix_mv_m)
{
    if (nband <= 0) {
        return;
    }
    std::fill(pre_matrix_mm_m.begin(), pre_matrix_mm_m.end(), 0);
    std::fill(pre_matrix_mv_m.begin(), pre_matrix_mv_m.end(), 0);
    int last_matrix_size = nband;
    int matrix_size = int(nband / 2);
    int divide_times = 0;
    std::vector<int> divide_points(nband);
    int res_nband = nband - matrix_size;
    while (matrix_size > 1)
    {
        int index = nband - matrix_size;
        if (divide_times == 0)
        {
            divide_points[0] = index;
            pre_matrix_mm_m[index] = matrix_size;
            if (res_nband == matrix_size) {
                pre_matrix_mv_m[index] = 1;
            } else {
                pre_matrix_mv_m[index] = 2;
}
            divide_times = 1;
        }
        else
        {
            for (int i = divide_times - 1; i >= 0; i--)
            {
                divide_points[i * 2] = divide_points[i] - matrix_size;
                divide_points[i * 2 + 1] = divide_points[i * 2] + last_matrix_size;
                pre_matrix_mm_m[divide_points[i * 2]] = matrix_size;
                pre_matrix_mm_m[divide_points[i * 2 + 1]] = matrix_size;
                if (res_nband == matrix_size)
                {
                    pre_matrix_mv_m[divide_points[i * 2]] = 1;
                    pre_matrix_mv_m[divide_points[i * 2 + 1]] = 1;
                }
                else
                {
                    pre_matrix_mv_m[divide_points[i * 2]] = 2;
                    pre_matrix_mv_m[divide_points[i * 2 + 1]] = 2;
                }
            }
            divide_times *= 2;
        }
        last_matrix_size = matrix_size;
        matrix_size = int(res_nband / 2);
        res_nband -= matrix_size;
    }
    // fill the pre_matrix_mv_m array
    pre_matrix_mv_m[0] = 1;
    for (int m = 1; m < nband; m++)
    {
        if (pre_matrix_mv_m[m] == 0)
        {
            pre_matrix_mv_m[m] = pre_matrix_mv_m[m - 1] + 1;
        }
    }
}


template <typename T, typename Device>
int DiagoDavid<T, Device>::diag(const HSOperator<T, Device>& op,
                                const int ld_psi,
                                T *psi_in,
                                Real* eigenvalue_in,
                                const std::vector<double>& ethr_band,
                                const int david_maxiter,
                                const int ntry_max,
                                const int notconv_max)
{
    /// record the times of trying iterative diagonalization
    int ntry = 0;
    this->notconv = 0;

    int sum_dav_iter = 0;
    do
    {
        sum_dav_iter += this->diag_once(op, dim, nband, ld_psi, psi_in, eigenvalue_in, ethr_band, david_maxiter);
        ++ntry;
    } while (!check_block_conv(ntry, this->notconv, ntry_max, notconv_max));

    if (notconv > std::max(5, nband / 4))
    {
        std::cout << "\n notconv = " << this->notconv;
        std::cout << "\n DiagoDavid::diag', too many bands are not converged! \n";
    }
    return sum_dav_iter;
}


template <typename T, typename Device>
inline bool DiagoDavid<T, Device>::check_block_conv(const int& ntry,
                                                    const int& notconv,
                                                    const int& ntry_max,
                                                    const int& notconv_max)
{
    // Allow at most 5 tries at diag. If more than 5 then exit loop.
    if(ntry > ntry_max)
    {
        return true;
    }
    // If notconv <= notconv_max allowed, set convergence to true and exit loop.
    if(notconv <= notconv_max)
    {
        return true;
    }
    // else return false, continue loop until either condition above is met.
    return false;
}

namespace hsolver {
template class DiagoDavid<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoDavid<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoDavid<std::complex<float>, base_device::DEVICE_GPU>;
template class DiagoDavid<std::complex<double>, base_device::DEVICE_GPU>;
#endif
#ifdef __LCAO
template class DiagoDavid<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoDavid<double, base_device::DEVICE_GPU>;
#endif
#endif
} // namespace hsolver
