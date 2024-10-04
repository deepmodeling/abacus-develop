#include "diago_dav_subspace.h"

#include "diago_iter_assist.h"
#include "module_base/memory.h"
#include "module_base/module_device/device.h"
#include "module_base/timer.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"

#include <vector>

using namespace hsolver;

template <typename T, typename Device>
Diago_DavSubspace<T, Device>::Diago_DavSubspace(const std::vector<Real>& precondition_in,
                                                const int& nband_in,
                                                const int& nbasis_in,
                                                const int& david_ndim_in,
                                                const double& diag_thr_in,
                                                const int& diag_nmax_in,
                                                const bool& need_subspace_in,
                                                const diag_comm_info& diag_comm_in)
    : n_band(nband_in), dim(nbasis_in), nbase_x(nband_in * david_ndim_in),
      diag_thr(diag_thr_in), iter_nmax(diag_nmax_in), is_subspace(need_subspace_in), diag_comm(diag_comm_in), precondition(precondition_in)
{
    this->device = base_device::get_device_type<Device>(this->ctx);

    this->one = &this->cs.one;
    this->zero = &this->cs.zero;
    this->neg_one = &this->cs.neg_one;

    assert(david_ndim_in > 1);
    assert(david_ndim_in * nband_in < nbasis_in * this->diag_comm.nproc);

    // TODO: Added memory usage statistics

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    resmem_complex_op()(this->ctx, this->psi_in_iter, this->nbase_x * this->dim, "DAV::psi_in_iter");
    setmem_complex_op()(this->ctx, this->psi_in_iter, 0, this->nbase_x * this->dim);

    // the product of H and psi in the reduced psi set
    resmem_complex_op()(this->ctx, this->hphi, this->nbase_x * this->dim, "DAV::hphi");
    setmem_complex_op()(this->ctx, this->hphi, 0, this->nbase_x * this->dim);

    // Hamiltonian on the reduced psi set
    resmem_complex_op()(this->ctx, this->hcc, this->nbase_x * this->nbase_x, "DAV::hcc");
    setmem_complex_op()(this->ctx, this->hcc, 0, this->nbase_x * this->nbase_x);

    // Overlap on the reduced psi set
    resmem_complex_op()(this->ctx, this->scc, this->nbase_x * this->nbase_x, "DAV::scc");
    setmem_complex_op()(this->ctx, this->scc, 0, this->nbase_x * this->nbase_x);

    // Eigenvectors
    resmem_complex_op()(this->ctx, this->vcc, this->nbase_x * this->nbase_x, "DAV::vcc");
    setmem_complex_op()(this->ctx, this->vcc, 0, this->nbase_x * this->nbase_x);
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        resmem_real_op()(this->ctx, this->d_precondition, nbasis_in);
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, this->d_precondition, this->precondition.data(), nbasis_in);
    }
#endif
}

template <typename T, typename Device>
Diago_DavSubspace<T, Device>::~Diago_DavSubspace()
{
    delmem_complex_op()(this->ctx, this->psi_in_iter);

    delmem_complex_op()(this->ctx, this->hphi);
    delmem_complex_op()(this->ctx, this->hcc);
    delmem_complex_op()(this->ctx, this->scc);
    delmem_complex_op()(this->ctx, this->vcc);

#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        delmem_real_op()(this->ctx, this->d_precondition);
    }
#endif
}

template <typename T, typename Device>
int Diago_DavSubspace<T, Device>::diag_once(const HPsiFunc& hpsi_func,
                                            T* psi_in,
                                            const int psi_in_dmax,
                                            Real* eigenvalue_in_hsolver,
                                            const std::vector<bool>& is_occupied)
{
    ModuleBase::timer::tick("Diago_DavSubspace", "diag_once");

    // the eigenvalues in dav iter
    std::vector<Real> eigenvalue_iter(this->nbase_x, 0.0);

    // convflag[m] = true if the m th band is convergent
    std::vector<bool> convflag(this->n_band, false);

    // unconv[m] store the index of the m th unconvergent band
    std::vector<int> unconv(this->n_band);

    // the dimension of the reduced psi set
    int nbase = 0;

    // the number of the unconvergent bands
    this->notconv = this->n_band;

    ModuleBase::timer::tick("Diago_DavSubspace", "first");

    // copy psi_in to psi_in_iter for the first n_band bands
    // NOTE: not the first nbase_x (nband * david_ndim) bands
    // since psi_in_iter is zero-initialized, the rest bands are zero
    for (int m = 0; m < this->n_band; m++)
    {
        unconv[m] = m;

        syncmem_complex_op()(this->ctx,
                             this->ctx,
                             this->psi_in_iter + m * this->dim,
                             psi_in + m * psi_in_dmax,
                             this->dim);
    }

    // compute h*psi_in_iter
    // NOTE: bands after the first n_band should yield zero
    // hphi[:, 0:nbase_x] = H * psi_in_iter[:, 0:nbase_x]
    hpsi_func(this->psi_in_iter, this->hphi, this->dim, this->nbase_x);

    // at this stage, notconv = n_band and nbase = 0
    // note that nbase of cal_elem is an inout parameter: nbase := nbase + notconv
    this->cal_elem(this->dim, nbase, this->notconv, this->psi_in_iter, this->hphi, this->hcc, this->scc);

    // generalized eigenvalue problem (hcc, scc) for the first n_band bands
    this->diag_zhegvx(nbase,
                      this->n_band,
                      this->hcc,
                      this->scc,
                      this->nbase_x,
                      &eigenvalue_iter,
                      this->vcc,
                      true,
                      this->is_subspace);

    for (size_t m = 0; m < this->n_band; m++)
    {
        eigenvalue_in_hsolver[m] = eigenvalue_iter[m];
    }

    ModuleBase::timer::tick("Diago_DavSubspace", "first");

    int dav_iter = 0;

    do
    {
        dav_iter++;

        // nbase is the effective subspace size
        this->cal_grad(hpsi_func,
                       this->dim,
                       nbase,
                       this->notconv,
                       this->psi_in_iter,
                       this->hphi,
                       this->vcc,
                       unconv.data(),
                       &eigenvalue_iter);

        this->cal_elem(this->dim, nbase, this->notconv, this->psi_in_iter, this->hphi, this->hcc, this->scc);

        this->diag_zhegvx(nbase,
                          this->n_band,
                          this->hcc,
                          this->scc,
                          this->nbase_x,
                          &eigenvalue_iter,
                          this->vcc,
                          false,
                          false);

        // check convergence and update eigenvalues
        ModuleBase::timer::tick("Diago_DavSubspace", "check_update");

        this->notconv = 0;
        for (int m = 0; m < this->n_band; m++)
        {
            if (is_occupied[m]) // always true
            {
                convflag[m] = (std::abs(eigenvalue_iter[m] - eigenvalue_in_hsolver[m]) < this->diag_thr);
            }
            else
            {
                const double empty_ethr = std::max(this->diag_thr * 5.0, 1e-5);
                convflag[m] = (std::abs(eigenvalue_iter[m] - eigenvalue_in_hsolver[m]) < empty_ethr);
            }

            if (!convflag[m])
            {
                unconv[this->notconv] = m;
                this->notconv++;
            }

            eigenvalue_in_hsolver[m] = eigenvalue_iter[m];
        }

        ModuleBase::timer::tick("Diago_DavSubspace", "check_update");

        if ((this->notconv == 0) || (nbase + this->notconv + 1 > this->nbase_x) || (dav_iter == this->iter_nmax))
        {
            ModuleBase::timer::tick("Diago_DavSubspace", "last");

            // updata eigenvectors of Hamiltonian
            setmem_complex_op()(this->ctx, psi_in, 0, n_band * psi_in_dmax);

            gemm_op<T, Device>()(this->ctx,
                                 'N',
                                 'N',
                                 this->dim,
                                 this->n_band,
                                 nbase,
                                 this->one,
                                 this->psi_in_iter,
                                 this->dim,
                                 this->vcc,
                                 this->nbase_x,
                                 this->zero,
                                 psi_in,
                                 psi_in_dmax);

            if (!this->notconv || (dav_iter == this->iter_nmax))
            {
                // overall convergence or last iteration: exit the iteration

                // update this->psi_in_iter according to psi_in
                for (size_t i = 0; i < this->n_band; i++)
                {
                    syncmem_complex_op()(this->ctx,
                                         this->ctx,
                                         this->psi_in_iter + i * this->dim,
                                         psi_in + i * psi_in_dmax,
                                         this->dim);
                }

                this->refresh(this->dim,
                              this->n_band,
                              nbase,
                              eigenvalue_in_hsolver,
                              this->psi_in_iter,
                              this->hphi,
                              this->hcc,
                              this->scc,
                              this->vcc);

                ModuleBase::timer::tick("Diago_DavSubspace", "last");
                break;
            }
            else
            {
                // if the dimension of the reduced basis set is becoming too large,
                // then replace the first N (=nband) basis vectors with the current
                // estimate of the eigenvectors and set the basis dimension to N;

                // update this->psi_in_iter according to psi_in
                for (size_t i = 0; i < this->n_band; i++)
                {
                    syncmem_complex_op()(this->ctx,
                                         this->ctx,
                                         this->psi_in_iter + i * this->dim,
                                         psi_in + i * psi_in_dmax,
                                         this->dim);
                }

                this->refresh(this->dim,
                              this->n_band,
                              nbase,
                              eigenvalue_in_hsolver,
                              this->psi_in_iter,
                              this->hphi,
                              this->hcc,
                              this->scc,
                              this->vcc);

                ModuleBase::timer::tick("Diago_DavSubspace", "last");
            }
        }

    } while (true);

    ModuleBase::timer::tick("Diago_DavSubspace", "diag_once");

    return dav_iter;
}

template <typename T, typename Device>
void Diago_DavSubspace<T, Device>::cal_grad(const HPsiFunc& hpsi_func,
                                            const int& dim,
                                            const int& nbase,
                                            const int& notconv,
                                            T* psi_iter,
                                            T* hphi,
                                            T* vcc,
                                            const int* unconv,
                                            std::vector<Real>* eigenvalue_iter)
{
    ModuleBase::timer::tick("Diago_DavSubspace", "cal_grad");

    for (size_t i = 0; i < notconv; i++)
    {
        if (unconv[i] != i)
        {
            syncmem_complex_op()(this->ctx, this->ctx, vcc + i * this->nbase_x, vcc + unconv[i] * this->nbase_x, nbase);
            (*eigenvalue_iter)[i] = (*eigenvalue_iter)[unconv[i]];
        }
    }

    // psi_iter[:, nbase:nbase+notconv]
    //  = psi_iter[:, :nbase] * vcc[:nbase, :notconv]
    gemm_op<T, Device>()(this->ctx,
                         'N',
                         'N',
                         this->dim,
                         notconv,
                         nbase,
                         this->one,
                         psi_iter,
                         this->dim,
                         vcc,
                         this->nbase_x,
                         this->zero,
                         psi_iter + nbase * this->dim,
                         this->dim);

    // for psi_iter[:, nbase:nbase+notconv],
    // each column is multiplied by the corresponding (minus) eigenvalue
    // NOTE: eigenvalue_iter[m] correspond to psi_iter[:, nbase+m] (to be verified)
    for (int m = 0; m < notconv; m++)
    {

        std::vector<Real> e_temp_cpu(this->dim, ((-1) * (*eigenvalue_iter)[m]));

        vector_mul_vector_op<T, Device>()(this->ctx,
                                          this->dim,
                                          psi_iter + (nbase + m) * this->dim,
                                          psi_iter + (nbase + m) * this->dim,
                                          e_temp_cpu.data());
    }

    // psi_iter[:, nbase:nbase+notconv] += hphi[:, :nbase] * vcc[:nbase, :notconv]
    //
    // in terms of input, psi_iter should be
    // psi_iter[:, nbase:nbase+notconv] :=
    // hpsi[:, :nbase] * vcc[:nbase, :notconv] - e[]*psi_iter[:, nbase:nbase+notconv]*vcc[:nbase, :notconv]
    // TODO two gemm operations can be combined into one
    gemm_op<T, Device>()(this->ctx,
                         'N',
                         'N',
                         this->dim,
                         notconv,
                         nbase,
                         this->one,
                         hphi,
                         this->dim,
                         vcc,
                         this->nbase_x,
                         this->one,
                         psi_iter + (nbase) * this->dim,
                         this->dim);

    // FIXME 
    // (QE dav solver -> isolver = 0)
    // In QE, for davidson, h_diag is computed by g2kin + v_of_0
    // this is passed to g_psi, where the residual vector is element-wise
    // multiplied by (h_diag - e * s_siag)
    //
    // There is a macro "TEST_NEW_PRECONDITIONING" in QE, which
    // takes x = h_diag - e * s_siag and define
    // denm = 0.5 * ( 1 + x + sqrt(1 + (x-1)^2) )
    // which then precondition the residual vector by psi := psi ./ denm
    // note that this preconditioning recovers x when x is large, and
    // approaches 0.5*(1+sqrt(2)) when x is close to 0.
    // 
    //
    // In the current version of ABACUS, "precondition" is simply g2kin,
    // and the tested preconditioning of QE is always active.
    //
    // "precondition!!!"
    std::vector<Real> pre(this->dim, 0.0);
    for (int m = 0; m < notconv; m++)
    {
        for (size_t i = 0; i < this->dim; i++)
        {
            //double x = this->precondition[i] - (*eigenvalue_iter)[m];
            double x = std::abs(this->precondition[i] - (*eigenvalue_iter)[m]);
            pre[i] = 0.5 * (1.0 + x + sqrt(1 + (x - 1.0) * (x - 1.0)));
        }
        vector_div_vector_op<T, Device>()(this->ctx,
                                          this->dim,
                                          psi_iter + (nbase + m) * this->dim,
                                          psi_iter + (nbase + m) * this->dim,
                                          pre.data());
    }

    // "normalize!!!" in order to improve numerical stability of subspace diagonalization

    // column-wise normalize psi_iter[:, nbase:nbase+notconv]
    std::vector<Real> psi_norm(notconv, 0.0);
    for (size_t i = 0; i < notconv; i++)
    {
        psi_norm[i] = dot_real_op<T, Device>()(this->ctx,
                                               this->dim,
                                               psi_iter + (nbase + i) * this->dim,
                                               psi_iter + (nbase + i) * this->dim,
                                               false); // FIXME why reduce is false?
        assert(psi_norm[i] > 0.0);
        psi_norm[i] = sqrt(psi_norm[i]);

        vector_div_constant_op<T, Device>()(this->ctx,
                                            this->dim,
                                            psi_iter + (nbase + i) * this->dim,
                                            psi_iter + (nbase + i) * this->dim,
                                            psi_norm[i]);
    }

    // update hpsi[:, nbase:nbase+notconv]
    // hpsi[:, nbase:nbase+notconv] = H * psi_iter[:, nbase:nbase+notconv]
    hpsi_func(psi_iter + nbase * dim, hphi + nbase * this->dim, this->dim, notconv);

    ModuleBase::timer::tick("Diago_DavSubspace", "cal_grad");
    return;
}

template <typename T, typename Device>
void Diago_DavSubspace<T, Device>::cal_elem(const int& dim,
                                            int& nbase,
                                            const int& notconv,
                                            const T* psi_iter,
                                            const T* hphi,
                                            T* hcc,
                                            T* scc)
{
    ModuleBase::timer::tick("Diago_DavSubspace", "cal_elem");

    // hcc[:, nbase:] = psi_iter.H * hpsi_iter[:, nbase:]
    // while psi_iter has nbase_x columns, only the first nbase + notconv columns are
    // involved in the matrix multiplication
    gemm_op<T, Device>()(this->ctx,
                         'C',
                         'N',
                         nbase + notconv,
                         notconv,
                         this->dim,
                         this->one,
                         psi_iter,
                         this->dim,
                         &hphi[nbase * this->dim],
                         this->dim,
                         this->zero,
                         &hcc[nbase * this->nbase_x],
                         this->nbase_x);

    // scc[:, :nbase] := psi_iter.H * psi_iter[:, nbase:]
    gemm_op<T, Device>()(this->ctx,
                         'C',
                         'N',
                         nbase + notconv,
                         notconv,
                         this->dim,
                         this->one,
                         psi_iter,
                         this->dim,
                         psi_iter + nbase * this->dim,
                         this->dim,
                         this->zero,
                         &scc[nbase * this->nbase_x],
                         this->nbase_x);

#ifdef __MPI
    if (this->diag_comm.nproc > 1)
    {
        auto* swap = new T[notconv * this->nbase_x];
        syncmem_complex_op()(this->ctx, this->ctx, swap, hcc + nbase * this->nbase_x, notconv * this->nbase_x);

        if (std::is_same<T, double>::value)
        {
            Parallel_Reduce::reduce_pool(hcc + nbase * this->nbase_x, notconv * this->nbase_x);
            Parallel_Reduce::reduce_pool(scc + nbase * this->nbase_x, notconv * this->nbase_x);
        }
        else
        {
            if (base_device::get_current_precision(swap) == "single")
            {
                MPI_Reduce(swap,
                           hcc + nbase * this->nbase_x,
                           notconv * this->nbase_x,
                           MPI_COMPLEX,
                           MPI_SUM,
                           0,
                           this->diag_comm.comm);
            }
            else
            {
                MPI_Reduce(swap,
                           hcc + nbase * this->nbase_x,
                           notconv * this->nbase_x,
                           MPI_DOUBLE_COMPLEX,
                           MPI_SUM,
                           0,
                           this->diag_comm.comm);
            }

            syncmem_complex_op()(this->ctx, this->ctx, swap, scc + nbase * this->nbase_x, notconv * this->nbase_x);

            if (base_device::get_current_precision(swap) == "single")
            {
                MPI_Reduce(swap,
                           scc + nbase * this->nbase_x,
                           notconv * this->nbase_x,
                           MPI_COMPLEX,
                           MPI_SUM,
                           0,
                           this->diag_comm.comm);
            }
            else
            {
                MPI_Reduce(swap,
                           scc + nbase * this->nbase_x,
                           notconv * this->nbase_x,
                           MPI_DOUBLE_COMPLEX,
                           MPI_SUM,
                           0,
                           this->diag_comm.comm);
            }
        }
        delete[] swap;
    }
#endif

    const size_t last_nbase = nbase; // init: last_nbase = 0
    nbase = nbase + notconv;


    // NOTE: nbase is an in & out parameter! nbase is updated to nbase + notconv

    for (size_t i = 0; i < nbase; i++)
    {
        if (i >= last_nbase)
        {
            // ensure the diagonal elements of hcc and scc are real
            // NOTE: set_real_tocomplex convert a real to itself, and a complex to a new complex
            // whose real part is the input's real part and imaginary part is zero
            hcc[i * this->nbase_x + i] = set_real_tocomplex(hcc[i * this->nbase_x + i]);
            scc[i * this->nbase_x + i] = set_real_tocomplex(scc[i * this->nbase_x + i]);
        }
        for (size_t j = std::max(i + 1, last_nbase); j < nbase; j++)
        {
            // ensure the hermicity of hcc and scc
            hcc[i * this->nbase_x + j] = get_conj(hcc[j * this->nbase_x + i]);
            scc[i * this->nbase_x + j] = get_conj(scc[j * this->nbase_x + i]);
        }
    }

    // make hcc[nbase:, nbase:] and scc[nbase:, nbase:] to be zero
    for (size_t i = nbase; i < this->nbase_x; i++)
    {
        for (size_t j = nbase; j < this->nbase_x; j++)
        {
            hcc[i * this->nbase_x + j] = cs.zero;
            scc[i * this->nbase_x + j] = cs.zero;
            hcc[j * this->nbase_x + i] = cs.zero;
            scc[j * this->nbase_x + i] = cs.zero;
        }
    }

    ModuleBase::timer::tick("Diago_DavSubspace", "cal_elem");
    return;
}

template <typename T, typename Device>
void Diago_DavSubspace<T, Device>::diag_zhegvx(const int& nbase,
                                               const int& nband,
                                               T* hcc,
                                               T* scc,
                                               const int& nbase_x,
                                               std::vector<Real>* eigenvalue_iter,
                                               T* vcc,
                                               bool init,
                                               bool is_subspace)
{
    ModuleBase::timer::tick("Diago_DavSubspace", "diag_zhegvx");

    if (is_subspace == false)
    {
        if (this->diag_comm.rank == 0)
        {
            assert(nbase_x >= std::max(1, nbase));

            std::vector<std::vector<T>> h_diag(nbase, std::vector<T>(nbase, cs.zero));
            std::vector<std::vector<T>> s_diag(nbase, std::vector<T>(nbase, cs.zero));

            for (size_t i = 0; i < nbase; i++)
            {
                for (size_t j = 0; j < nbase; j++)
                {
                    h_diag[i][j] = hcc[i * this->nbase_x + j];
                    s_diag[i][j] = scc[i * this->nbase_x + j];
                }
            }

            if (this->device == base_device::GpuDevice)
            {
#if defined(__CUDA) || defined(__ROCM)
                Real* eigenvalue_gpu = nullptr;
                resmem_real_op()(this->ctx, eigenvalue_gpu, this->nbase_x);

                syncmem_var_h2d_op()(this->ctx,
                                     this->cpu_ctx,
                                     eigenvalue_gpu,
                                     (*eigenvalue_iter).data(),
                                     this->nbase_x);

                dnevx_op<T, Device>()(this->ctx, nbase, this->nbase_x, this->hcc, nband, eigenvalue_gpu, this->vcc);

                syncmem_var_d2h_op()(this->cpu_ctx,
                                     this->ctx,
                                     (*eigenvalue_iter).data(),
                                     eigenvalue_gpu,
                                     this->nbase_x);

                delmem_real_op()(this->ctx, eigenvalue_gpu);
#endif
            }
            else
            {
                if (init)
                {
                    dnevx_op<T, Device>()(this->ctx,
                                          nbase,
                                          this->nbase_x,
                                          this->hcc,
                                          nband,
                                          (*eigenvalue_iter).data(),
                                          this->vcc);
                }
                else
                {

                    dngvx_op<T, Device>()(this->ctx,
                                          nbase,
                                          this->nbase_x,
                                          this->hcc,
                                          this->scc,
                                          nband,
                                          (*eigenvalue_iter).data(),
                                          this->vcc);
                }
            }

            // reset:
            for (size_t i = 0; i < nbase; i++)
            {
                for (size_t j = 0; j < nbase; j++)
                {
                    hcc[i * this->nbase_x + j] = h_diag[i][j];
                    scc[i * this->nbase_x + j] = s_diag[i][j];
                }

                for (size_t j = nbase; j < this->nbase_x; j++)
                {
                    hcc[i * this->nbase_x + j] = cs.zero;
                    hcc[j * this->nbase_x + i] = cs.zero;
                    scc[i * this->nbase_x + j] = cs.zero;
                    scc[j * this->nbase_x + i] = cs.zero;
                }
            }
        }

#ifdef __MPI
        if (this->diag_comm.nproc > 1)
        {
            // vcc: nbase * nband
            for (int i = 0; i < nband; i++)
            {
                MPI_Bcast(&vcc[i * this->nbase_x], nbase, MPI_DOUBLE_COMPLEX, 0, this->diag_comm.comm);
            }
            MPI_Bcast((*eigenvalue_iter).data(), nband, MPI_DOUBLE, 0, this->diag_comm.comm);
        }
#endif
    }
    else if (is_subspace == true)
    {
        for (size_t m = 0; m < nband; m++)
        {
            (*eigenvalue_iter)[m] = get_real(hcc[m * this->nbase_x + m]);

            vcc[m * this->nbase_x + m] = set_real_tocomplex(1.0);
        }

#ifdef __MPI
        if (this->diag_comm.nproc > 1)
        {
            MPI_Bcast((*eigenvalue_iter).data(), this->n_band, MPI_DOUBLE, 0, this->diag_comm.comm);
        }
#endif
    }

    ModuleBase::timer::tick("Diago_DavSubspace", "diag_zhegvx");
    return;
}

template <typename T, typename Device>
void Diago_DavSubspace<T, Device>::refresh(const int& dim,
                                           const int& nband,
                                           int& nbase,
                                           const Real* eigenvalue_in_hsolver,
                                           //    const psi::Psi<T, Device>& psi,
                                           T* psi_iter,
                                           T* hp,
                                           T* sp,
                                           T* hc,
                                           T* vc)
{
    ModuleBase::timer::tick("Diago_DavSubspace", "refresh");

    gemm_op<T, Device>()(this->ctx,
                         'N',
                         'N',
                         this->dim,
                         nband,
                         nbase,
                         this->one,
                         this->hphi,
                         this->dim,
                         this->vcc,
                         this->nbase_x,
                         this->zero,
                         psi_iter + nband * this->dim,
                         this->dim);

    // update hphi
    syncmem_complex_op()(this->ctx, this->ctx, hphi, psi_iter + nband * this->dim, this->dim * nband);

    nbase = nband;

    // set hcc/scc/vcc to 0
    for (size_t i = 0; i < nbase; i++)
    {
        setmem_complex_op()(this->ctx, &hcc[this->nbase_x * i], 0, nbase);
        setmem_complex_op()(this->ctx, &scc[this->nbase_x * i], 0, nbase);
        setmem_complex_op()(this->ctx, &vcc[this->nbase_x * i], 0, nbase);
    }

    if (this->device == base_device::GpuDevice)
    {
#if defined(__CUDA) || defined(__ROCM)
        T* hcc_cpu = nullptr;
        T* scc_cpu = nullptr;
        T* vcc_cpu = nullptr;
        base_device::memory::resize_memory_op<T, base_device::DEVICE_CPU>()(this->cpu_ctx,
                                                                            hcc_cpu,
                                                                            this->nbase_x * this->nbase_x,
                                                                            "DAV::hcc");
        base_device::memory::resize_memory_op<T, base_device::DEVICE_CPU>()(this->cpu_ctx,
                                                                            scc_cpu,
                                                                            this->nbase_x * this->nbase_x,
                                                                            "DAV::scc");
        base_device::memory::resize_memory_op<T, base_device::DEVICE_CPU>()(this->cpu_ctx,
                                                                            vcc_cpu,
                                                                            this->nbase_x * this->nbase_x,
                                                                            "DAV::vcc");

        syncmem_d2h_op()(this->cpu_ctx, this->ctx, hcc_cpu, hcc, this->nbase_x * this->nbase_x);
        syncmem_d2h_op()(this->cpu_ctx, this->ctx, scc_cpu, scc, this->nbase_x * this->nbase_x);
        syncmem_d2h_op()(this->cpu_ctx, this->ctx, vcc_cpu, vcc, this->nbase_x * this->nbase_x);

        for (int i = 0; i < nbase; i++)
        {
            hcc_cpu[i * this->nbase_x + i] = eigenvalue_in_hsolver[i];
            scc_cpu[i * this->nbase_x + i] = this->one[0];
            vcc_cpu[i * this->nbase_x + i] = this->one[0];
        }

        syncmem_h2d_op()(this->ctx, this->cpu_ctx, hcc, hcc_cpu, this->nbase_x * this->nbase_x);
        syncmem_h2d_op()(this->ctx, this->cpu_ctx, scc, scc_cpu, this->nbase_x * this->nbase_x);
        syncmem_h2d_op()(this->ctx, this->cpu_ctx, vcc, vcc_cpu, this->nbase_x * this->nbase_x);

        base_device::memory::delete_memory_op<T, base_device::DEVICE_CPU>()(this->cpu_ctx, hcc_cpu);
        base_device::memory::delete_memory_op<T, base_device::DEVICE_CPU>()(this->cpu_ctx, scc_cpu);
        base_device::memory::delete_memory_op<T, base_device::DEVICE_CPU>()(this->cpu_ctx, vcc_cpu);
#endif
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
    ModuleBase::timer::tick("Diago_DavSubspace", "refresh");

    return;
}

template <typename T, typename Device>
int Diago_DavSubspace<T, Device>::diag(const HPsiFunc& hpsi_func,
                                       T* psi_in,
                                       const int psi_in_dmax,
                                       Real* eigenvalue_in_hsolver,
                                       const std::vector<bool>& is_occupied,
                                       const bool& scf_type)
{
    /// record the times of trying iterative diagonalization
    this->notconv = 0;

    int sum_iter = 0;
    int ntry = 0;

    do
    {
        //printf("enter diag... is_subspace = %d, ntry = %d\n", this->is_subspace, ntry);
        if (this->is_subspace || ntry > 0)
        {
            this->diagH_subspace(psi_in, eigenvalue_in_hsolver, hpsi_func, this->n_band, this->dim, psi_in_dmax);
        }

        sum_iter += this->diag_once(hpsi_func, psi_in, psi_in_dmax, eigenvalue_in_hsolver, is_occupied);

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

    // if self consistent calculation, if not converged > 5,
    // using diagH_subspace and cg method again. ntry++
    const bool f3 = ((scf && (notconv > 5)));

    return (f1 && (f2 || f3));
}

template <typename T, typename Device>
void Diago_DavSubspace<T, Device>::diagH_subspace(T* psi_pointer, // [in] & [out] wavefunction
                                                  Real* en,       // [out] eigenvalues
                                                  const HPsiFunc hpsi_func,
                                                  const int n_band,
                                                  const int dmin,
                                                  const int dmax)
{
    ModuleBase::TITLE("Diago_DavSubspace", "subspace");
    ModuleBase::timer::tick("Diago_DavSubspace", "subspace");

    const int nstart = n_band;
    assert(n_band > 0);

    T* hcc = nullptr;
    T* scc = nullptr;
    T* vcc = nullptr;
    resmem_complex_op()(ctx, hcc, nstart * nstart, "DAV::hcc");
    resmem_complex_op()(ctx, scc, nstart * nstart, "DAV::scc");
    resmem_complex_op()(ctx, vcc, nstart * nstart, "DAV::vcc");
    setmem_complex_op()(ctx, hcc, 0, nstart * nstart);
    setmem_complex_op()(ctx, scc, 0, nstart * nstart);
    setmem_complex_op()(ctx, vcc, 0, nstart * nstart);

    T* hphi = nullptr;
    resmem_complex_op()(ctx, hphi, nstart * dmax, "DAV::hphi");
    setmem_complex_op()(ctx, hphi, 0, nstart * dmax);

    {
        // do hPsi for all bands
        // hphi[:, 0:nstart] = H * psi_pointer[:, 0:nstart]
        hpsi_func(psi_pointer, hphi, dmax, nstart);

        gemm_op<T, Device>()(ctx,
                             'C',
                             'N',
                             nstart,
                             nstart,
                             dmin,
                             this->one,
                             psi_pointer,
                             dmax,
                             hphi,
                             dmax,
                             this->zero,
                             hcc,
                             nstart);

        gemm_op<T, Device>()(ctx,
                             'C',
                             'N',
                             nstart,
                             nstart,
                             dmin,
                             this->one,
                             psi_pointer,
                             dmax,
                             psi_pointer,
                             dmax,
                             this->zero,
                             scc,
                             nstart);
    }

    if (GlobalV::NPROC_IN_POOL > 1)
    {
        Parallel_Reduce::reduce_pool(hcc, nstart * nstart);
        Parallel_Reduce::reduce_pool(scc, nstart * nstart);
    }

    // after generation of H and S matrix, diag them
    DiagoIterAssist<T, Device>::diagH_LAPACK(nstart, n_band, hcc, scc, nstart, en, vcc);

    { // code block to calculate evc
        gemm_op<T, Device>()(ctx,
                             'N',
                             'N',
                             dmin,
                             n_band,
                             nstart,
                             this->one,
                             psi_pointer, // dmin * nstart
                             dmax,
                             vcc, // nstart * n_band
                             nstart,
                             this->zero,
                             hphi,
                             dmin);
    }

    matrixSetToAnother<T, Device>()(ctx, n_band, hphi, dmin, psi_pointer, dmax);

    delmem_complex_op()(ctx, hphi);

    delmem_complex_op()(ctx, hcc);
    delmem_complex_op()(ctx, scc);
    delmem_complex_op()(ctx, vcc);

    ModuleBase::timer::tick("Diago_DavSubspace", "diagH_subspace");
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
