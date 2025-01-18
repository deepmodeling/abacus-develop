#include "diago_dav_subspace.h"

#include "diago_iter_assist.h"
#include "module_base/memory.h"
#include "module_base/module_device/device.h"
#include "module_base/timer.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_base/kernels/dsp/dsp_connector.h"

#include "module_hsolver/diag_hs_para.h"

#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

using namespace hsolver;

template <typename T, typename Device>
Diago_DavSubspace<T, Device>::Diago_DavSubspace(const std::vector<Real>& precondition_in,
                                                const int& nband_in,
                                                const int& nbasis_in,
                                                const int& david_ndim_in,
                                                const double& diag_thr_in,
                                                const int& diag_nmax_in,
                                                const bool& need_subspace_in,
                                                const diag_comm_info& diag_comm_in,
                                                const int diag_subspace_in,
                                                const int diago_subspace_bs_in)
    : precondition(precondition_in), n_band(nband_in), dim(nbasis_in), nbase_x(nband_in * david_ndim_in),
      diag_thr(diag_thr_in), iter_nmax(diag_nmax_in), is_subspace(need_subspace_in), diag_comm(diag_comm_in),
        diag_subspace(diag_subspace_in), diago_subspace_bs(diago_subspace_bs_in)
{
    this->device = base_device::get_device_type<Device>(this->ctx);

    this->one = &one_;
    this->zero = &zero_;
    this->neg_one = &neg_one_;

    assert(david_ndim_in > 1);
    assert(david_ndim_in * nband_in < nbasis_in * this->diag_comm.nproc);
    assert(diag_subspace >= 0 && diag_subspace < 3);

    // TODO: Added memory usage statistics

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    resmem_complex_op()(this->psi_in_iter, this->nbase_x * this->dim, "DAV::psi_in_iter");
    setmem_complex_op()(this->psi_in_iter, 0, this->nbase_x * this->dim);

    // the product of H and psi in the reduced psi set
    resmem_complex_op()(this->hphi, this->nbase_x * this->dim, "DAV::hphi");
    setmem_complex_op()(this->hphi, 0, this->nbase_x * this->dim);

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
        // syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, this->d_precondition, this->precondition.data(), nbasis_in);
    }
#endif
}

template <typename T, typename Device>
Diago_DavSubspace<T, Device>::~Diago_DavSubspace()
{
    delmem_complex_op()(this->psi_in_iter);

    delmem_complex_op()(this->hphi);
    delmem_complex_op()(this->hcc);
    delmem_complex_op()(this->scc);
    delmem_complex_op()(this->vcc);

#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        delmem_real_op()(this->d_precondition);
    }
#endif
}

template <typename T, typename Device>
int Diago_DavSubspace<T, Device>::diag_once(const HPsiFunc& hpsi_func,
                                            T* psi_in,
                                            const int psi_in_dmax,
                                            Real* eigenvalue_in_hsolver,
                                            const std::vector<double>& ethr_band)
{
    ModuleBase::timer::tick("Diago_DavSubspace", "diag_once");

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

    ModuleBase::timer::tick("Diago_DavSubspace", "first");

    for (int m = 0; m < this->n_band; m++)
    {
        unconv[m] = m;

        syncmem_complex_op()(this->psi_in_iter + m * this->dim,
                             psi_in + m * psi_in_dmax,
                             this->dim);
    }

    // compute h*psi_in_iter
    // NOTE: bands after the first n_band should yield zero
    // hphi[:, 0:nbase_x] = H * psi_in_iter[:, 0:nbase_x]
    hpsi_func(this->psi_in_iter, this->hphi, this->dim, this->notconv);

    this->cal_elem(this->dim, nbase, this->notconv, this->psi_in_iter, this->hphi, this->hcc, this->scc);

    this->diag_zhegvx(nbase, this->notconv, this->hcc, this->scc, this->nbase_x, &eigenvalue_iter, this->vcc);

    for (size_t m = 0; m < this->n_band; m++)
    {
        eigenvalue_in_hsolver[m] = eigenvalue_iter[m];
    }

    ModuleBase::timer::tick("Diago_DavSubspace", "first");

    int dav_iter = 0;

    do
    {
        dav_iter++;

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

        this->diag_zhegvx(nbase, this->n_band, this->hcc, this->scc, this->nbase_x, &eigenvalue_iter, this->vcc);

        // check convergence and update eigenvalues
        ModuleBase::timer::tick("Diago_DavSubspace", "check_update");

        this->notconv = 0;
        for (int m = 0; m < this->n_band; m++)
        {
            convflag[m] = (std::abs(eigenvalue_iter[m] - eigenvalue_in_hsolver[m]) < ethr_band[m]);

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
            setmem_complex_op()(psi_in, 0, n_band * psi_in_dmax);

#ifdef __DSP
    gemm_op_mt<T, Device>()  // In order to not coding another whole template, using this method to minimize the code change.
#else
    gemm_op<T, Device>()
#endif
                                (this->ctx,
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
                    syncmem_complex_op()(this->psi_in_iter + i * this->dim,
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
            syncmem_complex_op()(vcc + i * this->nbase_x, vcc + unconv[i] * this->nbase_x, nbase);
            (*eigenvalue_iter)[i] = (*eigenvalue_iter)[unconv[i]];
        }
    }

#ifdef __DSP
    gemm_op_mt<T, Device>()
#else
    gemm_op<T, Device>()
#endif
                        (this->ctx,
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
                         this->zero,
                         psi_iter + (nbase) * this->dim,
                         this->dim);

    std::vector<Real> e_temp_cpu(nbase, 0);
    Real* e_temp_hd = e_temp_cpu.data();
    if(this->device == base_device::GpuDevice)
    {
        e_temp_hd = nullptr;
        resmem_real_op()(e_temp_hd, nbase);
    }
    for (int m = 0; m < notconv; m++)
    {
        e_temp_cpu.assign(nbase, (-1.0 * (*eigenvalue_iter)[m]));
        if (this->device == base_device::GpuDevice)
        {
            syncmem_var_h2d_op()(e_temp_hd, e_temp_cpu.data(), nbase);
        }
        vector_mul_vector_op<T, Device>()(this->ctx,
                                                nbase,
                                                vcc + m * this->nbase_x,
                                                vcc + m * this->nbase_x,
                                                e_temp_hd);
    }
    if(this->device == base_device::GpuDevice)
    {
        delmem_real_op()(e_temp_hd);
    }

#ifdef __DSP
    gemm_op_mt<T, Device>()
#else
    gemm_op<T, Device>()
#endif                  
                        (this->ctx,
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
                         this->one,
                         psi_iter + nbase * this->dim,
                         this->dim);

    // "precondition!!!"
    std::vector<Real> pre(this->dim, 0.0);
    for (int m = 0; m < notconv; m++)
    {
        for (size_t i = 0; i < this->dim; i++)
        {
            // pre[i] = std::abs(this->precondition[i] - (*eigenvalue_iter)[m]);
            double x = std::abs(this->precondition[i] - (*eigenvalue_iter)[m]);
            pre[i] = 0.5 * (1.0 + x + sqrt(1 + (x - 1.0) * (x - 1.0)));
        }
#if defined(__CUDA) || defined(__ROCM)
        if (this->device == base_device::GpuDevice)
        {
            syncmem_var_h2d_op()(this->d_precondition, pre.data(), this->dim);
            vector_div_vector_op<T, Device>()(this->ctx,
                                              this->dim,
                                              psi_iter + (nbase + m) * this->dim,
                                              psi_iter + (nbase + m) * this->dim,
                                              this->d_precondition);
        }
        else
#endif
        {
            vector_div_vector_op<T, Device>()(this->ctx,
                                              this->dim,
                                              psi_iter + (nbase + m) * this->dim,
                                              psi_iter + (nbase + m) * this->dim,
                                              pre.data());
        }
    }

    // "normalize!!!" in order to improve numerical stability of subspace diagonalization
    std::vector<Real> psi_norm(notconv, 0.0);
    for (size_t i = 0; i < notconv; i++)
    {
        psi_norm[i] = dot_real_op<T, Device>()(this->ctx,
                                               this->dim,
                                               psi_iter + (nbase + i) * this->dim,
                                               psi_iter + (nbase + i) * this->dim,
                                               true);
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

#ifdef __DSP
    gemm_op_mt<T, Device>()
#else
    gemm_op<T, Device>()
#endif 
                        (this->ctx,
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

#ifdef __DSP
    gemm_op_mt<T, Device>()
#else
    gemm_op<T, Device>()
#endif
                        (this->ctx,
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
#ifdef __DSP
        // Only on dsp hardware need an extra space to reduce data
        dsp_dav_subspace_reduce(hcc, scc, nbase, this->nbase_x, this->notconv, this->diag_comm.comm);
#else
        auto* swap = new T[notconv * this->nbase_x];

        syncmem_complex_op()(swap, hcc + nbase * this->nbase_x, notconv * this->nbase_x);

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

            syncmem_complex_op()(swap, scc + nbase * this->nbase_x, notconv * this->nbase_x);

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
#endif
    }
#endif

    const size_t last_nbase = nbase; // init: last_nbase = 0
    nbase = nbase + notconv;

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
                                               T* vcc)
{
    ModuleBase::timer::tick("Diago_DavSubspace", "diag_zhegvx");
    assert(nbase_x >= std::max(1, nbase));

    if (this->device == base_device::GpuDevice)
    {
#if defined(__CUDA) || defined(__ROCM)
        if (this->diag_comm.rank == 0)
        {
            Real* eigenvalue_gpu = nullptr;
            resmem_real_op()(eigenvalue_gpu, this->nbase_x);

            syncmem_var_h2d_op()(eigenvalue_gpu, (*eigenvalue_iter).data(), this->nbase_x);

            T* hcc_gpu = nullptr;
            T* scc_gpu = nullptr;
            T* vcc_gpu = nullptr;
            base_device::memory::resize_memory_op<T, Device>()(hcc_gpu, nbase * nbase);
            base_device::memory::resize_memory_op<T, Device>()(scc_gpu, nbase * nbase);
            base_device::memory::resize_memory_op<T, Device>()(vcc_gpu, nbase * nbase);
            for(int i=0;i<nbase;i++)
            {
                base_device::memory::synchronize_memory_op<T, Device, Device>()(hcc_gpu + i * nbase, hcc + i * nbase_x, nbase);
                base_device::memory::synchronize_memory_op<T, Device, Device>()(scc_gpu + i * nbase, scc + i * nbase_x, nbase);
            }
            dngvd_op<T, Device>()(this->ctx, nbase, nbase, hcc_gpu, scc_gpu, eigenvalue_gpu, vcc_gpu);
            for(int i=0;i<nbase;i++)
            {
                base_device::memory::synchronize_memory_op<T, Device, Device>()(vcc + i * nbase_x, vcc_gpu + i * nbase, nbase);
            }
            delmem_complex_op()(hcc_gpu);
            delmem_complex_op()(scc_gpu);
            delmem_complex_op()(vcc_gpu);

            syncmem_var_d2h_op()((*eigenvalue_iter).data(), eigenvalue_gpu, this->nbase_x);

            delmem_real_op()(eigenvalue_gpu);
        }
#endif
    }
    else
    {
        if (this->diag_subspace == 0)
        {
            if (this->diag_comm.rank == 0)
            {
                std::vector<std::vector<T>> h_diag(nbase, std::vector<T>(nbase, *this->zero));
                std::vector<std::vector<T>> s_diag(nbase, std::vector<T>(nbase, *this->zero));

                for (size_t i = 0; i < nbase; i++)
                {
                    for (size_t j = 0; j < nbase; j++)
                    {
                        h_diag[i][j] = hcc[i * this->nbase_x + j];
                        s_diag[i][j] = scc[i * this->nbase_x + j];
                    }
                }
                dngvx_op<T, Device>()(this->ctx,
                                      nbase,
                                      this->nbase_x,
                                      this->hcc,
                                      this->scc,
                                      nband,
                                      (*eigenvalue_iter).data(),
                                      this->vcc);
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
                        hcc[i * this->nbase_x + j] = *this->zero;
                        hcc[j * this->nbase_x + i] = *this->zero;
                        scc[i * this->nbase_x + j] = *this->zero;
                        scc[j * this->nbase_x + i] = *this->zero;
                    }
                }
            }
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
            MPI_Bcast(&vcc[i * this->nbase_x], nbase, MPI_DOUBLE_COMPLEX, 0, this->diag_comm.comm);
        }
        MPI_Bcast((*eigenvalue_iter).data(), nband, MPI_DOUBLE, 0, this->diag_comm.comm);
    }
#endif

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

#ifdef __DSP
    gemm_op_mt<T, Device>()
#else
    gemm_op<T, Device>()
#endif
                        (this->ctx,
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
    syncmem_complex_op()(hphi, psi_iter + nband * this->dim, this->dim * nband);

    nbase = nband;

    // set hcc/scc/vcc to 0
    for (size_t i = 0; i < nbase; i++)
    {
        setmem_complex_op()(&hcc[this->nbase_x * i], 0, nbase);
        setmem_complex_op()(&scc[this->nbase_x * i], 0, nbase);
        setmem_complex_op()(&vcc[this->nbase_x * i], 0, nbase);
    }

    if (this->device == base_device::GpuDevice)
    {
#if defined(__CUDA) || defined(__ROCM)
        T* hcc_cpu = nullptr;
        T* scc_cpu = nullptr;
        T* vcc_cpu = nullptr;
        base_device::memory::resize_memory_op<T, base_device::DEVICE_CPU>()(hcc_cpu,
                                                                            this->nbase_x * this->nbase_x,
                                                                            "DAV::hcc");
        base_device::memory::resize_memory_op<T, base_device::DEVICE_CPU>()(scc_cpu,
                                                                            this->nbase_x * this->nbase_x,
                                                                            "DAV::scc");
        base_device::memory::resize_memory_op<T, base_device::DEVICE_CPU>()(vcc_cpu,
                                                                            this->nbase_x * this->nbase_x,
                                                                            "DAV::vcc");

        syncmem_d2h_op()(hcc_cpu, hcc, this->nbase_x * this->nbase_x);
        syncmem_d2h_op()(scc_cpu, scc, this->nbase_x * this->nbase_x);
        syncmem_d2h_op()(vcc_cpu, vcc, this->nbase_x * this->nbase_x);

        for (int i = 0; i < nbase; i++)
        {
            hcc_cpu[i * this->nbase_x + i] = eigenvalue_in_hsolver[i];
            scc_cpu[i * this->nbase_x + i] = this->one[0];
            vcc_cpu[i * this->nbase_x + i] = this->one[0];
        }

        syncmem_h2d_op()(hcc, hcc_cpu, this->nbase_x * this->nbase_x);
        syncmem_h2d_op()(scc, scc_cpu, this->nbase_x * this->nbase_x);
        syncmem_h2d_op()(vcc, vcc_cpu, this->nbase_x * this->nbase_x);

        base_device::memory::delete_memory_op<T, base_device::DEVICE_CPU>()(hcc_cpu);
        base_device::memory::delete_memory_op<T, base_device::DEVICE_CPU>()(scc_cpu);
        base_device::memory::delete_memory_op<T, base_device::DEVICE_CPU>()(vcc_cpu);
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
                                       const std::vector<double>& ethr_band,
                                       const bool& scf_type)
{
    /// record the times of trying iterative diagonalization
    this->notconv = 0;

    int sum_iter = 0;
    int ntry = 0;

    do
    {

        sum_iter += this->diag_once(hpsi_func, psi_in, psi_in_dmax, eigenvalue_in_hsolver, ethr_band);

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
