#include "diago_dav_subspace.h"

#include "diago_iter_assist.h"
#include "module_base/memory.h"
#include "module_base/module_device/device.h"
#include "module_base/timer.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"

using namespace hsolver;

template <typename T, typename Device>
Diago_DavSubspace<T, Device>::Diago_DavSubspace(const Real* precondition_in,
                                                const int david_ndim_in,
                                                const double diag_thr_in,
                                                const int diag_nmax_in,
                                                const bool need_subspace_in,
                                                const diag_comm_info& diag_comm_in)
    : david_ndim_(david_ndim_in), diag_thr_(diag_thr_in), diag_nmax_(diag_nmax_in), is_subspace(need_subspace_in),
      diag_comm(diag_comm_in)
{
    this->device = base_device::get_device_type<Device>(this->ctx);
    this->precondition = precondition_in;

    this->one = &this->cs.one;
    this->zero = &this->cs.zero;
    this->neg_one = &this->cs.neg_one;
}

template <typename T, typename Device>
Diago_DavSubspace<T, Device>::~Diago_DavSubspace()
{
    delmem_complex_op()(this->ctx, this->hphi);
    delmem_complex_op()(this->ctx, this->hcc);
    delmem_complex_op()(this->ctx, this->scc);
    delmem_complex_op()(this->ctx, this->vcc);

    delmem_real_h_op()(this->cpu_ctx, this->eigenvalue_in_dav);

#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        delmem_real_op()(this->ctx, this->d_precondition);
    }
#endif
}

template <typename T, typename Device>
int Diago_DavSubspace<T, Device>::diag_once(hamilt::Hamilt<T, Device>* phm_in,
                                            psi::Psi<T, Device>& psi,
                                            Real* eigenvalue_in_hsolver,
                                            const std::vector<bool>& is_occupied)
{
    ModuleBase::timer::tick("Diago_DavSubspace", "diag_once");

    assert(this->david_ndim_ > 1);
    assert(this->david_ndim_ * psi.get_nbands() < psi.get_current_nbas() * this->diag_comm.nproc);

    this->dim = psi.get_k_first() ? psi.get_current_nbas() : psi.get_nk() * psi.get_nbasis();
    this->n_band = psi.get_nbands();

    // maximum dimension of the reduced basis set
    this->nbase_x = this->david_ndim_ * this->n_band;

    psi::Psi<T, Device> basis(1, this->nbase_x, this->dim, &(psi.get_ngk(0)));
    ModuleBase::Memory::record("DAV::basis", this->nbase_x * this->dim * sizeof(T));

    // the eigenvalues in dav iter
    resmem_real_h_op()(this->cpu_ctx, this->eigenvalue_in_dav, this->nbase_x, "DAV::eig");
    setmem_real_h_op()(this->cpu_ctx, this->eigenvalue_in_dav, 0, this->nbase_x);

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // the product of H and psi in the reduced basis set
    resmem_complex_op()(this->ctx, this->hphi, this->nbase_x * this->dim, "DAV::hphi");
    setmem_complex_op()(this->ctx, this->hphi, 0, this->nbase_x * this->dim);

    // Hamiltonian on the reduced basis set
    resmem_complex_op()(this->ctx, this->hcc, this->nbase_x * this->nbase_x, "DAV::hcc");
    setmem_complex_op()(this->ctx, this->hcc, 0, this->nbase_x * this->nbase_x);

    // Overlap on the reduced basis set
    resmem_complex_op()(this->ctx, this->scc, this->nbase_x * this->nbase_x, "DAV::scc");
    setmem_complex_op()(this->ctx, this->scc, 0, this->nbase_x * this->nbase_x);

    // Eigenvectors
    resmem_complex_op()(this->ctx, this->vcc, this->nbase_x * this->nbase_x, "DAV::vcc");
    setmem_complex_op()(this->ctx, this->vcc, 0, this->nbase_x * this->nbase_x);
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // convflag[m] = true if the m th band is convergent
    std::vector<bool> convflag(this->n_band, false);

    // unconv[m] store the number of the m th unconvergent band
    std::vector<int> unconv(this->n_band);

    // the dimension of the reduced basis set
    int nbase = 0;

    // the number of the unconvergent bands
    this->notconv = this->n_band;

    ModuleBase::timer::tick("Diago_DavSubspace", "first");

    for (int m = 0; m < this->n_band; m++)
    {
        unconv[m] = m;

        syncmem_complex_op()(this->ctx,
                             this->ctx,
                             &basis(m, 0),
                             psi.get_k_first() ? &psi(m, 0) : &psi(m, 0, 0),
                             this->dim);
    }

    // calculate H|psi>
    hpsi_info dav_hpsi_in(&basis, psi::Range(1, 0, 0, this->n_band - 1), this->hphi);
    phm_in->ops->hPsi(dav_hpsi_in);

    this->cal_elem(this->dim, nbase, this->notconv, basis, this->hphi, this->hcc, this->scc);

    this->diag_zhegvx(nbase,
                      this->n_band,
                      this->hcc,
                      this->scc,
                      this->nbase_x,
                      this->eigenvalue_in_dav,
                      this->vcc,
                      true,
                      this->is_subspace);

    for (size_t m = 0; m < this->n_band; m++)
    {
        eigenvalue_in_hsolver[m] = this->eigenvalue_in_dav[m];
    }

    ModuleBase::timer::tick("Diago_DavSubspace", "first");

    int dav_iter = 0;

    do
    {
        dav_iter++;

        this->cal_grad(phm_in,
                       this->dim,
                       nbase,
                       this->notconv,
                       basis,
                       this->hphi,
                       this->vcc,
                       unconv.data(),
                       this->eigenvalue_in_dav);

        this->cal_elem(this->dim, nbase, this->notconv, basis, this->hphi, this->hcc, this->scc);

        this->diag_zhegvx(nbase,
                          this->n_band,
                          this->hcc,
                          this->scc,
                          this->nbase_x,
                          this->eigenvalue_in_dav,
                          this->vcc,
                          false,
                          false);

        // check convergence and update eigenvalues
        ModuleBase::timer::tick("Diago_DavSubspace", "check_update");

        this->notconv = 0;
        for (int m = 0; m < this->n_band; m++)
        {
            if (is_occupied[m])
            {
                convflag[m] = (std::abs(this->eigenvalue_in_dav[m] - eigenvalue_in_hsolver[m]) < this->diag_thr_);
            }
            else
            {
                const double empty_ethr = std::max(this->diag_thr_ * 5.0, this->dav_large_thr);
                convflag[m] = (std::abs(this->eigenvalue_in_dav[m] - eigenvalue_in_hsolver[m]) < empty_ethr);
            }

            if (!convflag[m])
            {
                unconv[this->notconv] = m;
                this->notconv++;
            }

            eigenvalue_in_hsolver[m] = this->eigenvalue_in_dav[m];
        }

        ModuleBase::timer::tick("Diago_DavSubspace", "check_update");

        if ((this->notconv == 0) || (nbase + this->notconv + 1 > this->nbase_x) || (dav_iter == this->diag_nmax_))
        {
            ModuleBase::timer::tick("Diago_DavSubspace", "last");

            // updata eigenvectors of Hamiltonian
            setmem_complex_op()(this->ctx, psi.get_pointer(), 0, n_band * psi.get_nbasis());
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            // haozhihan repalce 2022-10-18
            gemm_op<T, Device>()(this->ctx,
                                 'N',
                                 'N',
                                 this->dim,    // m: row of A,C
                                 this->n_band, // n: col of B,C
                                 nbase,        // k: col of A, row of B
                                 this->one,
                                 basis.get_pointer(), // A dim * nbase
                                 this->dim,
                                 this->vcc, // B nbase * n_band
                                 this->nbase_x,
                                 this->zero,
                                 psi.get_pointer(), // C dim * n_band
                                 psi.get_nbasis());

            if (!this->notconv || (dav_iter == this->diag_nmax_))
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

                this->refresh(this->dim,
                              this->n_band,
                              nbase,
                              eigenvalue_in_hsolver,
                              psi,
                              basis,
                              this->hphi,
                              this->hcc,
                              this->scc,
                              this->vcc);
                ModuleBase::timer::tick("Diago_DavSubspace", "last");
            }
        }

    } while (1);

    ModuleBase::timer::tick("Diago_DavSubspace", "diag_once");

    return dav_iter;
}

template <typename T, typename Device>
void Diago_DavSubspace<T, Device>::cal_grad(hamilt::Hamilt<T, Device>* phm_in,
                                            const int& dim,
                                            const int& nbase,
                                            const int& notconv,
                                            psi::Psi<T, Device>& basis,
                                            T* hphi,
                                            T* vcc,
                                            const int* unconv,
                                            Real* eigenvalue_in_dav)
{
    ModuleBase::timer::tick("Diago_DavSubspace", "cal_grad");

    for (size_t i = 0; i < notconv; i++)
    {
        if (unconv[i] != i)
        {
            syncmem_complex_op()(this->ctx, this->ctx, vcc + i * this->nbase_x, vcc + unconv[i] * this->nbase_x, nbase);
            this->eigenvalue_in_dav[i] = this->eigenvalue_in_dav[unconv[i]];
        }
    }

    gemm_op<T, Device>()(this->ctx,
                         'N',
                         'N',
                         this->dim,        // m: row of A,C
                         notconv,          // n: col of B,C
                         nbase,            // k: col of A, row of B
                         this->one,        // alpha
                         &basis(0, 0),     // A
                         this->dim,        // LDA
                         vcc,              // B
                         this->nbase_x,    // LDB
                         this->zero,       // belta
                         &basis(nbase, 0), // C dim * notconv
                         this->dim         // LDC
    );

    for (int m = 0; m < notconv; m++)
    {

        std::vector<Real> e_temp_cpu(this->dim, (-this->eigenvalue_in_dav[m]));

        vector_mul_vector_op<T, Device>()(this->ctx,
                                          this->dim,
                                          &basis(nbase + m, 0),
                                          &basis(nbase + m, 0),
                                          e_temp_cpu.data());
    }

    gemm_op<T, Device>()(this->ctx,
                         'N',
                         'N',
                         this->dim,        // m: row of A,C
                         notconv,          // n: col of B,C
                         nbase,            // k: col of A, row of B
                         this->one,        // alpha
                         hphi,             // A dim * nbase
                         this->dim,        // LDA
                         vcc,              // B nbase * notconv
                         this->nbase_x,    // LDB
                         this->one,        // belta
                         &basis(nbase, 0), // C dim * notconv
                         this->dim         // LDC
    );

    // "precondition!!!"
    std::vector<Real> pre(this->dim, 0.0);
    for (int m = 0; m < notconv; m++)
    {
        for (size_t i = 0; i < this->dim; i++)
        {
            double x = this->precondition[i] - this->eigenvalue_in_dav[m];
            pre[i] = 0.5 * (1.0 + x + sqrt(1 + (x - 1.0) * (x - 1.0)));
        }
        vector_div_vector_op<T, Device>()(this->ctx, this->dim, &basis(nbase + m, 0), &basis(nbase + m, 0), pre.data());
    }

    // "normalize!!!" in order to improve numerical stability of subspace diagonalization
    std::vector<Real> psi_norm(notconv, 0.0);
    for (size_t i = 0; i < notconv; i++)
    {
        psi_norm[i] = dot_real_op<T, Device>()(this->ctx, this->dim, &basis(nbase + i, 0), &basis(nbase + i, 0), false);
        assert(psi_norm[i] > 0.0);
        psi_norm[i] = sqrt(psi_norm[i]);

        vector_div_constant_op<T, Device>()(this->ctx,
                                            this->dim,
                                            &basis(nbase + i, 0),
                                            &basis(nbase + i, 0),
                                            psi_norm[i]);
    }

    // "calculate H|psi>" for not convergence bands
    hpsi_info dav_hpsi_in(&basis, psi::Range(1, 0, nbase, nbase + notconv - 1), &hphi[nbase * this->dim]);
    phm_in->ops->hPsi(dav_hpsi_in);

    ModuleBase::timer::tick("Diago_DavSubspace", "cal_grad");
    return;
}

template <typename T, typename Device>
void Diago_DavSubspace<T, Device>::cal_elem(const int& dim,
                                            int& nbase,
                                            const int& notconv,
                                            const psi::Psi<T, Device>& basis,
                                            const T* hphi,
                                            T* hcc,
                                            T* scc)
{
    ModuleBase::timer::tick("Diago_DavSubspace", "cal_elem");

    gemm_op<T, Device>()(this->ctx,
                         'C',
                         'N',
                         nbase + notconv,
                         notconv,
                         this->dim,
                         this->one,
                         &basis(0, 0),
                         this->dim,
                         &hphi[nbase * this->dim],
                         this->dim,
                         this->zero,
                         &hcc[nbase * this->nbase_x],
                         this->nbase_x);

    gemm_op<T, Device>()(this->ctx,
                         'C',
                         'N',
                         nbase + notconv,
                         notconv,
                         this->dim,
                         this->one,
                         &basis(0, 0),
                         this->dim,
                         &basis(nbase, 0),
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

    for (size_t i = 0; i < nbase; i++)
    {
        if (i >= last_nbase)
        {
            hcc[i * this->nbase_x + i] = set_real_tocomplex(hcc[i * this->nbase_x + i]);
            scc[i * this->nbase_x + i] = set_real_tocomplex(scc[i * this->nbase_x + i]);
        }
        for (size_t j = std::max(i + 1, last_nbase); j < nbase; j++)
        {
            hcc[i * this->nbase_x + j] = get_conj(hcc[j * this->nbase_x + i]);
            scc[i * this->nbase_x + j] = get_conj(scc[j * this->nbase_x + i]);
        }
    }

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
                                               Real* eigenvalue_in_dav,
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

                syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, eigenvalue_gpu, this->eigenvalue_in_dav, this->nbase_x);

                dnevx_op<T, Device>()(this->ctx, nbase, this->nbase_x, this->hcc, nband, eigenvalue_gpu, this->vcc);

                syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, this->eigenvalue_in_dav, eigenvalue_gpu, this->nbase_x);

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
                                          this->eigenvalue_in_dav,
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
                                          this->eigenvalue_in_dav,
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
            MPI_Bcast(this->eigenvalue_in_dav, nband, MPI_DOUBLE, 0, this->diag_comm.comm);
        }
#endif
    }
    else if (is_subspace == true)
    {
        for (size_t m = 0; m < nband; m++)
        {
            this->eigenvalue_in_dav[m] = get_real(hcc[m * this->nbase_x + m]);

            vcc[m * this->nbase_x + m] = set_real_tocomplex(1.0);
        }

#ifdef __MPI
        if (this->diag_comm.nproc > 1)
        {
            MPI_Bcast(this->eigenvalue_in_dav, this->n_band, MPI_DOUBLE, 0, this->diag_comm.comm);
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
                                           const psi::Psi<T, Device>& psi,
                                           psi::Psi<T, Device>& basis,
                                           T* hp,
                                           T* sp,
                                           T* hc,
                                           T* vc)
{
    ModuleBase::timer::tick("Diago_DavSubspace", "refresh");

    // update basis
    for (size_t i = 0; i < nband; i++)
    {
        syncmem_complex_op()(this->ctx, this->ctx, &basis(i, 0), &psi(i, 0), this->dim);
    }
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
                         &basis(nband, 0),
                         this->dim);

    // update hphi
    syncmem_complex_op()(this->ctx, this->ctx, hphi, &basis(nband, 0), this->dim * nband);

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
int Diago_DavSubspace<T, Device>::diag(hamilt::Hamilt<T, Device>* phm_in,
                                       psi::Psi<T, Device>& psi,
                                       Real* eigenvalue_in_hsolver,
                                       const std::vector<bool>& is_occupied,
                                       const bool& scf_type)
{
    /// record the times of trying iterative diagonalization
    this->notconv = 0;

#if defined(__CUDA) || defined(__ROCM)
    if (this->device == base_device::GpuDevice)
    {
        resmem_real_op()(this->ctx, this->d_precondition, psi.get_nbasis());
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, this->d_precondition, this->precondition, psi.get_nbasis());
    }
#endif

    bool outputscc = false;
    bool outputeigenvalue = false;

    if (outputscc)
    {
        this->nbase_x = 2 * psi.get_nbands();
        resmem_complex_op()(this->ctx, this->scc, this->nbase_x * this->nbase_x, "DAV::scc");
        setmem_complex_op()(this->ctx, this->scc, 0, this->nbase_x * this->nbase_x);

        std::cout << "before dav 111" << std::endl;
        gemm_op<T, Device>()(this->ctx,
                             'C',
                             'N',
                             psi.get_nbands(),
                             psi.get_nbands(),
                             psi.get_current_nbas(),
                             this->one,
                             &psi(0, 0),
                             psi.get_current_nbas(),
                             &psi(0, 0),
                             psi.get_current_nbas(),
                             this->zero,
                             this->scc,
                             this->nbase_x);
        // output scc
        for (size_t i = 0; i < psi.get_nbands(); i++)
        {
            for (size_t j = 0; j < psi.get_nbands(); j++)
            {
                std::cout << this->scc[i * this->nbase_x + j] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    if (outputeigenvalue)
    {
        // output: eigenvalue_in_hsolver
        std::cout << "before dav, output eigenvalue_in_hsolver" << std::endl;
        for (size_t i = 0; i < psi.get_nbands(); i++)
        {
            std::cout << eigenvalue_in_hsolver[i] << "\t";
        }
        std::cout << std::endl;
    }

    int sum_iter = 0;
    int ntry = 0;
    do
    {
        if (this->is_subspace || ntry > 0)
        {
            DiagoIterAssist<T, Device>::diagH_subspace(phm_in, psi, psi, eigenvalue_in_hsolver, psi.get_nbands());
        }

        sum_iter += this->diag_once(phm_in, psi, eigenvalue_in_hsolver, is_occupied);

        ++ntry;

    } while (this->test_exit_cond(ntry, this->notconv, scf_type));

    if (notconv > std::max(5, psi.get_nbands() / 4))
    {
        std::cout << "\n notconv = " << this->notconv;
        std::cout << "\n Diago_DavSubspace::diag', too many bands are not converged! \n";
    }

    if (outputeigenvalue)
    {
        // output: eigenvalue_in_hsolver
        std::cout << "after dav, output eigenvalue_in_hsolver" << std::endl;
        for (size_t i = 0; i < psi.get_nbands(); i++)
        {
            std::cout << eigenvalue_in_hsolver[i] << "\t";
        }
        std::cout << std::endl;
    }

    if (outputscc)
    {
        std::cout << "after dav 222 " << std::endl;
        gemm_op<T, Device>()(this->ctx,
                             'C',
                             'N',
                             psi.get_nbands(),
                             psi.get_nbands(),
                             psi.get_current_nbas(),
                             this->one,
                             &psi(0, 0),
                             psi.get_current_nbas(),
                             &psi(0, 0),
                             psi.get_current_nbas(),
                             this->zero,
                             this->scc,
                             this->nbase_x);
        // output scc
        for (size_t i = 0; i < psi.get_nbands(); i++)
        {
            for (size_t j = 0; j < psi.get_nbands(); j++)
            {
                std::cout << this->scc[i * this->nbase_x + j] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
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
