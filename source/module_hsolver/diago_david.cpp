#include "diago_david.h"

#include "diago_iter_assist.h"
#include "module_base/blas_connector.h"
#include "module_base/constants.h"
#include "module_base/lapack_connector.h"
#include "module_base/memory.h"
#include "module_base/parallel_common.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"

#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif

using namespace hsolver;

template <typename T, typename Device> DiagoDavid<T, Device>::DiagoDavid(const Real* precondition_in)
{
    this->device = psi::device::get_device_type<Device>(this->ctx);
    this->precondition = precondition_in;

    test_david = 2;
    this->one = &this->cs.one;
    this->zero = &this->cs.zero;
    this->neg_one = &this->cs.neg_one;
    // 1: check which function is called and which step is executed
    // 2: check the eigenvalues of the result of each iteration
    // 3: check the eigenvalues and errors of the last result
    // default: no check
}

template <typename T, typename Device> DiagoDavid<T, Device>::~DiagoDavid()
{
    delmem_complex_op()(this->ctx, this->hphi);
    delmem_complex_op()(this->ctx, this->sphi);
    delmem_complex_op()(this->ctx, this->hcc);
    delmem_complex_op()(this->ctx, this->scc);
    delmem_complex_op()(this->ctx, this->vcc);
    delmem_complex_op()(this->ctx, this->lagrange_matrix);
    psi::memory::delete_memory_op<Real, psi::DEVICE_CPU>()(this->cpu_ctx, this->eigenvalue);
    if (this->device == psi::GpuDevice) {
        delmem_var_op()(this->ctx, this->d_precondition);
    }
}

template <typename T, typename Device>
void DiagoDavid<T, Device>::diag_mock(hamilt::Hamilt<T, Device>* phm_in,
                                           psi::Psi<T, Device>& psi,
                                           Real* eigenvalue_in)
{
    if (test_david == 1)
        ModuleBase::TITLE("DiagoDavid", "diag_mock");
    ModuleBase::timer::tick("DiagoDavid", "diag_mock");

    assert(DiagoDavid::PW_DIAG_NDIM > 1);
    assert(DiagoDavid::PW_DIAG_NDIM * psi.get_nbands() < psi.get_current_nbas() * GlobalV::NPROC_IN_POOL);
    // qianrui change it 2021-7-25.
    // In strictly speaking, it shoule be PW_DIAG_NDIM*nband < npw sum of all pools. We roughly estimate it here.
    // However, in most cases, total number of plane waves should be much larger than nband*PW_DIAG_NDIM

    /// initialize variables
    /// k_first = 0 means that nks is more like a dimension of "basis" to be contracted in "HPsi".In LR-TDDFT the formula writes :
    /// $$\sum_{ jb\mathbf{k}'}A^I_{ia\mathbf{k}, jb\mathbf{k}' }X ^ I_{ jb\mathbf{k}'}$$
    /// In the code :
    /// - "H" means A
    /// - "Psi" means X
    /// - "band" means the superscript I : the number of excited states to be solved
    /// - k : k-points, the same meaning as the ground state
    /// - "basis" : number of occupied ks-orbitals(subscripts i,j) * number of unoccupied ks-orbitals(subscripts a,b), corresponding to "bands" of the ground state
    this->dim = psi.get_k_first() ? psi.get_current_nbas() : psi.get_nk() * psi.get_nbasis();
    this->dmx = psi.get_k_first() ? psi.get_nbasis() : psi.get_nk() * psi.get_nbasis();
    this->n_band = psi.get_nbands();
    this->nbase_x = DiagoDavid::PW_DIAG_NDIM * this->n_band; // maximum dimension of the reduced basis set

    // the lowest N eigenvalues
    psi::memory::resize_memory_op<Real, psi::DEVICE_CPU>()(this->cpu_ctx, this->eigenvalue, this->nbase_x, "DAV::eig");
    psi::memory::set_memory_op<Real, psi::DEVICE_CPU>()(this->cpu_ctx, this->eigenvalue, 0, this->nbase_x);

    psi::Psi<T, Device> basis(1,
                                                 this->nbase_x,
                                                 this->dim,
                                                 &(psi.get_ngk(0))); // the reduced basis set
    ModuleBase::Memory::record("DAV::basis", this->nbase_x * this->dim * sizeof(T));

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // ModuleBase::ComplexMatrix hp(nbase_x, this->dim); // the product of H and psi in the reduced basis set
    resmem_complex_op()(this->ctx, this->hphi, this->nbase_x * this->dim, "DAV::hphi");
    setmem_complex_op()(this->ctx, this->hphi, 0, this->nbase_x * this->dim);

    // ModuleBase::ComplexMatrix sp(nbase_x, this->dim); // the Product of S and psi in the reduced basis set
    resmem_complex_op()(this->ctx, this->sphi, this->nbase_x * this->dim, "DAV::sphi");
    setmem_complex_op()(this->ctx, this->sphi, 0, this->nbase_x * this->dim);

    // ModuleBase::ComplexMatrix hc(this->nbase_x, this->nbase_x); // Hamiltonian on the reduced basis
    resmem_complex_op()(this->ctx, this->hcc, this->nbase_x * this->nbase_x, "DAV::hcc");
    setmem_complex_op()(this->ctx, this->hcc, 0, this->nbase_x * this->nbase_x);

    // ModuleBase::ComplexMatrix sc(this->nbase_x, this->nbase_x); // Overlap on the reduced basis
    resmem_complex_op()(this->ctx, this->scc, this->nbase_x * this->nbase_x, "DAV::scc");
    setmem_complex_op()(this->ctx, this->scc, 0, this->nbase_x * this->nbase_x);

    // ModuleBase::ComplexMatrix vc(this->nbase_x, this->nbase_x); // Eigenvectors of hc
    resmem_complex_op()(this->ctx, this->vcc, this->nbase_x * this->nbase_x, "DAV::vcc");
    setmem_complex_op()(this->ctx, this->vcc, 0, this->nbase_x * this->nbase_x);
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // convflag[m] = true if the m th band is convergent
    std::vector<bool> convflag(this->n_band, false);
    // unconv[m] store the number of the m th unconvergent band
    std::vector<int> unconv(this->n_band);

    int nbase = 0; // the dimension of the reduced basis set

    this->notconv = this->n_band; // the number of the unconvergent bands

    for (int m = 0; m < this->n_band; m++)
        unconv[m] = m;

    ModuleBase::timer::tick("DiagoDavid", "first");

    // orthogonalise the initial trial psi(0~nband-1)

    // ModuleBase::ComplexMatrix lagrange_matrix(this->n_band, this->n_band);
    resmem_complex_op()(this->ctx, this->lagrange_matrix, this->n_band * this->n_band);
    setmem_complex_op()(this->ctx, this->lagrange_matrix, 0, this->n_band * this->n_band);

    // plan for SchmitOrth
    std::vector<int> pre_matrix_mm_m(this->n_band, 0);
    std::vector<int> pre_matrix_mv_m(this->n_band, 1);
    this->planSchmitOrth(this->n_band, pre_matrix_mm_m.data(), pre_matrix_mv_m.data());

    for (int m = 0; m < this->n_band; m++)
    {
        if(GlobalV::use_paw)
        {
#ifdef USE_PAW
#ifdef __DEBUG
            assert(psi.get_k_first());
#endif 
            GlobalC::paw_cell.paw_nl_psi(1, reinterpret_cast<const std::complex<double>*> (&psi(m, 0)),
                reinterpret_cast<std::complex<double>*>(&this->sphi[m * this->dim]));
#endif
        }
        else
        {
            phm_in->sPsi(psi.get_k_first() ? &psi(m, 0) : &psi(m, 0, 0),
                         &this->sphi[m * this->dim],
                         this->dim,
                         this->dim,
                         1);
        }
    }
    // begin SchmitOrth
    for (int m = 0; m < this->n_band; m++)
    {
        // haozhihan replace 2022-10-23
        syncmem_complex_op()(this->ctx, this->ctx, &basis(m, 0), psi.get_k_first() ? &psi(m, 0) : &psi(m, 0, 0), this->dim);

        this->SchmitOrth(this->dim,
                         this->n_band,
                         m,
                         basis,
                         this->sphi,
                         &this->lagrange_matrix[m * this->n_band],
                         pre_matrix_mm_m[m],
                         pre_matrix_mv_m[m]);
        if(GlobalV::use_paw)
        {
#ifdef USE_PAW
            GlobalC::paw_cell.paw_nl_psi(1,reinterpret_cast<const std::complex<double>*> (&basis(m, 0)),
                reinterpret_cast<std::complex<double>*>(&this->sphi[m * this->dim]));
#endif
        }
        else
        {
            phm_in->sPsi(&basis(m, 0), &this->sphi[m * this->dim], this->dim, this->dim, 1);
        }
    }

    // end of SchmitOrth and calculate H|psi>
    hpsi_info dav_hpsi_in(&basis, psi::Range(1, 0, 0, this->n_band - 1), this->hphi);
    phm_in->ops->hPsi(dav_hpsi_in);

    this->cal_elem(this->dim, nbase, this->notconv, basis, this->hphi, this->sphi, this->hcc, this->scc);

    this->diag_zhegvx(nbase, this->n_band, this->hcc, this->scc, this->nbase_x, this->eigenvalue, this->vcc);

    for (int m = 0; m < this->n_band; m++)
    {
        eigenvalue_in[m] = this->eigenvalue[m];
    }

    ModuleBase::timer::tick("DiagoDavid", "first");

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
                       this->sphi,
                       this->vcc,
                       unconv.data(),
                       this->eigenvalue);

        this->cal_elem(this->dim, nbase, this->notconv, basis, this->hphi, this->sphi, this->hcc, this->scc);

        this->diag_zhegvx(nbase, this->n_band, this->hcc, this->scc, this->nbase_x, this->eigenvalue, this->vcc);

        // check convergence and update eigenvalues
        ModuleBase::timer::tick("DiagoDavid", "check_update");

        this->notconv = 0;
        for (int m = 0; m < this->n_band; m++)
        {
            convflag[m] = (std::abs(this->eigenvalue[m] - eigenvalue_in[m]) < DiagoIterAssist<T, Device>::PW_DIAG_THR);

            if (!convflag[m])
            {
                unconv[this->notconv] = m;
                this->notconv++;
            }

            eigenvalue_in[m] = this->eigenvalue[m];
        }

        ModuleBase::timer::tick("DiagoDavid", "check_update");
        if (!this->notconv || (nbase + this->notconv > this->nbase_x)
            || (dav_iter == DiagoIterAssist<T, Device>::PW_DIAG_NMAX))
        {
            ModuleBase::timer::tick("DiagoDavid", "last");

            // updata eigenvectors of Hamiltonian

            // ModuleBase::GlobalFunc::ZEROS(psi.get_pointer(), n_band * this->dmx);
            setmem_complex_op()(this->ctx, psi.get_pointer(), 0, n_band * this->dmx);
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            // haozhihan repalce 2022-10-18
            gemm_op<T, Device>()(this->ctx,
                                      'N',
                                      'N',
                                      this->dim,           // m: row of A,C
                                      this->n_band,        // n: col of B,C
                                      nbase,               // k: col of A, row of B
                                      this->one,
                                      basis.get_pointer(), // A dim * nbase
                                      this->dim,
                                      this->vcc,           // B nbase * n_band
                                      this->nbase_x,
                                      this->zero,
                                      psi.get_pointer(),   // C dim * n_band
                                      this->dmx
            );

            if (!this->notconv || (dav_iter == DiagoIterAssist<T, Device>::PW_DIAG_NMAX))
            {
                // overall convergence or last iteration: exit the iteration

                ModuleBase::timer::tick("DiagoDavid", "last");
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
                              eigenvalue_in,
                              psi,
                              basis,
                              this->hphi,
                              this->sphi,
                              this->hcc,
                              this->scc,
                              this->vcc);
                ModuleBase::timer::tick("DiagoDavid", "last");
            }

        } // end of if

    } while (1);

    DiagoIterAssist<T, Device>::avg_iter += static_cast<double>(dav_iter);

    ModuleBase::timer::tick("DiagoDavid", "diag_mock");

    return;
}

template <typename T, typename Device>
void DiagoDavid<T, Device>::cal_grad(hamilt::Hamilt<T, Device>* phm_in,
                                          const int& dim,
                                          const int& nbase, // current dimension of the reduced basis
                                          const int& notconv,
                                          psi::Psi<T, Device>& basis,
                                          T* hphi,
                                          T* sphi,
                                          const T* vcc,
                                          const int* unconv,
                                          const Real* eigenvalue)
{
    if (test_david == 1)
        ModuleBase::TITLE("DiagoDavid", "cal_grad");
    if (notconv == 0)
        return;
    ModuleBase::timer::tick("DiagoDavid", "cal_grad");

    // use template pointer for accelerate
    std::complex<double>* spsi;
    std::complex<double>* ppsi;

    // expand the reduced basis set with the new basis vectors P|Real(psi)>...
    // in which psi are the last eigenvectors
    // we define |Real(psi)> as (H-ES)*|Psi>, E = <psi|H|psi>/<psi|S|psi>

    // ModuleBase::ComplexMatrix vc_ev_vector(notconv, nbase);
    T* vc_ev_vector = nullptr;
    resmem_complex_op()(this->ctx, vc_ev_vector, notconv * nbase);
    setmem_complex_op()(this->ctx, vc_ev_vector, 0, notconv * nbase);

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // for (int m = 0; m < notconv; m++)
    // {
    //     for (int i = 0; i < nbase; i++)
    //     {
    //         // vc_ev_vector(m, i) = vc(i, unconv[m]);
    //         vc_ev_vector[m * nbase + i] = vcc[i * this->nbase_x + unconv[m]];
    //     }
    // }
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // replace by haozhihan
    for (int m = 0; m < notconv; m++)
    {
        syncmem_complex_op()(this->ctx,
                             this->ctx,
                             vc_ev_vector + m * nbase,
                             vcc + unconv[m] * this->nbase_x,
                             nbase);
    }


    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // haozhihan repalce 2022-10-18
    gemm_op<T, Device>()(this->ctx,
                              'N',
                              'N',
                              this->dim, // m: row of A,C
                              notconv, // n: col of B,C
                              nbase, // k: col of A, row of B
                              this->one, // alpha
                              hphi, // A dim * nbase
                              this->dim, // LDA: if(N) max(1,m) if(T) max(1,k)
                              vc_ev_vector, // B nbase * notconv
                              nbase, // LDB: if(N) max(1,k) if(T) max(1,n)
                              this->zero, // belta
                              &basis(nbase, 0), // C dim * notconv
                              this->dim // LDC: if(N) max(1, m)
    );
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // for (int m = 0; m < notconv; m++)
    // {
    //     for (int i = 0; i < nbase; i++)
    //     {
    //         vc_ev_vector[m * nbase + i] *= -1 * this->eigenvalue[unconv[m]];
    //     }
    // }
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // haozhihan replace 2022.11.18
    for (int m = 0; m < notconv; m++)
    {
        std::vector<Real> e_temp_cpu(nbase, (-1.0 * this->eigenvalue[unconv[m]]));

        if (this->device == psi::GpuDevice)
        {
#if defined(__CUDA) || defined(__ROCM)
            Real* e_temp_gpu = nullptr;
            resmem_var_op()(this->ctx, e_temp_gpu, nbase);
            syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, e_temp_gpu, e_temp_cpu.data(), nbase);
            vector_mul_vector_op<T, Device>()(this->ctx,
                                                   nbase,
                                                   vc_ev_vector + m * nbase,
                                                   vc_ev_vector + m * nbase,
                                                   e_temp_gpu);
            delmem_var_op()(this->ctx, e_temp_gpu);
#endif
        }
        else
        {
            vector_mul_vector_op<T, Device>()(this->ctx,
                                                   nbase,
                                                   vc_ev_vector + m * nbase,
                                                   vc_ev_vector + m * nbase,
                                                   e_temp_cpu.data());
        }
    }
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // haozhihan repalce 2022-10-18
    gemm_op<T, Device>()(this->ctx,
                              'N',
                              'N',
                              this->dim, // m: row of A,C
                              notconv, // n: col of B,C
                              nbase, // k: col of A, row of B
                              this->one, // alpha
                              sphi, // A
                              this->dim, // LDA: if(N) max(1,m) if(T) max(1,k)
                              vc_ev_vector, // B
                              nbase, // LDB: if(N) max(1,k) if(T) max(1,n)
                              this->one, // belta
                              &basis(nbase, 0), // C dim * notconv
                              this->dim // LDC: if(N) max(1, m)
    );
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    for (int m = 0; m < notconv; m++)
    {
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // haozhihan replace 2022-10-18
        if (this->device == psi::GpuDevice)
        {
#if defined(__CUDA) || defined(__ROCM)
            vector_div_vector_op<T, Device>()(this->ctx,
                                                   this->dim,
                                                   &basis(nbase + m, 0),
                                                   &basis(nbase + m, 0),
                                                   this->d_precondition);
#endif
        }
        else
        {
            vector_div_vector_op<T, Device>()(this->ctx,
                                                   this->dim,
                                                   &basis(nbase + m, 0),
                                                   &basis(nbase + m, 0),
                                                   this->precondition);
        }
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // for (int ig = 0; ig < this->dim; ig++)
        // {
        //     ppsi[ig] /= this->precondition[ig];
        // }
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    }

    // there is a nbase to nbase + notconv band orthogonalise
    // plan for SchmitOrth
    T* lagrange = nullptr;
    resmem_complex_op()(this->ctx, lagrange, notconv * (nbase + notconv));
    setmem_complex_op()(this->ctx, lagrange, 0, notconv * (nbase + notconv));

    std::vector<int> pre_matrix_mm_m(notconv, 0);
    std::vector<int> pre_matrix_mv_m(notconv, 1);
    this->planSchmitOrth(notconv, pre_matrix_mm_m.data(), pre_matrix_mv_m.data());
    for (int m = 0; m < notconv; m++)
    {
        if(GlobalV::use_paw)
        {
#ifdef USE_PAW
            GlobalC::paw_cell.paw_nl_psi(1,reinterpret_cast<const std::complex<double>*> (&basis(nbase + m, 0)),
                reinterpret_cast<std::complex<double>*>(&sphi[(nbase + m) * this->dim]));
#endif
        }
        else
        {
            phm_in->sPsi(&basis(nbase + m, 0), &sphi[(nbase + m) * this->dim], this->dim, this->dim, 1);
        }
    }
    // first nbase bands psi* dot notconv bands spsi to prepare lagrange_matrix

    // calculate the square matrix for future lagranges
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //  haozhihan replace 2022-10-18
    gemm_op<T, Device>()(this->ctx,
                              'C',
                              'N',
                              nbase, // m: row of A,C
                              notconv, // n: col of B,C
                              this->dim, // k: col of A, row of B
                              this->one, // alpha
                              &basis(0, 0), // A
                              this->dim, // LDA: if(N) max(1,m) if(T) max(1,k)
                              &sphi[nbase * this->dim], // B
                              this->dim, // LDB: if(N) max(1,k) if(T) max(1,n)
                              this->zero, // belta
                              lagrange, // C
                              nbase + notconv // LDC: if(N) max(1, m)
    );

    for (int m = 0; m < notconv; m++)
    {
        this->SchmitOrth(this->dim,
                         nbase + notconv,
                         nbase + m,
                         basis,
                         sphi,
                         &lagrange[m * (nbase + notconv)],
                         pre_matrix_mm_m[m],
                         pre_matrix_mv_m[m]);
        if(GlobalV::use_paw)
        {
#ifdef USE_PAW
            GlobalC::paw_cell.paw_nl_psi(1,reinterpret_cast<const std::complex<double>*> (&basis(nbase + m, 0)),
                reinterpret_cast<std::complex<double>*>(&sphi[(nbase + m) * this->dim]));
#endif
        }
        else
        {
            phm_in->sPsi(&basis(nbase + m, 0), &sphi[(nbase + m) * this->dim], this->dim, this->dim, 1);
        }
    }
    // calculate H|psi> for not convergence bands
    hpsi_info dav_hpsi_in(&basis,
                          psi::Range(1, 0, nbase, nbase + notconv - 1),
                          &hphi[nbase * this->dim]); // &hp(nbase, 0)
    phm_in->ops->hPsi(dav_hpsi_in);

    delmem_complex_op()(this->ctx, lagrange);
    delmem_complex_op()(this->ctx, vc_ev_vector);

    ModuleBase::timer::tick("DiagoDavid", "cal_grad");
    return;
}

template <typename T, typename Device>
void DiagoDavid<T, Device>::cal_elem(const int& dim,
                                          int& nbase, // current dimension of the reduced basis
                                          const int& notconv, // number of newly added basis vectors
                                          const psi::Psi<T, Device>& basis,
                                          const T* hphi,
                                          const T* sphi,
                                          T* hcc,
                                          T* scc)
{
    if (test_david == 1)
        ModuleBase::TITLE("DiagoDavid", "cal_elem");

    if (notconv == 0)
        return;
    ModuleBase::timer::tick("DiagoDavid", "cal_elem");

    gemm_op<T, Device>()(this->ctx,
                              'C',
                              'N',
                              notconv,
                              nbase + notconv,
                              this->dim,
                              this->one,
                              &basis(nbase, 0),   // this->dim * notconv
                              this->dim,
                              hphi,               // this->dim * (nbase + notconv)
                              this->dim,
                              this->zero,
                              hcc + nbase,        // notconv * (nbase + notconv)
                              this->nbase_x);

    gemm_op<T, Device>()(this->ctx,
                              'C',
                              'N',
                              notconv,
                              nbase + notconv,
                              this->dim,
                              this->one,
                              &basis(nbase, 0),   // this->dim * notconv
                              this->dim,
                              sphi,               // this->dim * (nbase + notconv)
                              this->dim,
                              this->zero,
                              scc + nbase,        // notconv * (nbase + notconv)
                              this->nbase_x);


#ifdef __MPI
    if (GlobalV::NPROC_IN_POOL > 1)
    {
        matrixTranspose_op<T, Device>()(this->ctx, this->nbase_x, this->nbase_x, hcc, hcc);
        matrixTranspose_op<T, Device>()(this->ctx, this->nbase_x, this->nbase_x, scc, scc);

        auto* swap = new T[notconv * this->nbase_x];
        syncmem_complex_op()(this->ctx, this->ctx, swap, hcc + nbase * this->nbase_x, notconv * this->nbase_x);
        if (std::is_same<T, double>::value)
        {
            Parallel_Reduce::reduce_pool(hcc + nbase * this->nbase_x, notconv * this->nbase_x);
        }
        else
        {
            if (psi::device::get_current_precision(swap) == "single") {
                MPI_Reduce(swap, hcc + nbase * this->nbase_x, notconv * this->nbase_x, MPI_COMPLEX, MPI_SUM, 0, POOL_WORLD);
            }
            else {
                MPI_Reduce(swap, hcc + nbase * this->nbase_x, notconv * this->nbase_x, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, POOL_WORLD);
            }
            syncmem_complex_op()(this->ctx, this->ctx, swap, scc + nbase * this->nbase_x, notconv * this->nbase_x);
            if (psi::device::get_current_precision(swap) == "single") {
                MPI_Reduce(swap, scc + nbase * this->nbase_x, notconv * this->nbase_x, MPI_COMPLEX, MPI_SUM, 0, POOL_WORLD);
            }
            else {
                MPI_Reduce(swap, scc + nbase * this->nbase_x, notconv * this->nbase_x, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, POOL_WORLD);
            }
        }
        delete[] swap;

        // Parallel_Reduce::reduce_complex_double_pool( hcc + nbase * this->nbase_x, notconv * this->nbase_x );
        // Parallel_Reduce::reduce_complex_double_pool( scc + nbase * this->nbase_x, notconv * this->nbase_x );

        matrixTranspose_op<T, Device>()(this->ctx, this->nbase_x, this->nbase_x, hcc, hcc);
        matrixTranspose_op<T, Device>()(this->ctx, this->nbase_x, this->nbase_x, scc, scc);
    }
#endif

    nbase += notconv;
    ModuleBase::timer::tick("DiagoDavid", "cal_elem");
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
                                             const T* scc,
                                             const int& nbase_x,
                                             Real* eigenvalue, // in CPU
                                             T* vcc)
{
    //	ModuleBase::TITLE("DiagoDavid","diag_zhegvx");
    ModuleBase::timer::tick("DiagoDavid", "diag_zhegvx");
    if (GlobalV::RANK_IN_POOL == 0)
    {
        assert(nbase_x >= std::max(1, nbase));

        if (this->device == psi::GpuDevice)
        {
#if defined(__CUDA) || defined(__ROCM)
            Real* eigenvalue_gpu = nullptr;
            resmem_var_op()(this->ctx, eigenvalue_gpu, this->nbase_x);
            syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, eigenvalue_gpu, this->eigenvalue, this->nbase_x);

            dnevx_op<T, Device>()(this->ctx, nbase, this->nbase_x, this->hcc, nband, eigenvalue_gpu, this->vcc);

            syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, this->eigenvalue, eigenvalue_gpu, this->nbase_x);
            delmem_var_op()(this->ctx, eigenvalue_gpu);
#endif
        }
        else
        {
            dnevx_op<T, Device>()(this->ctx, nbase, this->nbase_x, this->hcc, nband, this->eigenvalue, this->vcc);
        }
    }

#ifdef __MPI
    if (GlobalV::NPROC_IN_POOL > 1)
    {
        // vcc: nbase * nband
        for (int i = 0; i < nband; i++)
        {
            MPI_Bcast(&vcc[i * this->nbase_x], nbase, MPI_DOUBLE_COMPLEX, 0, POOL_WORLD);
        }
        MPI_Bcast(this->eigenvalue, nband, MPI_DOUBLE, 0, POOL_WORLD);
    }
#endif

    ModuleBase::timer::tick("DiagoDavid", "diag_zhegvx");
    return;
}

template <typename T, typename Device>
void DiagoDavid<T, Device>::refresh(const int& dim,
                                         const int& nband,
                                         int& nbase,
                                         const Real* eigenvalue_in,
                                         const psi::Psi<T, Device>& psi,
                                         psi::Psi<T, Device>& basis,
                                         T* hp,
                                         T* sp,
                                         T* hc,
                                         T* sc,
                                         T* vc)
{
    if (test_david == 1)
        ModuleBase::TITLE("DiagoDavid", "refresh");
    ModuleBase::timer::tick("DiagoDavid", "refresh");

    // update hp,sp
    basis.zero_out();
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // haozhihan repalce 2022-10-18
    gemm_op<T, Device>()(this->ctx,
                              'N',
                              'N',
                              this->dim,            // m: row of A,C
                              nband,                // n: col of B,C
                              nbase,                // k: col of A, row of B
                              this->one,
                              this->hphi,           // A dim * nbase
                              this->dim,
                              this->vcc,            // B nbase * nband
                              this->nbase_x,
                              this->zero,
                              basis.get_pointer(),  // C dim * nband
                              this->dim
    );

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // haozhihan repalce 2022-10-18
    gemm_op<T, Device>()(this->ctx,
                              'N',
                              'N',
                              this->dim,                // m: row of A,C
                              nband,                    // n: col of B,C
                              nbase,                    // k: col of A, row of B
                              this->one,
                              this->sphi,               // A dim * nbase
                              this->dim,
                              this->vcc,                // B nbase * nband
                              this->nbase_x,
                              this->zero,
                              &basis(nband, 0),         // C dim * nband
                              this->dim
    );

    syncmem_complex_op()(this->ctx, this->ctx, hphi, &basis(0, 0), this->dim * nband);
    syncmem_complex_op()(this->ctx, this->ctx, sphi, &basis(nband, 0), this->dim * nband);
    /*for (int m = 0; m < nband; m++)
    {
        for (int ig = 0; ig < this->dim; ig++)
        {
            hp(m, ig) = basis(m, ig);
            sp(m, ig) = basis(m + nband, ig);
        }
    }*/

    // update basis
    basis.zero_out();
    for (int m = 0; m < nband; m++)
    {
        syncmem_complex_op()(this->ctx, this->ctx, &basis(m, 0), psi.get_k_first() ? &psi(m, 0) : &psi(m, 0, 0), this->dim);
        /*for (int ig = 0; ig < npw; ig++)
            basis(m, ig) = psi(m, ig);*/
    }

    // updata the reduced Hamiltonian
    nbase = nband;

    // hc.zero_out();
    setmem_complex_op()(this->ctx, hcc, 0, this->nbase_x * this->nbase_x);

    // sc.zero_out();
    setmem_complex_op()(this->ctx, scc, 0, this->nbase_x * this->nbase_x);

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    if (this->device == psi::GpuDevice)
    {
#if defined(__CUDA) || defined(__ROCM)
        T* hcc_cpu = nullptr;
        T* scc_cpu = nullptr;
        T* vcc_cpu = nullptr;
        psi::memory::resize_memory_op<T, psi::DEVICE_CPU>()(this->cpu_ctx,
                                                                               hcc_cpu,
                                                                               this->nbase_x * this->nbase_x,
                                                                               "DAV::hcc");
        psi::memory::resize_memory_op<T, psi::DEVICE_CPU>()(this->cpu_ctx,
                                                                               scc_cpu,
                                                                               this->nbase_x * this->nbase_x,
                                                                               "DAV::scc");
        psi::memory::resize_memory_op<T, psi::DEVICE_CPU>()(this->cpu_ctx,
                                                                               vcc_cpu,
                                                                               this->nbase_x * this->nbase_x,
                                                                               "DAV::vcc");

        syncmem_d2h_op()(this->cpu_ctx, this->ctx, hcc_cpu, hcc, this->nbase_x * this->nbase_x);
        syncmem_d2h_op()(this->cpu_ctx, this->ctx, scc_cpu, scc, this->nbase_x * this->nbase_x);
        syncmem_d2h_op()(this->cpu_ctx, this->ctx, vcc_cpu, vcc, this->nbase_x * this->nbase_x);

        for (int i = 0; i < nbase; i++)
        {
            hcc_cpu[i * this->nbase_x + i] = eigenvalue_in[i];
            scc_cpu[i * this->nbase_x + i] = this->one[0];
            vcc_cpu[i * this->nbase_x + i] = this->one[0];
        }

        syncmem_h2d_op()(this->ctx, this->cpu_ctx, hcc, hcc_cpu, this->nbase_x * this->nbase_x);
        syncmem_h2d_op()(this->ctx, this->cpu_ctx, scc, scc_cpu, this->nbase_x * this->nbase_x);
        syncmem_h2d_op()(this->ctx, this->cpu_ctx, vcc, vcc_cpu, this->nbase_x * this->nbase_x);

        psi::memory::delete_memory_op<T, psi::DEVICE_CPU>()(this->cpu_ctx, hcc_cpu);
        psi::memory::delete_memory_op<T, psi::DEVICE_CPU>()(this->cpu_ctx, scc_cpu);
        psi::memory::delete_memory_op<T, psi::DEVICE_CPU>()(this->cpu_ctx, vcc_cpu);
#endif
    }
    else
    {
        for (int i = 0; i < nbase; i++)
        {
            hcc[i * this->nbase_x + i] = eigenvalue_in[i];
            // sc(i, i) = this->one;
            scc[i * this->nbase_x + i] = this->one[0];
            // vc(i, i) = this->one;
            vcc[i * this->nbase_x + i] = this->one[0];
        }
    }
    ModuleBase::timer::tick("DiagoDavid", "refresh");
    return;
}

template <typename T, typename Device>
void DiagoDavid<T, Device>::SchmitOrth(const int& dim,
                                            const int nband,
                                            const int m,
                                            psi::Psi<T, Device>& basis,
                                            const T* sphi,
                                            T* lagrange_m,
                                            const int mm_size,
                                            const int mv_size)
{
    //	if(test_david == 1) ModuleBase::TITLE("DiagoDavid","SchmitOrth");
    ModuleBase::timer::tick("DiagoDavid", "SchmitOrth");

    // orthogonalize starting eigenfunction to those already calculated
    // psi_m orthogonalize to psi(0) ~ psi(m-1)
    // Attention, the orthogonalize here read as
    // psi(m) -> psi(m) - \sum_{i < m} \langle psi(i)|S|psi(m) \rangle psi(i)
    // so the orthogonalize is performed about S.

    assert(basis.get_nbands() >= nband);
    assert(m >= 0);
    assert(m < nband);

    T* psi_m = &basis(m, 0);

    // std::complex<double> *lagrange = new std::complex<double>[m + 1];
    // ModuleBase::GlobalFunc::ZEROS(lagrange, m + 1);

    // calculate the square matrix for future lagranges
    if (mm_size != 0)
    {
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // haozhihan repalce 2022-10-16
        gemm_op<T, Device>()(this->ctx,
                                  'C',
                                  'N',
                                  mm_size, // m: row of A,C
                                  mm_size, // n: col of B,C
                                  this->dim, // k: col of A, row of B
                                  this->one, // alpha
                                  &basis(m - mv_size + 1 - mm_size, 0), // A
                                  this->dim, // LDA: if(N) max(1,m) if(T) max(1,k)
                                  &sphi[m * this->dim], // B
                                  this->dim, // LDB: if(N) max(1,k) if(T) max(1,n)
                                  this->zero, // belta
                                  &lagrange_m[m - mv_size + 1 - mm_size], // C
                                  nband // LDC: if(N) max(1, m)
        );
    }
    // calculate other lagranges for this band
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //  haozhihan repalce 2022-10-16
    gemv_op<T, Device>()(this->ctx,
                              'C',
                              this->dim,
                              mv_size,
                              this->one,
                              &basis(m - mv_size + 1, 0),
                              this->dim,
                              &sphi[m * this->dim],
                              1,
                              this->zero,
                              &lagrange_m[m - mv_size + 1],
                              1);
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    Parallel_Reduce::reduce_pool(lagrange_m, m + 1);

    T var = *this->zero;
    syncmem_d2h_op()(this->cpu_ctx, this->ctx, &var, lagrange_m + m, 1);
    double psi_norm = get_real(var);

    assert(psi_norm > 0.0);

    // haozhihan replace 2022-10-24
    gemv_op<T, Device>()(this->ctx,
                              'N',
                              this->dim,
                              m,
                              this->neg_one,
                              &basis(0, 0),
                              this->dim,
                              lagrange_m,
                              1,
                              this->one,
                              psi_m,
                              1);

    psi_norm -= dot_real_op<T, Device>()(this->ctx, m, lagrange_m, lagrange_m, false);

    // for (int j = 0; j < m; j++)
    // {
    //     const std::complex<double> alpha = std::complex<double>(-1, 0) * lagrange_m[j];
    //     zaxpy_(&npw, &alpha, &psi(j,0), &inc, psi_m, &inc);
    //     /*for (int ig = 0; ig < npw; ig++)
    //     {
    //         psi_m[ig] -= lagrange[j] * psi(j, ig);
    //     }*/
    //     psi_norm -= (conj(lagrange_m[j]) * lagrange_m[j]).real();
    // }

    assert(psi_norm > 0.0);

    psi_norm = sqrt(psi_norm);

    if (psi_norm < 1.0e-12)
    {
        std::cout << "DiagoDavid::SchmitOrth:aborted for psi_norm <1.0e-12" << std::endl;
        std::cout << "nband = " << nband << std::endl;
        std::cout << "m = " << m << std::endl;
        exit(0);
    }
    else
    {
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // haozhihan repalce 2022-10-16
        vector_div_constant_op<T, Device>()(this->ctx, this->dim, psi_m, psi_m, psi_norm);
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // for (int i = 0; i < npw; i++)
        // {
        //     psi_m[i] /= psi_norm;
        // }
    }

    // delete[] lagrange;
    ModuleBase::timer::tick("DiagoDavid", "SchmitOrth");
    return;
}

template <typename T, typename Device>
void DiagoDavid<T, Device>::planSchmitOrth(const int nband, int* pre_matrix_mm_m, int* pre_matrix_mv_m)
{
    if (nband <= 0)
        return;
    ModuleBase::GlobalFunc::ZEROS(pre_matrix_mm_m, nband);
    ModuleBase::GlobalFunc::ZEROS(pre_matrix_mv_m, nband);
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
            if (res_nband == matrix_size)
                pre_matrix_mv_m[index] = 1;
            else
                pre_matrix_mv_m[index] = 2;
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
void DiagoDavid<T, Device>::diag(hamilt::Hamilt<T, Device>* phm_in,
                                      psi::Psi<T, Device>& psi,
                                      Real* eigenvalue_in)
{
    /// record the times of trying iterative diagonalization
    int ntry = 0;
    this->notconv = 0;

#if defined(__CUDA) || defined(__ROCM)
    if (this->device == psi::GpuDevice)
    {
        resmem_var_op()(this->ctx, this->d_precondition, psi.get_nbasis());
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, this->d_precondition, this->precondition, psi.get_nbasis());
    }
#endif

    do
    {
        this->diag_mock(phm_in, psi, eigenvalue_in);
        ++ntry;
    } while (DiagoIterAssist<T, Device>::test_exit_cond(ntry, this->notconv));

    if (notconv > std::max(5, psi.get_nbands() / 4))
    {
        std::cout << "\n notconv = " << this->notconv;
        std::cout << "\n DiagoDavid::diag', too many bands are not converged! \n";
    }
    return;
}

namespace hsolver {
template class DiagoDavid<std::complex<float>, psi::DEVICE_CPU>;
template class DiagoDavid<std::complex<double>, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoDavid<std::complex<float>, psi::DEVICE_GPU>;
template class DiagoDavid<std::complex<double>, psi::DEVICE_GPU>;
#endif
#ifdef __LCAO
template class DiagoDavid<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoDavid<double, psi::DEVICE_GPU>;
#endif
#endif
} // namespace hsolver
