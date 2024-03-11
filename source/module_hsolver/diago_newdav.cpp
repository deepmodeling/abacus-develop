#include "diago_newdav.h"

#include <algorithm>
#include <type_traits>

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

using namespace hsolver;

inline double get_real(const double& x)
{
    return x;
}
inline double get_real(const std::complex<double>& x)
{
    return x.real();
}
inline float get_real(const std::complex<float>& x)
{
    return x.real();
}

inline void set_value(double& input, double x, double y)
{
    input = x;
}
inline void set_value(std::complex<double>& input, double x, double y)
{
    input = {x, y};
}
inline void set_value(std::complex<float>& input, double x, double y)
{
    input = {static_cast<float>(x), static_cast<float>(y)};
}

inline std::complex<double> set_real_tocomplex(const std::complex<double>& x)
{
    return {x.real(), 0.0};
}

inline std::complex<float> set_real_tocomplex(const std::complex<float>& x)
{
    return {x.real(), 0.0};
}

inline double set_real_tocomplex(const double& x)
{
    return x;
}

inline std::complex<double> get_conj(const std::complex<double>& x)
{
    return {x.real(), -x.imag()};
}

inline std::complex<float> get_conj(const std::complex<float>& x)
{
    return {x.real(), -x.imag()};
}

inline double get_conj(const double& x)
{
    return x;
}

template <typename T, typename Device>
Diago_NewDav<T, Device>::Diago_NewDav(const Real* precondition_in)
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

template <typename T, typename Device>
Diago_NewDav<T, Device>::~Diago_NewDav()
{
    delmem_complex_op()(this->ctx, this->hphi);
    delmem_complex_op()(this->ctx, this->sphi);
    delmem_complex_op()(this->ctx, this->hcc);
    delmem_complex_op()(this->ctx, this->scc);
    delmem_complex_op()(this->ctx, this->vcc);
    delmem_complex_op()(this->ctx, this->lagrange_matrix);
    psi::memory::delete_memory_op<Real, psi::DEVICE_CPU>()(this->cpu_ctx, this->eigenvalue);
    if (this->device == psi::GpuDevice)
    {
        delmem_var_op()(this->ctx, this->d_precondition);
    }
}

template <typename T, typename Device>
void Diago_NewDav<T, Device>::diag_mock(hamilt::Hamilt<T, Device>* phm_in,
                                        psi::Psi<T, Device>& psi,
                                        Real* eigenvalue_in)
{
    if (test_david == 1)
        ModuleBase::TITLE("Diago_NewDav", "diag_mock");
    ModuleBase::timer::tick("Diago_NewDav", "diag_mock");

    assert(Diago_NewDav::PW_DIAG_NDIM > 1);
    assert(Diago_NewDav::PW_DIAG_NDIM * psi.get_nbands() < psi.get_current_nbas() * GlobalV::NPROC_IN_POOL);
    // qianrui change it 2021-7-25.
    // In strictly speaking, it shoule be PW_DIAG_NDIM*nband < npw sum of all pools. We roughly estimate it here.
    // However, in most cases, total number of plane waves should be much larger than nband*PW_DIAG_NDIM

    /// initialize variables
    /// k_first = 0 means that nks is more like a dimension of "basis" to be contracted in "HPsi".In LR-TDDFT the
    /// formula writes :
    /// $$\sum_{ jb\mathbf{k}'}A^I_{ia\mathbf{k}, jb\mathbf{k}' }X ^ I_{ jb\mathbf{k}'}$$
    /// In the code :
    /// - "H" means A
    /// - "Psi" means X
    /// - "band" means the superscript I : the number of excited states to be solved
    /// - k : k-points, the same meaning as the ground state
    /// - "basis" : number of occupied ks-orbitals(subscripts i,j) * number of unoccupied ks-orbitals(subscripts a,b),
    /// corresponding to "bands" of the ground state
    this->dim = psi.get_k_first() ? psi.get_current_nbas() : psi.get_nk() * psi.get_nbasis();
    this->dmx = psi.get_k_first() ? psi.get_nbasis() : psi.get_nk() * psi.get_nbasis();
    this->n_band = psi.get_nbands();
    // this->nbase_x = Diago_NewDav::PW_DIAG_NDIM * this->n_band; // maximum dimension of the reduced basis set
    this->nbase_x = 2 * this->n_band;

    // the lowest N eigenvalues
    psi::memory::resize_memory_op<Real, psi::DEVICE_CPU>()(this->cpu_ctx, this->eigenvalue, this->nbase_x, "DAV::eig");
    psi::memory::set_memory_op<Real, psi::DEVICE_CPU>()(this->cpu_ctx, this->eigenvalue, 0, this->nbase_x);

    psi::Psi<T, Device> basis(1, this->nbase_x, this->dim,
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

    std::cout << std::fixed << std::setprecision(7);

    ModuleBase::timer::tick("Diago_NewDav", "first");

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

    this->cal_elem(this->dim, nbase, this->notconv, basis, this->hphi, this->sphi, this->hcc, this->scc, true);

    this->diag_zhegvx(nbase, this->n_band, this->hcc, this->scc, this->nbase_x, this->eigenvalue, this->vcc, true);

    ModuleBase::timer::tick("Diago_NewDav", "first");

    int dav_iter = 0;
    do
    {
        dav_iter++;

        // std::cout << "dav iter: " << dav_iter << " nbase == " << nbase << std::endl;

        // for (int m = 0; m < nbase; m++)
        // {
        //     phm_in->sPsi(basis.get_k_first() ? &basis(m, 0) : &basis(m, 0, 0),
        //                  &this->sphi[m * this->dim],
        //                  this->dim,
        //                  this->dim,
        //                  1);
        // }

        this->cal_grad(phm_in,
                       this->dim,
                       nbase,
                       this->notconv,
                       basis,
                       this->hphi,
                       this->vcc,
                       unconv.data(),
                       this->eigenvalue);

        this->cal_elem(this->dim, nbase, this->notconv, basis, this->hphi, this->sphi, this->hcc, this->scc, false);

        this->diag_zhegvx(nbase, this->n_band, this->hcc, this->scc, this->nbase_x, this->eigenvalue, this->vcc, false);

        // check convergence and update eigenvalues
        ModuleBase::timer::tick("Diago_NewDav", "check_update");

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

        ModuleBase::timer::tick("Diago_NewDav", "check_update");
        if (this->notconv == 0 || (nbase + this->notconv + 1 > this->nbase_x)
            || (dav_iter == DiagoIterAssist<T, Device>::PW_DIAG_NMAX))
        {
            ModuleBase::timer::tick("Diago_NewDav", "last");

            // updata eigenvectors of Hamiltonian

            // ModuleBase::GlobalFunc::ZEROS(psi.get_pointer(), n_band * this->dmx);
            setmem_complex_op()(this->ctx, psi.get_pointer(), 0, n_band * this->dmx);
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
                                 this->dmx);

            if (!this->notconv || (dav_iter == DiagoIterAssist<T, Device>::PW_DIAG_NMAX))
            {
                // overall convergence or last iteration: exit the iteration

                ModuleBase::timer::tick("Diago_NewDav", "last");
                break;
            }
            else
            {
                // if the dimension of the reduced basis set is becoming too large,
                // then replace the first N (=nband) basis vectors with the current
                // estimate of the eigenvectors and set the basis dimension to N;

                // std::cout << dav_iter << " refresh" << std::endl;

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
                ModuleBase::timer::tick("Diago_NewDav", "last");

                // std::cout << "after refresh output scc: " << std::endl;
                // for (size_t i = 0; i < this->nbase_x; i++)
                // {
                //     for (size_t j = 0; j < this->nbase_x; j++)
                //     {
                //         std::cout << scc[i * this->nbase_x + j] << "\t";
                //     }
                //     std::cout << std::endl;
                // }
                // std::cout << std::endl;
            }

        } // end of if

    } while (1);

    // std::cout << "dav_iter == " << dav_iter << std::endl;

    DiagoIterAssist<T, Device>::avg_iter += static_cast<double>(dav_iter);

    ModuleBase::timer::tick("Diago_NewDav", "diag_mock");

    return;
}

// /*

template <typename T, typename Device>
void Diago_NewDav<T, Device>::cal_grad(hamilt::Hamilt<T, Device>* phm_in,
                                       const int& dim,
                                       const int& nbase, // current dimension of the reduced basis
                                       const int& notconv,
                                       psi::Psi<T, Device>& basis,
                                       T* hphi,
                                       T* vcc,
                                       const int* unconv,
                                       Real* eigenvalue)
{
    for (size_t i = 0; i < notconv; i++)
    {
        if (unconv[i] != i)
        {
            // std::cout << "unconv[i]: " << unconv[i] << ", i = " << i << std::endl;
            syncmem_complex_op()(this->ctx, this->ctx, vcc + i * this->nbase_x, vcc + unconv[i] * this->nbase_x, nbase);
            this->eigenvalue[i] = this->eigenvalue[unconv[i]];
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

        std::vector<Real> e_temp_cpu(this->dim, (-this->eigenvalue[m]));

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

    // precondition!!!
    // std::vector<Real> pre(this->dim, 0.0);
    for (int m = 0; m < notconv; m++)
    {
        // for (size_t i = 0; i < this->dim; i++)
        // {
        //     std::cout << "this->precondition[i] i==" << i << "  == " << this->precondition[i] << std::endl;
        //     double x = this->precondition[i] - this->eigenvalue[m];
        //     pre[i] = 0.5 * (1.0 + x + sqrt(1 + (x - 1.0) * (x - 1.0)));
        // }

        vector_div_vector_op<T, Device>()(this->ctx,
                                          this->dim,
                                          &basis(nbase + m, 0),
                                          &basis(nbase + m, 0),
                                          this->precondition);
    }

    // /*
    //  * "normalize" correction vectors psi(:, nbase+1:nbase+notcnv)
    //  * in order to improve numerical stability of subspace diagonalization
    //  *
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

    // calculate H|psi> for not convergence bands
    hpsi_info dav_hpsi_in(&basis, psi::Range(1, 0, nbase, nbase + notconv - 1), &hphi[nbase * this->dim]);
    phm_in->ops->hPsi(dav_hpsi_in);

    return;
}

// */

/*
template <typename T, typename Device>
void Diago_NewDav<T, Device>::cal_grad(hamilt::Hamilt<T, Device>* phm_in,
                                       const int& dim,
                                       const int& nbase, // current dimension of the reduced basis
                                       const int& notconv,
                                       psi::Psi<T, Device>& basis,
                                       T* hphi,
                                       T* sphi,
                                       T* vcc,
                                       const int* unconv,
                                       Real* eigenvalue)
{
    if (test_david == 1)
        ModuleBase::TITLE("Diago_NewDav", "cal_grad");
    if (notconv == 0)
        return;
    ModuleBase::timer::tick("Diago_NewDav", "cal_grad");

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
        syncmem_complex_op()(this->ctx, this->ctx, vc_ev_vector + m * nbase, vcc + unconv[m] * this->nbase_x, nbase);
    }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // haozhihan repalce 2022-10-18
    gemm_op<T, Device>()(this->ctx,
                         'N',
                         'N',
                         this->dim,        // m: row of A,C
                         notconv,          // n: col of B,C
                         nbase,            // k: col of A, row of B
                         this->one,        // alpha
                         hphi,             // A dim * nbase
                         this->dim,        // LDA: if(N) max(1,m) if(T) max(1,k)
                         vc_ev_vector,     // B nbase * notconv
                         nbase,            // LDB: if(N) max(1,k) if(T) max(1,n)
                         this->zero,       // belta
                         &basis(nbase, 0), // C dim * notconv
                         this->dim         // LDC: if(N) max(1, m)
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
                         this->dim,        // m: row of A,C
                         notconv,          // n: col of B,C
                         nbase,            // k: col of A, row of B
                         this->one,        // alpha
                         sphi,             // A
                         this->dim,        // LDA: if(N) max(1,m) if(T) max(1,k)
                         vc_ev_vector,     // B
                         nbase,            // LDB: if(N) max(1,k) if(T) max(1,n)
                         this->one,        // belta
                         &basis(nbase, 0), // C dim * notconv
                         this->dim         // LDC: if(N) max(1, m)
    );
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // from there, basis(0,0) - basis(nbase, 0) is init guess
    // basis(nbase, 0) - basis(nbase_x, dim) is new data
    // do precondition operator
    for (int m = 0; m < notconv; m++)
    {
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

    // /*
    //  * "normalize" correction vectors psi(:, nbase+1:nbase+notcnv)
    //  * in order to improve numerical stability of subspace diagonalization
    std::vector<Real> psi_norm(notconv, 0.0);
    for (size_t i = 0; i < notconv; i++)
    {
        psi_norm[i] = dot_real_op<T, Device>()(this->ctx, this->dim, &basis(nbase + i, 0), &basis(nbase + i, 0), false);
        assert(psi_norm[i] > 0.0);
        psi_norm[i] = sqrt(psi_norm[i]);

        vector_div_constant_op<T, Device>()(this->ctx, this->dim, &basis(nbase + i, 0), &basis(nbase + i, 0),
psi_norm[i]);
    }

    // calculate H|psi> for not convergence bands
    hpsi_info dav_hpsi_in(&basis,
                          psi::Range(1, 0, nbase, nbase + notconv - 1),
                          &hphi[nbase * this->dim]); // &hp(nbase, 0)
    phm_in->ops->hPsi(dav_hpsi_in);

    // delmem_complex_op()(this->ctx, lagrange);
    delmem_complex_op()(this->ctx, vc_ev_vector);

    ModuleBase::timer::tick("Diago_NewDav", "cal_grad");
    return;
}

*/

template <typename T, typename Device>
void Diago_NewDav<T, Device>::cal_elem(const int& dim,
                                       int& nbase,         // current dimension of the reduced basis
                                       const int& notconv, // number of newly added basis vectors
                                       const psi::Psi<T, Device>& basis,
                                       const T* hphi,
                                       const T* sphi,
                                       T* hcc,
                                       T* scc,
                                       bool init)
{
    if (test_david == 1)
        ModuleBase::TITLE("Diago_NewDav", "cal_elem");

    if (notconv == 0)
        return;
    ModuleBase::timer::tick("Diago_NewDav", "cal_elem");

    if (init)
    {
        gemm_op<T, Device>()(this->ctx,
                            'C',
                            'N',
                            notconv,
                            nbase + notconv,
                            this->dim,
                            this->one,
                            &basis(nbase, 0), // this->dim * notconv
                            this->dim,
                            hphi, // this->dim * (nbase + notconv)
                            this->dim,
                            this->zero,
                            hcc + nbase, // notconv * (nbase + notconv)
                            this->nbase_x);

        gemm_op<T, Device>()(this->ctx,
                            'C',
                            'N',
                            notconv,
                            nbase + notconv,
                            this->dim,
                            this->one,
                            &basis(nbase, 0), // this->dim * notconv
                            this->dim,
                            &basis(nbase, 0), // this->dim * (nbase + notconv)
                            this->dim,
                            this->zero,
                            scc + nbase, // notconv * (nbase + notconv)
                            this->nbase_x);
    } 
    else
    {
        gemm_op<T, Device>()(this->ctx,
                            'C',
                            'N',
                            notconv,
                            nbase + notconv,
                            this->dim,
                            this->one,
                            &hphi[nbase * this->dim], // this->dim * notconv
                            this->dim,
                            &basis(0, 0), // this->dim * (nbase + notconv)
                            this->dim,
                            this->zero,
                            hcc + nbase, // notconv * (nbase + notconv)
                            this->nbase_x);

        gemm_op<T, Device>()(this->ctx,
                            'C',
                            'N',
                            notconv,
                            nbase + notconv,
                            this->dim,
                            this->one,
                            &basis(nbase, 0), // this->dim * notconv
                            this->dim,
                            &basis(0, 0), // this->dim * (nbase + notconv)
                            this->dim,
                            this->zero,
                            scc + nbase, // notconv * (nbase + notconv)
                            this->nbase_x);
    }

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
            if (psi::device::get_current_precision(swap) == "single")
            {
                MPI_Reduce(swap,
                           hcc + nbase * this->nbase_x,
                           notconv * this->nbase_x,
                           MPI_COMPLEX,
                           MPI_SUM,
                           0,
                           POOL_WORLD);
            }
            else
            {
                MPI_Reduce(swap,
                           hcc + nbase * this->nbase_x,
                           notconv * this->nbase_x,
                           MPI_DOUBLE_COMPLEX,
                           MPI_SUM,
                           0,
                           POOL_WORLD);
            }
            syncmem_complex_op()(this->ctx, this->ctx, swap, scc + nbase * this->nbase_x, notconv * this->nbase_x);
            if (psi::device::get_current_precision(swap) == "single")
            {
                MPI_Reduce(swap,
                           scc + nbase * this->nbase_x,
                           notconv * this->nbase_x,
                           MPI_COMPLEX,
                           MPI_SUM,
                           0,
                           POOL_WORLD);
            }
            else
            {
                MPI_Reduce(swap,
                           scc + nbase * this->nbase_x,
                           notconv * this->nbase_x,
                           MPI_DOUBLE_COMPLEX,
                           MPI_SUM,
                           0,
                           POOL_WORLD);
            }
        }
        delete[] swap;

        // Parallel_Reduce::reduce_complex_double_pool( hcc + nbase * this->nbase_x, notconv * this->nbase_x );
        // Parallel_Reduce::reduce_complex_double_pool( scc + nbase * this->nbase_x, notconv * this->nbase_x );

        matrixTranspose_op<T, Device>()(this->ctx, this->nbase_x, this->nbase_x, hcc, hcc);
        matrixTranspose_op<T, Device>()(this->ctx, this->nbase_x, this->nbase_x, scc, scc);
    }
#endif

    int nb1 = nbase;

    nbase += notconv;

    // reset:
    if (init)
    {
        for (size_t i = 0; i < nbase; i++)
        {
            if (i >= nb1)
            {
                hcc[i * this->nbase_x + i] = set_real_tocomplex(hcc[i * this->nbase_x + i]);
                scc[i * this->nbase_x + i] = cs.one;
            }
            for (size_t j = std::max(i + 1, (size_t)nb1); j < nbase; j++)
            {
                hcc[j * this->nbase_x + i] = cs.zero;
                hcc[i * this->nbase_x + j] = cs.zero;
                scc[j * this->nbase_x + i] = cs.zero;
                scc[i * this->nbase_x + j] = cs.zero;
            }
        }   
    } 
    else
    {
        for (size_t i = 0; i < nbase; i++)
        {
            if (i >= nb1)
            {
                hcc[i * this->nbase_x + i] = set_real_tocomplex(hcc[i * this->nbase_x + i]);
                scc[i * this->nbase_x + i] = set_real_tocomplex(scc[i * this->nbase_x + i]);
            }
            for (size_t j = std::max(i + 1, (size_t)nb1); j < nbase; j++)
            {
                hcc[j * this->nbase_x + i] = get_conj(hcc[i * this->nbase_x + j]);
                scc[j * this->nbase_x + i] = get_conj(scc[i * this->nbase_x + j]);
            }
        }   
    }

    // if (!init)
    // {
    //     std::cout << "output hcc: " << std::endl;
    //     for (size_t i = 0; i < this->nbase_x; i++)
    //     {
    //         for (size_t j = 0; j < this->nbase_x; j++)
    //         {
    //             std::cout << hcc[i * this->nbase_x + j] << "\t";
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << std::endl;

    //     std::cout << "output scc: " << std::endl;
    //     for (size_t i = 0; i < this->nbase_x; i++)
    //     {
    //         for (size_t j = 0; j < this->nbase_x; j++)
    //         {
    //             std::cout << scc[i * this->nbase_x + j] << "\t";
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << std::endl;
    // }

    ModuleBase::timer::tick("Diago_NewDav", "cal_elem");
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
void Diago_NewDav<T, Device>::diag_zhegvx(const int& nbase,
                                          const int& nband,
                                          T* hcc,
                                          T* scc,
                                          const int& nbase_x,
                                          Real* eigenvalue, // in CPU
                                          T* vcc,
                                          bool init)
{
    //	ModuleBase::TITLE("Diago_NewDav","diag_zhegvx");
    ModuleBase::timer::tick("Diago_NewDav", "diag_zhegvx");
    if (GlobalV::RANK_IN_POOL == 0)
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
            if (init)
            {
                dnevx_op<T, Device>()(this->ctx, nbase, this->nbase_x, this->hcc, nband, this->eigenvalue, this->vcc);
            }
            else
            {

                dngvx_op<T, Device>()(this->ctx,
                                      nbase,
                                      this->nbase_x,
                                      this->hcc,
                                      this->scc,
                                      nband,
                                      this->eigenvalue,
                                      this->vcc);
                // dngvd_op<T, Device>()(this->ctx,
                //                       nbase,
                //                       this->nbase_x,
                //                       this->hcc,
                //                       this->scc,
                //                       //   nband,
                //                       this->eigenvalue,
                //                       this->vcc);
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

    // std::cout << "after diago output hcc: " << std::endl;
    // for (size_t i = 0; i < this->nbase_x; i++)
    // {
    //     for (size_t j = 0; j < this->nbase_x; j++)
    //     {
    //         std::cout << hcc[i * this->nbase_x + j] << "\t";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    ModuleBase::timer::tick("Diago_NewDav", "diag_zhegvx");
    return;
}

template <typename T, typename Device>
void Diago_NewDav<T, Device>::refresh(const int& dim,
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
        ModuleBase::TITLE("Diago_NewDav", "refresh");
    ModuleBase::timer::tick("Diago_NewDav", "refresh");

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
            scc[i * this->nbase_x + i] = this->one[0];
            vcc[i * this->nbase_x + i] = this->one[0];
        }
    }
    ModuleBase::timer::tick("Diago_NewDav", "refresh");

    // cancan hcc
    // gemm_op<T, Device>()(this->ctx,
    //                     'C',
    //                     'N',
    //                     this->nbase_x,
    //                     this->nbase_x,
    //                     this->dim,
    //                     this->one,
    //                     hphi,
    //                     this->dim,
    //                     &basis(0, 0),
    //                     this->dim,
    //                     this->zero,
    //                     hcc,
    //                     this->nbase_x);

    // std::cout << "refresh hcc: " << std::endl;
    // for (size_t i = 0; i < this->nbase_x; i++)
    // {
    //     for (size_t j = 0; j < this->nbase_x; j++)
    //     {
    //         std::cout << hcc[i * this->nbase_x + j] << "\t";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    return;
}

/*

// template <typename T, typename Device>
// void Diago_NewDav<T, Device>::SchmitOrth(const int& dim,
//                                          const int nband,
//                                          const int m,
//                                          psi::Psi<T, Device>& basis,
//                                          const T* sphi,
//                                          T* lagrange_m,
//                                          const int mm_size,
//                                          const int mv_size)
// {
//     //	if(test_david == 1) ModuleBase::TITLE("Diago_NewDav","SchmitOrth");
//     ModuleBase::timer::tick("Diago_NewDav", "SchmitOrth");

//     // orthogonalize starting eigenfunction to those already calculated
//     // psi_m orthogonalize to psi(0) ~ psi(m-1)
//     // Attention, the orthogonalize here read as
//     // psi(m) -> psi(m) - \sum_{i < m} \langle psi(i)|S|psi(m) \rangle psi(i)
//     // so the orthogonalize is performed about S.

//     assert(basis.get_nbands() >= nband);
//     assert(m >= 0);
//     assert(m < nband);

//     T* psi_m = &basis(m, 0);

//     // std::complex<double> *lagrange = new std::complex<double>[m + 1];
//     // ModuleBase::GlobalFunc::ZEROS(lagrange, m + 1);

//     // calculate the square matrix for future lagranges
//     if (mm_size != 0)
//     {
//         //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//         // haozhihan repalce 2022-10-16
//         gemm_op<T, Device>()(this->ctx,
//                              'C',
//                              'N',
//                              mm_size,                                // m: row of A,C
//                              mm_size,                                // n: col of B,C
//                              this->dim,                              // k: col of A, row of B
//                              this->one,                              // alpha
//                              &basis(m - mv_size + 1 - mm_size, 0),   // A
//                              this->dim,                              // LDA: if(N) max(1,m) if(T) max(1,k)
//                              &sphi[m * this->dim],                   // B
//                              this->dim,                              // LDB: if(N) max(1,k) if(T) max(1,n)
//                              this->zero,                             // belta
//                              &lagrange_m[m - mv_size + 1 - mm_size], // C
//                              nband                                   // LDC: if(N) max(1, m)
//         );
//     }
//     // calculate other lagranges for this band
//     //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//     //  haozhihan repalce 2022-10-16
//     gemv_op<T, Device>()(this->ctx,
//                          'C',
//                          this->dim,
//                          mv_size,
//                          this->one,
//                          &basis(m - mv_size + 1, 0),
//                          this->dim,
//                          &sphi[m * this->dim],
//                          1,
//                          this->zero,
//                          &lagrange_m[m - mv_size + 1],
//                          1);
//     //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

//     Parallel_Reduce::reduce_pool(lagrange_m, m + 1);

//     T var = *this->zero;
//     syncmem_d2h_op()(this->cpu_ctx, this->ctx, &var, lagrange_m + m, 1);
//     double psi_norm = get_real(var);

//     assert(psi_norm > 0.0);

//     // haozhihan replace 2022-10-24
//     gemv_op<T, Device>()(this->ctx,
//                          'N',
//                          this->dim,
//                          m,
//                          this->neg_one,
//                          &basis(0, 0),
//                          this->dim,
//                          lagrange_m,
//                          1,
//                          this->one,
//                          psi_m,
//                          1);

//     psi_norm -= dot_real_op<T, Device>()(this->ctx, m, lagrange_m, lagrange_m, false);

//     // for (int j = 0; j < m; j++)
//     // {
//     //     const std::complex<double> alpha = std::complex<double>(-1, 0) * lagrange_m[j];
//     //     zaxpy_(&npw, &alpha, &psi(j,0), &inc, psi_m, &inc);
//     //     /*for (int ig = 0; ig < npw; ig++)
//     //     {
//     //         psi_m[ig] -= lagrange[j] * psi(j, ig);
//     //     }*/
//     //     psi_norm -= (conj(lagrange_m[j]) * lagrange_m[j]).real();
//     // }

//     assert(psi_norm > 0.0);

//     psi_norm = sqrt(psi_norm);

//     if (psi_norm < 1.0e-12)
//     {
//         std::cout << "Diago_NewDav::SchmitOrth:aborted for psi_norm <1.0e-12" << std::endl;
//         std::cout << "nband = " << nband << std::endl;
//         std::cout << "m = " << m << std::endl;
//         exit(0);
//     }
//     else
//     {
//         //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//         // haozhihan repalce 2022-10-16
//         vector_div_constant_op<T, Device>()(this->ctx, this->dim, psi_m, psi_m, psi_norm);
//         //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//         // for (int i = 0; i < npw; i++)
//         // {
//         //     psi_m[i] /= psi_norm;
//         // }
//     }

//     // delete[] lagrange;
//     ModuleBase::timer::tick("Diago_NewDav", "SchmitOrth");
//     return;
// }
/*
 */
// template <typename T, typename Device>
// void Diago_NewDav<T, Device>::planSchmitOrth(const int nband, int* pre_matrix_mm_m, int* pre_matrix_mv_m)
// {
//     if (nband <= 0)
//         return;
//     ModuleBase::GlobalFunc::ZEROS(pre_matrix_mm_m, nband);
//     ModuleBase::GlobalFunc::ZEROS(pre_matrix_mv_m, nband);
//     int last_matrix_size = nband;
//     int matrix_size = int(nband / 2);
//     int divide_times = 0;
//     std::vector<int> divide_points(nband);
//     int res_nband = nband - matrix_size;
//     while (matrix_size > 1)
//     {
//         int index = nband - matrix_size;
//         if (divide_times == 0)
//         {
//             divide_points[0] = index;
//             pre_matrix_mm_m[index] = matrix_size;
//             if (res_nband == matrix_size)
//                 pre_matrix_mv_m[index] = 1;
//             else
//                 pre_matrix_mv_m[index] = 2;
//             divide_times = 1;
//         }
//         else
//         {
//             for (int i = divide_times - 1; i >= 0; i--)
//             {
//                 divide_points[i * 2] = divide_points[i] - matrix_size;
//                 divide_points[i * 2 + 1] = divide_points[i * 2] + last_matrix_size;
//                 pre_matrix_mm_m[divide_points[i * 2]] = matrix_size;
//                 pre_matrix_mm_m[divide_points[i * 2 + 1]] = matrix_size;
//                 if (res_nband == matrix_size)
//                 {
//                     pre_matrix_mv_m[divide_points[i * 2]] = 1;
//                     pre_matrix_mv_m[divide_points[i * 2 + 1]] = 1;
//                 }
//                 else
//                 {
//                     pre_matrix_mv_m[divide_points[i * 2]] = 2;
//                     pre_matrix_mv_m[divide_points[i * 2 + 1]] = 2;
//                 }
//             }
//             divide_times *= 2;
//         }
//         last_matrix_size = matrix_size;
//         matrix_size = int(res_nband / 2);
//         res_nband -= matrix_size;
//     }
//     // fill the pre_matrix_mv_m array
//     pre_matrix_mv_m[0] = 1;
//     for (int m = 1; m < nband; m++)
//     {
//         if (pre_matrix_mv_m[m] == 0)
//         {
//             pre_matrix_mv_m[m] = pre_matrix_mv_m[m - 1] + 1;
//         }
//     }
// }

template <typename T, typename Device>
void Diago_NewDav<T, Device>::diag(hamilt::Hamilt<T, Device>* phm_in,
                                   psi::Psi<T, Device>& psi,
                                   Real* eigenvalue_in,
                                   bool need_subspace)
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

        // if (need_subspace || ntry > 0)
        {
            DiagoIterAssist<T, Device>::diagH_subspace(phm_in, psi, psi, eigenvalue_in, psi.get_nbands());
        }

        this->diag_mock(phm_in, psi, eigenvalue_in);

        ++ntry;

    } while (DiagoIterAssist<T, Device>::test_exit_cond(ntry, this->notconv));

    // std::cout << "in new_dav, ntry == " << ntry << std::endl;

    if (notconv > std::max(5, psi.get_nbands() / 4))
    {
        std::cout << "\n notconv = " << this->notconv;
        std::cout << "\n Diago_NewDav::diag', too many bands are not converged! \n";
    }
    return;
}

namespace hsolver
{
template class Diago_NewDav<std::complex<float>, psi::DEVICE_CPU>;
template class Diago_NewDav<std::complex<double>, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Diago_NewDav<std::complex<float>, psi::DEVICE_GPU>;
template class Diago_NewDav<std::complex<double>, psi::DEVICE_GPU>;
#endif
#ifdef __LCAO
template class Diago_NewDav<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Diago_NewDav<double, psi::DEVICE_GPU>;
#endif
#endif
} // namespace hsolver