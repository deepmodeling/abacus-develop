#include <module_hsolver/diago_cg.h>

#include <module_base/memory.h>
#include <module_base/parallel_reduce.h>
#include <module_base/timer.h>
#include <module_base/constants.h>

#include <ATen/kernels/lapack.h>
#include <ATen/kernels/memory.h>
#include <ATen/core/tensor_map.h>
#include <ATen/core/tensor_utils.h>

#include <ATen/ops/einsum_op.h>
#include <ATen/ops/linalg_op.h>

using namespace hsolver;

template<typename T, typename Device>
DiagoCG<T, Device>::DiagoCG(
    const std::string& basis_type,
    const std::string& calculation)
{
    basis_type_ = basis_type;
    calculation_ = calculation;
    this->one_ = new T(static_cast<T>(1.0));
    this->zero_ = new T(static_cast<T>(0.0));
    this->neg_one_ = new T(static_cast<T>(-1.0));
}

template<typename T, typename Device>
DiagoCG<T, Device>::DiagoCG(
    const std::string& basis_type,
    const std::string& calculation,
    const bool& need_subspace,
    const Func& subspace_func,
    const Real& pw_diag_thr,
    const int& pw_diag_nmax,
    const int& nproc_in_pool)
{
    basis_type_ = basis_type;
    calculation_ = calculation;
    need_subspace_ = need_subspace;
    subspace_func_ = subspace_func;
    pw_diag_thr_ = pw_diag_thr;
    pw_diag_nmax_ = pw_diag_nmax;
    nproc_in_pool_ = nproc_in_pool;
    this->one_ = new T(static_cast<T>(1.0));
    this->zero_ = new T(static_cast<T>(0.0));
    this->neg_one_ = new T(static_cast<T>(-1.0));
}

template<typename T, typename Device>
DiagoCG<T, Device>::~DiagoCG()
{
    delete this->one_;
    delete this->zero_;
    delete this->neg_one_;
}

template<typename T, typename Device>
void DiagoCG<T, Device>::diag_mock(const ct::Tensor& prec_in, ct::Tensor& psi, ct::Tensor& eigen)
{
    ModuleBase::TITLE("DiagoCG", "diag_once");
    ModuleBase::timer::tick("DiagoCG", "diag_once");
    
    /// out : record for states of convergence
    this->notconv_ = 0;
    /// initialize variables
    this->n_band_ = psi.shape().dim_size(0);
    this->n_basis_ = psi.shape().dim_size(1);

    /// record for how many loops in cg convergence
    int avg = 0;
    //-------------------------------------------------------------------
    // "poor man" iterative diagonalization of a complex hermitian matrix
    // through preconditioned conjugate grad algorithm
    // Band-by-band algorithm with minimal use of memory
    // Calls hPhi and sPhi to calculate H|phi> and S|phi>
    // Works for generalized eigenvalue problem (US pseudopotentials) as well
    //-------------------------------------------------------------------
    // phi_m = new psi::Psi<T, Device>(phi, 1, 1);
    auto phi_m = std::move(ct::Tensor(
        ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, {this->n_basis_}));
    // hphi.resize(this->n_basis_max_, ModuleBase::ZERO);
    auto hphi = std::move(ct::Tensor(
        ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, {this->n_basis_}));
    // sphi.resize(this->n_basis_max_, ModuleBase::ZERO);
    auto sphi = std::move(ct::Tensor(
        ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, {this->n_basis_}));
    // pphi.resize(this->n_basis_max_, ModuleBase::ZERO);
    auto pphi = std::move(ct::Tensor(
        ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, {this->n_basis_}));

    // cg = new psi::Psi<T, Device>(phi, 1, 1);
    auto cg = std::move(ct::Tensor(
        ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, {this->n_basis_}));
    // scg.resize(this->n_basis_max_, ModuleBase::ZERO);
    auto scg = std::move(ct::Tensor(
        ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, {this->n_basis_}));

    // grad.resize(this->n_basis_max_, ModuleBase::ZERO);
    auto grad = std::move(ct::Tensor(
        ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, {this->n_basis_}));
    // g0.resize(this->n_basis_max_, ModuleBase::ZERO);
    auto g0 = std::move(ct::Tensor(
        ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, {this->n_basis_}));
    // lagrange.resize(this->n_band, ModuleBase::ZERO);
    auto lagrange = std::move(ct::Tensor(
        ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, {this->n_band_}));
    
    auto prec = prec_in;
    if (prec.NumElements() == 0) {
        prec = ct::Tensor(
            ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, {this->n_basis_});
        prec.set_value(static_cast<Real>(1.0));
    }

    ModuleBase::Memory::record("DiagoCG", this->n_basis_ * 10);

    eigen.zero();
    auto eigen_pack = eigen.accessor<Real, 1>();
    for (int m = 0; m < this->n_band_; m++)
    {
        phi_m.sync(psi[m]);
        //copy psi_in into internal psi, m=0 has been done in Constructor
        this->spsi_func_(phi_m, sphi); // sphi = S|psi(m)>
        this->schmit_orth(m, psi, sphi, phi_m);
        this->spsi_func_(phi_m, sphi); // sphi = S|psi(m)>
        this->hpsi_func_(phi_m, hphi); // hphi = H|psi(m)>

        eigen_pack[m]  = 
            dot_real_op()(ctx_, this->n_basis_, phi_m.data<T>(), hphi.data<T>());
        
        int  iter      = 0;
        Real gg_last   = 0.0;
        Real cg_norm   = 0.0;
        Real theta     = 0.0;
        bool converged = false;

        do 
        {
            this->calc_grad(prec, grad, hphi, sphi, pphi);
            this->orth_grad(psi, m, grad, scg, lagrange);
            this->calc_gamma_cg(
                iter,  // const int&
                cg_norm, theta, // const Real&
                prec, scg, grad, phi_m, // const Tensor&
                gg_last, // Real&
                g0, cg); // Tensor&
            
            this->hpsi_func_(cg, pphi);
            this->spsi_func_(cg, scg);
            
            converged = this->update_psi(
                pphi, cg, scg, // const Tensor&
                cg_norm, theta, eigen_pack[m], // Real&
                phi_m, sphi, hphi); // Tensor&

        } while (!converged && ++iter < pw_diag_nmax_);

        psi[m].sync(phi_m);
        if (!converged) {
            ++this->notconv_;
        }
        avg += static_cast<Real>(iter) + 1.00;

        // reorder eigenvalue if they are not in the right order
        // (this CAN and WILL happen in not-so-special cases)
        if (m > 0) {
            ModuleBase::GlobalFunc::NOTE("reorder bands!");
            if (eigen_pack[m] - eigen_pack[m - 1] < -2.0 * pw_diag_thr_) {
                // if the last calculated eigenvalue is not the largest...
                int ii = 0;
                for (ii = m - 2; ii >= 0; ii--) {
                    if (eigen_pack[m] - eigen_pack[ii] > 2.0 * pw_diag_thr_)
                        break;
                }
                ii++;
                // last calculated eigenvalue should be in the ii-th position: reorder
                Real e0 = eigen_pack[m];
                // ModuleBase::GlobalFunc::COPYARRAY(psi_temp, pphi, this->n_basis_);
                pphi.sync(psi[m]);

                for (int jj = m; jj >= ii + 1; jj--) {
                    eigen_pack[jj] = eigen_pack[jj - 1];
                    psi[jj].sync(psi[jj - 1]);
                }
                eigen_pack[ii] = e0;
                psi[ii].sync(pphi);
            } // endif
        } // end reorder
    } // end m

    avg /= this->n_band_;
    avg_iter_ += avg;

    ModuleBase::timer::tick("DiagoCG", "diag_once");
} // end subroutine ccgdiagg

template<typename T, typename Device>
void DiagoCG<T, Device>::calc_grad(
    const ct::Tensor& prec,
    ct::Tensor& grad,
    ct::Tensor& hphi,
    ct::Tensor& sphi,
    ct::Tensor& pphi)
{
    // for (int i = 0; i < this->n_basis_; i++)
    // {
    //     //(2) PH|psi>
    //     grad.data<T>()[i] = this->hphi[i] / this->precondition[i];
    //     //(3) PS|psi>
    //     this->pphi[i] = this->sphi[i] / this->precondition[i];
    // }
    // denghui replace this at 20221106
    // TODO: use GPU precondition to initialize CG class
    vector_div_vector_op<T, Device>()(ctx_, this->n_basis_, grad.data<T>(), hphi.data<T>(), prec.data<Real>());
    vector_div_vector_op<T, Device>()(ctx_, this->n_basis_, pphi.data<T>(), sphi.data<T>(), prec.data<Real>());

    // Update lambda !
    // (4) <psi|SPH|psi >
    const Real eh = hsolver::dot_real_op<T, Device>()(ctx_, this->n_basis_, sphi.data<T>(), grad.data<T>());
    // (5) <psi|SPS|psi >
    const Real es = hsolver::dot_real_op<T, Device>()(ctx_, this->n_basis_, sphi.data<T>(), pphi.data<T>());
    const Real lambda = eh / es;

    // Update g!
    // for (int i = 0; i < this->n_basis_; i++)
    // {
    //     //               <psi|SPH|psi>
    //     // (6) PH|psi> - ------------- * PS |psi>
    //     //               <psi|SPS|psi>
    //     //
    //     // So here we get the gradient.
    //     grad.data<T>()[i] -= lambda * this->pphi[i];
    // }
    // haozhihan replace this 2022-10-6
    constantvector_addORsub_constantVector_op<T, Device>()(ctx_, this->n_basis_, grad.data<T>(), grad.data<T>(), 1.0, pphi.data<T>(), (-lambda));
}

template<typename T, typename Device>
void DiagoCG<T, Device>::orth_grad(
    const ct::Tensor& psi, 
    const int& m, 
    ct::Tensor& grad, 
    ct::Tensor& scg,
    ct::Tensor& lagrange)
{
    this->spsi_func_(grad, scg); // scg = S|grad>
    gemv_op<T, Device>()(
        ctx_,
        'C',
        this->n_basis_,
        m,
        this->one_,
        psi.data<T>(),
        this->n_basis_,
        scg.data<T>(),
        1,
        this->zero_,
        lagrange.data<T>(),
        1);

    Parallel_Reduce::reduce_pool(lagrange.data<T>(), m);

    // (3) orthogonal |g> and |scg> to all states (0~m-1)
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // haozhihan replace 2022-10-07
    gemv_op<T, Device>()(
        ctx_,
        'N',
        this->n_basis_,
        m,
        this->neg_one_,
        psi.data<T>(),
        this->n_basis_,
        lagrange.data<T>(),
        1,
        this->one_,
        grad.data<T>(),
        1);

    gemv_op<T, Device>()(
        ctx_,
        'N',
        this->n_basis_,
        m,
        this->neg_one_,
        psi.data<T>(),
        this->n_basis_,
        lagrange.data<T>(),
        1,
        this->one_,
        scg.data<T>(),
        1);
}

template<typename T, typename Device>
void DiagoCG<T, Device>::calc_gamma_cg(
    const int& iter,
    const Real& cg_norm, 
    const Real& theta,
    const ct::Tensor& prec,
    const ct::Tensor& scg,
    const ct::Tensor& grad,
    const ct::Tensor& phi_m,
    Real& gg_last,
    ct::Tensor& g0,
    ct::Tensor& cg)
{
    Real gg_inter;
    if (iter > 0)
    {
        // (1) Update gg_inter!
        // gg_inter = <g|g0>
        // Attention : the 'g' in g0 is getted last time
        gg_inter
            = hsolver::dot_real_op<T, Device>()(ctx_, this->n_basis_, grad.data<T>(), g0.data<T>()); // b means before
    }

    // (2) Update for g0!
    // two usage:
    // firstly, for now, calculate: gg_now
    // secondly, prepare for the next iteration: gg_inter
    // |g0> = P | scg >
    // for (int i = 0; i < this->n_basis_; i++)
    // {
    //     g0[i] = this->precondition[i] * this->scg[i];
    // }
    // denghui replace this 20221106
    // TODO: use GPU precondition instead
    vector_mul_vector_op<T, Device>()(ctx_, this->n_basis_, g0.data<T>(), scg.data<T>(), prec.data<Real>());

    // (3) Update gg_now!
    // gg_now = < g|P|scg > = < g|g0 >
    const Real gg_now = hsolver::dot_real_op<T, Device>()(ctx_, this->n_basis_, grad.data<T>(), g0.data<T>());

    if (iter == 0)
    {
        // (40) gg_last first value : equal gg_now
        gg_last = gg_now;
        // (50) cg direction first value : |g>
        // |cg> = |g>
        cg.sync(grad);
    }
    else
    {
        // (4) Update gamma !
        REQUIRES_OK(gg_last != 0.0,
            "DiagoCG_New::calc_gamma_cg: gg_last is zero, which is not allowed!");
        const Real gamma = (gg_now - gg_inter) / gg_last;

        // (5) Update gg_last !
        gg_last = gg_now;
        // (6) Update cg direction !(need gamma and |go> ):
        // for (int i = 0; i < this->n_basis_; i++)
        // {
        //     pcg[i] = gamma * pcg[i] + grad.data<T>()[i];
        // }
        // haozhihan replace this 2022-10-6
        constantvector_addORsub_constantVector_op<T, Device>()(ctx_, this->n_basis_, cg.data<T>(), cg.data<T>(), gamma, grad.data<T>(), 1.0);

        const Real norma = gamma * cg_norm * sin(theta);
        T znorma = static_cast<T>(norma * -1);

        // haozhihan replace this 2022-10-6
        // const int one = 1;
        // zaxpy_(&this->n_basis_, &znorma, pphi_m, &one, pcg, &one);
        /*for (int i = 0; i < this->n_basis_; i++)
        {
            pcg[i] -= norma * pphi_m[i];
        }*/
        axpy_op<T, Device>()(ctx_, this->n_basis_, &znorma, phi_m.data<T>(), 1, cg.data<T>(), 1);
    }
}

template<typename T, typename Device>
bool DiagoCG<T, Device>::update_psi(
    const ct::Tensor& pphi,
    const ct::Tensor& cg,
    const ct::Tensor& scg,
    Real &cg_norm, 
    Real &theta, 
    Real &eigen,
    ct::Tensor& phi_m,
    ct::Tensor& sphi,
    ct::Tensor& hphi)
{
    cg_norm = sqrt(hsolver::dot_real_op<T, Device>()(ctx_, this->n_basis_, cg.data<T>(), scg.data<T>()));

    if (cg_norm < 1.0e-10)
        return true;

    const Real a0
        = hsolver::dot_real_op<T, Device>()(ctx_, this->n_basis_, phi_m.data<T>(), pphi.data<T>()) * 2.0 / cg_norm;
    const Real b0
        = hsolver::dot_real_op<T, Device>()(ctx_, this->n_basis_, cg.data<T>(), pphi.data<T>()) / (cg_norm * cg_norm);

    const Real e0 = eigen;
    theta = atan(a0 / (e0 - b0)) / 2.0;

    const Real new_e = (e0 - b0) * cos(2.0 * theta) + a0 * sin(2.0 * theta);

    const Real e1 = (e0 + b0 + new_e) / 2.0;
    const Real e2 = (e0 + b0 - new_e) / 2.0;

    if (e1 > e2) {
        theta += ModuleBase::PI_HALF;
    }

    eigen = std::min(e1, e2);

    const Real cost = cos(theta);
    const Real sint_norm = sin(theta) / cg_norm;

    // for (int i = 0; i < this->n_basis_; i++)
    // {
    //     phi_m_pointer[i] = phi_m_pointer[i] * cost + sint_norm * pcg[i];
    // }
    
    // haozhihan replace this 2022-10-6
    constantvector_addORsub_constantVector_op<T, Device>()(ctx_, this->n_basis_, phi_m.data<T>(), phi_m.data<T>(), cost, cg.data<T>(), sint_norm);

    if (std::abs(eigen - e0) < pw_diag_thr_) {
        // ModuleBase::timer::tick("DiagoCG","update");
        return true;
    }
    else {
        // for (int i = 0; i < this->n_basis_; i++)
        // {
        //     this->sphi[i] = this->sphi[i] * cost + sint_norm * this->scg[i];
        //     this->hphi[i] = this->hphi[i] * cost + sint_norm * this->pphi[i];
        // }

        // haozhihan replace this 2022-10-6
        constantvector_addORsub_constantVector_op<T, Device>()(ctx_, this->n_basis_, sphi.data<T>(), sphi.data<T>(), cost, scg.data<T>(),  sint_norm);
        constantvector_addORsub_constantVector_op<T, Device>()(ctx_, this->n_basis_, hphi.data<T>(), hphi.data<T>(), cost, pphi.data<T>(), sint_norm);
        return false;
    }
}


template<typename T, typename Device>
void DiagoCG<T, Device>::schmit_orth(
    const int& m, 
    const ct::Tensor& psi, 
    const ct::Tensor& sphi, 
    ct::Tensor& phi_m)
{
    //	ModuleBase::TITLE("DiagoCG","schmit_orth");
    // ModuleBase::timer::tick("DiagoCG","schmit_orth");
    // orthogonalize starting eigenfunction to those already calculated
    // phi_m orthogonalize to psi(start) ~ psi(m-1)
    // Attention, the orthogonalize here read as
    // psi(m) -> psi(m) - \sum_{i < m} < psi(i) | S | psi(m) > psi(i)
    // so the orthogonalize is performed about S.
    REQUIRES_OK(m >= 0,
        "DiagoCG_New::schmit_orth: m < 0");
    REQUIRES_OK(this->n_band_ >= m,
        "DiagoCG_New::schmit_orth: n_band < m");

    ct::Tensor lagrange_so = ct::Tensor(
        ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, {m + 1});

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // haozhihan replace 2022-10-6
    int inc = 1;
    gemv_op<T, Device>()(
        ctx_,
        'C',
        this->n_basis_,
        m + 1,
        this->one_,
        psi.data<T>(),
        this->n_basis_,
        sphi.data<T>(),
        inc,
        this->zero_,
        lagrange_so.data<T>(),
        inc);

    // be careful , here reduce m+1
    Parallel_Reduce::reduce_pool(lagrange_so.data<T>(), m + 1);

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // haozhihan replace 2022-10-6
    gemv_op<T, Device>()(
        ctx_,
        'N',
        this->n_basis_,
        m,
        this->neg_one_,
        psi.data<T>(),
        this->n_basis_,
        lagrange_so.data<T>(),
        inc,
        this->one_,
        phi_m.data<T>(),
        inc);

    //======================================================================
    /*for (int j = 0; j < m; j++)
    {
        for (int ig =0; ig < dim; ig++)
        {
            phi_m[ig] -= lagrange[j] * psi(j, ig);
        }
        psi_norm -= ( conj(lagrange[j]) * lagrange[j] ).real();
    }*/
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    auto psi_norm = ct::extract<Real>(lagrange_so[m])
        - dot_real_op()(ctx_, m, lagrange_so.data<T>(), lagrange_so.data<T>(), false);

    if (psi_norm <= 0.0)
    {
        std::cout << " m = " << m << std::endl;
        for (int j = 0; j <= m; ++j)
        {
            std::cout << "j = " << j << " lagrange norm = " 
                      << ct::extract<Real>(lagrange_so[j] * lagrange_so[j])
                      << std::endl;
        }
        std::cout << " in DiagoCG, psi norm = " << psi_norm << std::endl;
        std::cout << " If you use GNU compiler, it may due to the zdotc is unavailable." << std::endl;
        ModuleBase::WARNING_QUIT("schmit_orth", "psi_norm <= 0.0");
    }

    psi_norm = sqrt(psi_norm);

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // haozhihan replace 2022-10-6
    // scal_op<Real, Device>()(ctx_, this->n_basis_, &psi_norm, pphi_m, 1);
    //======================================================================
    // for (int ig = 0; ig < this->n_basis_; ig++)
    // {
    //     pphi_m[ig] /= psi_norm;
    // }
    vector_div_constant_op<T, Device>()(ctx_, this->n_basis_, phi_m.data<T>(), phi_m.data<T>(), psi_norm);

    // ModuleBase::timer::tick("DiagoCG","schmit_orth");
}


template<typename T, typename Device>
bool DiagoCG<T, Device>::test_exit_cond(
    const int& ntry,
    const int& notconv) const
{
    const bool scf = calculation_ != "nscf";
    // If ntry <=5, try to do it better, if ntry > 5, exit.
    const bool f1 = ntry <= 5;
    // In non-self consistent calculation, do until totally converged.
    const bool f2 = !scf && notconv > 0;
    // if self consistent calculation, if not converged > 5,
    // using diagH_subspace and cg method again. ntry++
    const bool f3 = scf && notconv > 5;
    return f1 && (f2 || f3);
}

template<typename T, typename Device>
void DiagoCG<T, Device>::diag(
    const Func& hpsi_func, 
    const Func& spsi_func, 
    ct::Tensor& psi,
    ct::Tensor& eigen,
    const ct::Tensor& prec)
{
    /// record the times of trying iterative diagonalization
    int ntry = 0;
    this->notconv_ = 0;
    hpsi_func_ = hpsi_func;
    spsi_func_ = spsi_func;

    do
    {
        if (need_subspace_ || ntry > 0) {
            this->subspace_func_(psi, psi);
        }

        ++ntry;
        avg_iter_ += 1.0;
        this->diag_mock(prec, psi, eigen);
    } while (this->test_exit_cond(ntry, this->notconv_));

    if (this->notconv_ > std::max(5, this->n_band_ / 4)) {
        std::cout << "\n notconv = " << this->notconv_;
        std::cout << "\n DiagoCG::diag', too many bands are not converged! \n";
    }
}

namespace hsolver {
template class DiagoCG<std::complex<float>, psi::DEVICE_CPU>;
template class DiagoCG<std::complex<double>, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoCG<std::complex<float>, psi::DEVICE_GPU>;
template class DiagoCG<std::complex<double>, psi::DEVICE_GPU>;
#endif 

#ifdef __LCAO
template class DiagoCG<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoCG<double, psi::DEVICE_GPU>;
#endif
#endif
} // namespace hsolver