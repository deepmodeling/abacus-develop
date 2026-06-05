#include <source_hsolver/diago_cg_mixed.h>

#include <algorithm>
#include <cmath>

#include <ATen/core/tensor_map.h>
#include <ATen/core/tensor_utils.h>
#include <ATen/kernels/memory.h>
#include <source_base/constants.h>
#include <source_base/memory_recorder.h>
#include <source_base/parallel_reduce.h>
#include <source_base/timer.h>
#include <source_base/tool_title.h>
#include <source_base/global_function.h>

namespace hsolver {

template <typename T, typename Device>
DiagoCGMixed<T, Device>::DiagoCGMixed(const std::string& basis_type, const std::string& calculation)
{
    basis_type_ = basis_type; calculation_ = calculation;
}

template <typename T, typename Device>
DiagoCGMixed<T, Device>::DiagoCGMixed(const std::string& basis_type, const std::string& calculation,
                                       const bool& need_subspace, const Func& subspace_func,
                                       const Real& pw_diag_thr, const int& pw_diag_nmax, const int& nproc_in_pool)
{
    basis_type_ = basis_type; calculation_ = calculation;
    need_subspace_ = need_subspace; subspace_func_ = subspace_func;
    pw_diag_thr_ = pw_diag_thr; pw_diag_nmax_ = pw_diag_nmax; nproc_in_pool_ = nproc_in_pool;
}

template <typename T, typename Device>
void DiagoCGMixed<T, Device>::convert_d2f(const ct::Tensor& d_src, ct::Tensor& f_dst)
{
    const int n = d_src.NumElements();
    const T* d = d_src.data<T>();
    T_float* f = f_dst.data<T_float>();
    for (int i = 0; i < n; i++) f[i] = static_cast<T_float>(d[i]);
}

template <typename T, typename Device>
void DiagoCGMixed<T, Device>::convert_f2d(const ct::Tensor& f_src, ct::Tensor& d_dst)
{
    const int n = f_src.NumElements();
    const T_float* f = f_src.data<T_float>();
    T* d = d_dst.data<T>();
    for (int i = 0; i < n; i++) d[i] = static_cast<T>(f[i]);
}

template <typename T, typename Device>
void DiagoCGMixed<T, Device>::diag_once(const ct::Tensor& prec_in, ct::Tensor& psi,
                                         ct::Tensor& eigen, const std::vector<double>& ethr_band)
{
    ModuleBase::TITLE("DiagoCGMixed", "diag_once");
    ModuleBase::timer::tick("DiagoCGMixed", "diag_once");

    notconv_ = 0;
    n_band_ = psi.shape().dim_size(0);
    n_basis_ = psi.shape().dim_size(1);
    int avg = 0;

    auto dt = ct::DataTypeToEnum<T>::value;
    auto dtf = ct::DataTypeToEnum<T_float>::value;
    auto dev = ct::DeviceTypeToEnum<ct_Device>::value;
    auto dtr = ct::DataTypeToEnum<Real>::value;
    auto dtfr = ct::DataTypeToEnum<Real_float>::value;

    auto phi_m = ct::Tensor(dt, dev, {n_basis_}), hphi = ct::Tensor(dt, dev, {n_basis_});
    auto sphi = ct::Tensor(dt, dev, {n_basis_}), pphi = ct::Tensor(dt, dev, {n_basis_});
    auto cg = ct::Tensor(dt, dev, {n_basis_}), scg = ct::Tensor(dt, dev, {n_basis_});
    auto grad = ct::Tensor(dt, dev, {n_basis_}), g0 = ct::Tensor(dt, dev, {n_basis_});
    auto lagrange = ct::Tensor(dt, dev, {n_band_});

    auto phi_m_f = ct::Tensor(dtf, dev, {n_basis_}), hphi_f = ct::Tensor(dtf, dev, {n_basis_});
    auto sphi_f = ct::Tensor(dtf, dev, {n_basis_});

    auto prec = prec_in;
    if (prec.NumElements() == 0) {
        prec = ct::Tensor(ct::DataTypeToEnum<Real>::value, ct::DeviceTypeToEnum<ct_Device>::value, {n_basis_});
        prec.set_value(static_cast<Real>(1.0));
    }
    auto prec_f = ct::Tensor(ct::DataTypeToEnum<Real_float>::value, ct::DeviceTypeToEnum<ct_Device>::value, {n_basis_});
    {
        auto* p = prec.data<Real>();
        auto* q = prec_f.data<Real_float>();
        for (int i = 0; i < n_basis_; i++) q[i] = static_cast<Real_float>(p[i]);
    }

    ModuleBase::Memory::record("DiagoCGMixed", n_basis_ * 12);
    eigen.zero();
    auto eig = eigen.accessor<Real, 1>();

    for (int m = 0; m < n_band_; m++)
    {
        phi_m.sync(psi[m]);
        spsi_func_(phi_m, sphi);
        schmit_orth(m, psi, sphi, phi_m);

        convert_d2f(phi_m, phi_m_f);
        hpsi_func_(phi_m_f, hphi_f);
        spsi_func_(phi_m_f, sphi_f);
        convert_f2d(hphi_f, hphi);
        convert_f2d(sphi_f, sphi);

        eig[m] = ModuleBase::dot_real_op<T, Device>()(n_basis_, phi_m.data<T>(), hphi.data<T>());

        int iter = 0;
        Real gg_last = 0, cg_norm = 0, theta = 0;
        bool converged = false;

        do {
            calc_grad(prec_f, grad, hphi, sphi, pphi);
            orth_grad(psi, m, grad, scg, lagrange);
            calc_gamma_cg(iter, cg_norm, theta, prec_f, scg, grad, phi_m, gg_last, g0, cg);

            {
                auto cg_f = ct::Tensor(dtf, dev, {n_basis_}), pphi_f = ct::Tensor(dtf, dev, {n_basis_}), scg_f = ct::Tensor(dtf, dev, {n_basis_});
                convert_d2f(cg, cg_f);
                hpsi_func_(cg_f, pphi_f);
                spsi_func_(cg_f, scg_f);
                convert_f2d(pphi_f, pphi);
                convert_f2d(scg_f, scg);
            }

            converged = update_psi(pphi, cg, scg, ethr_band[m], cg_norm, theta, eig[m], phi_m, sphi, hphi);
        } while (!converged && ++iter < pw_diag_nmax_);

        psi[m].sync(phi_m);
        if (!converged) ++notconv_;
        avg += static_cast<Real>(iter) + 1;

        if (m > 0 && eig[m] - eig[m - 1] < -2.0 * pw_diag_thr_) {
            int ii = m - 2;
            while (ii >= 0 && eig[m] - eig[ii] <= 2.0 * pw_diag_thr_) ii--;
            ii++;
            Real e0 = eig[m]; pphi.sync(psi[m]);
            for (int jj = m; jj > ii; jj--) { eig[jj] = eig[jj - 1]; psi[jj].sync(psi[jj - 1]); }
            eig[ii] = e0; psi[ii].sync(pphi);
        }
    }
    avg /= n_band_; avg_iter_ += avg;
    ModuleBase::timer::tick("DiagoCGMixed", "diag_once");
}

template <typename T, typename Device>
void DiagoCGMixed<T, Device>::calc_grad(const ct::Tensor& prec_f, ct::Tensor& grad,
                                         ct::Tensor& hphi, ct::Tensor& sphi, ct::Tensor& pphi)
{
    auto dtf = ct::DataTypeToEnum<T_float>::value;
    auto dev = ct::DeviceTypeToEnum<ct_Device>::value;
    auto hphi_f = ct::Tensor(dtf, dev, {n_basis_}), sphi_f = ct::Tensor(dtf, dev, {n_basis_});
    auto grad_f = ct::Tensor(dtf, dev, {n_basis_}), pphi_f = ct::Tensor(dtf, dev, {n_basis_});
    convert_d2f(hphi, hphi_f);
    convert_d2f(sphi, sphi_f);

    ModuleBase::vector_div_vector_op<T_float, Device>()(n_basis_, grad_f.data<T_float>(), hphi_f.data<T_float>(), prec_f.data<Real_float>());
    ModuleBase::vector_div_vector_op<T_float, Device>()(n_basis_, pphi_f.data<T_float>(), sphi_f.data<T_float>(), prec_f.data<Real_float>());

    convert_f2d(grad_f, grad);
    convert_f2d(pphi_f, pphi);

    const Real eh = ModuleBase::dot_real_op<T, Device>()(n_basis_, sphi.data<T>(), grad.data<T>());
    const Real es = ModuleBase::dot_real_op<T, Device>()(n_basis_, sphi.data<T>(), pphi.data<T>());
    ModuleBase::vector_add_vector_op<T, Device>()(n_basis_, grad.data<T>(), grad.data<T>(), 1.0, pphi.data<T>(), -(eh / es));
}

template <typename T, typename Device>
void DiagoCGMixed<T, Device>::orth_grad(const ct::Tensor& psi, const int& m,
                                         ct::Tensor& grad, ct::Tensor& scg, ct::Tensor& lagrange)
{
    const T one(1.0), zero(0.0), neg_one(-1.0);
    spsi_func_(grad, scg);

    ModuleBase::gemv_op<T, Device>()('C', n_basis_, m, &one, psi.data<T>(), n_basis_, scg.data<T>(), 1, &zero, lagrange.data<T>(), 1);
    Parallel_Reduce::reduce_pool(lagrange.data<T>(), m);

    ModuleBase::gemv_op<T, Device>()('N', n_basis_, m, &neg_one, psi.data<T>(), n_basis_, lagrange.data<T>(), 1, &one, grad.data<T>(), 1);
    ModuleBase::gemv_op<T, Device>()('N', n_basis_, m, &neg_one, psi.data<T>(), n_basis_, lagrange.data<T>(), 1, &one, scg.data<T>(), 1);
}

template <typename T, typename Device>
void DiagoCGMixed<T, Device>::calc_gamma_cg(const int& iter, const Real& cg_norm, const Real& theta,
                                             const ct::Tensor& prec_f, const ct::Tensor& scg,
                                             const ct::Tensor& grad, const ct::Tensor& phi_m,
                                             Real& gg_last, ct::Tensor& g0, ct::Tensor& cg)
{
    Real gg_inter = 0.0;
    if (iter > 0) gg_inter = ModuleBase::dot_real_op<T, Device>()(n_basis_, grad.data<T>(), g0.data<T>());

    auto dtf = ct::DataTypeToEnum<T_float>::value;
    auto dev = ct::DeviceTypeToEnum<ct_Device>::value;
    auto scg_f = ct::Tensor(dtf, dev, {n_basis_}), g0_f = ct::Tensor(dtf, dev, {n_basis_});
    convert_d2f(scg, scg_f);
    ModuleBase::vector_mul_vector_op<T_float, Device>()(n_basis_, g0_f.data<T_float>(), scg_f.data<T_float>(), prec_f.data<Real_float>());
    convert_f2d(g0_f, g0);

    const Real gg_now = ModuleBase::dot_real_op<T, Device>()(n_basis_, grad.data<T>(), g0.data<T>());

    if (iter == 0) {
        gg_last = gg_now;
        cg.sync(grad);
    } else {
        const Real gamma = (gg_now - gg_inter) / gg_last;
        gg_last = gg_now;
        ModuleBase::vector_add_vector_op<T, Device>()(n_basis_, cg.data<T>(), cg.data<T>(), gamma, grad.data<T>(), 1.0);
        T znorma = static_cast<T>(-gamma * cg_norm * std::sin(theta));
        ModuleBase::axpy_op<T, Device>()(n_basis_, &znorma, phi_m.data<T>(), 1, cg.data<T>(), 1);
    }
}

template <typename T, typename Device>
bool DiagoCGMixed<T, Device>::update_psi(const ct::Tensor& pphi, const ct::Tensor& cg, const ct::Tensor& scg,
                                          const double& ethreshold, Real& cg_norm, Real& theta, Real& eigen,
                                          ct::Tensor& phi_m, ct::Tensor& sphi, ct::Tensor& hphi)
{
    cg_norm = std::sqrt(ModuleBase::dot_real_op<T, Device>()(n_basis_, cg.data<T>(), scg.data<T>()));
    if (cg_norm < 1e-10) return true;

    const Real a0 = ModuleBase::dot_real_op<T, Device>()(n_basis_, phi_m.data<T>(), pphi.data<T>()) * 2.0 / cg_norm;
    const Real b0 = ModuleBase::dot_real_op<T, Device>()(n_basis_, cg.data<T>(), pphi.data<T>()) / (cg_norm * cg_norm);
    const Real e0 = eigen;
    const Real denom = e0 - b0;
    // Guard against near-degenerate case where e0 ≈ b0
    if (std::abs(denom) < 1e-30) {
        // For degenerate case, diagonalize directly
        const Real disc = std::sqrt(a0 * a0);
        eigen = std::min(e0 + b0 - disc, e0 + b0 + disc) / 2.0;
        theta = (a0 > 0 ? 1.0 : -1.0) * ModuleBase::PI / 4.0;
    } else {
        theta = std::atan(a0 / denom) / 2.0;
        const Real new_e = denom * std::cos(2.0 * theta) + a0 * std::sin(2.0 * theta);
        const Real e1 = (e0 + b0 + new_e) / 2.0, e2 = (e0 + b0 - new_e) / 2.0;
        if (e1 > e2) theta += ModuleBase::PI_HALF;
        eigen = std::min(e1, e2);
    }

    const Real cost = std::cos(theta), sint_norm = std::sin(theta) / cg_norm;
    ModuleBase::vector_add_vector_op<T, Device>()(n_basis_, phi_m.data<T>(), phi_m.data<T>(), cost, cg.data<T>(), sint_norm);

    if (std::abs(eigen - e0) < ethreshold) return true;

    ModuleBase::vector_add_vector_op<T, Device>()(n_basis_, sphi.data<T>(), sphi.data<T>(), cost, scg.data<T>(), sint_norm);
    ModuleBase::vector_add_vector_op<T, Device>()(n_basis_, hphi.data<T>(), hphi.data<T>(), cost, pphi.data<T>(), sint_norm);
    return false;
}

template <typename T, typename Device>
void DiagoCGMixed<T, Device>::schmit_orth(const int& m, const ct::Tensor& psi,
                                           const ct::Tensor& sphi, ct::Tensor& phi_m)
{
    const T one(1.0), zero(0.0), neg_one(-1.0);
    ct::Tensor lagrange_so(ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, {m + 1});

    ModuleBase::gemv_op<T, Device>()('C', n_basis_, m + 1, &one, psi.data<T>(), n_basis_, sphi.data<T>(), 1, &zero, lagrange_so.data<T>(), 1);
    Parallel_Reduce::reduce_pool(lagrange_so.data<T>(), m + 1);
    ModuleBase::gemv_op<T, Device>()('N', n_basis_, m, &neg_one, psi.data<T>(), n_basis_, lagrange_so.data<T>(), 1, &one, phi_m.data<T>(), 1);

    auto psi_norm = ct::extract<Real>(lagrange_so[m])
                    - ModuleBase::dot_real_op<T, Device>()(m, lagrange_so.data<T>(), lagrange_so.data<T>(), false);
    ModuleBase::vector_div_constant_op<T, Device>()(n_basis_, phi_m.data<T>(), phi_m.data<T>(), std::sqrt(psi_norm));
}

template <typename T, typename Device>
bool DiagoCGMixed<T, Device>::test_exit_cond(const int& ntry, const int& notconv) const
{
    return ntry <= 5 && (calculation_ != "nscf" ? notconv > 5 : notconv > 0);
}

template <typename T, typename Device>
void DiagoCGMixed<T, Device>::diag(const Func& hpsi_func, const Func& spsi_func,
                                    ct::Tensor& psi, ct::Tensor& eigen,
                                    const std::vector<double>& ethr_band, const ct::Tensor& prec)
{
    int ntry = 0; notconv_ = 0;
    hpsi_func_ = hpsi_func; spsi_func_ = spsi_func;
    const int nbasis = int(psi.shape().dim_size(1));
    const int nprec = prec.NumElements() > 0 ? int(prec.shape().dim_size(0)) : nbasis;
    ct::Tensor psi_temp = psi.slice({0, 0}, {int(psi.shape().dim_size(0)), nprec});
    do {
        if (need_subspace_ || ntry > 0) {
            ct::TensorMap psi_map = ct::TensorMap(psi.data(), psi_temp);
            subspace_func_(psi_temp, psi_map);
            psi_temp.sync(psi_map);
        }
        ++ntry; avg_iter_ += 1.0;
        diag_once(prec, psi_temp, eigen, ethr_band);
    } while (test_exit_cond(ntry, notconv_));
    psi.zero(); psi.sync(psi_temp);
}

template class DiagoCGMixed<std::complex<double>, base_device::DEVICE_CPU>;
template class DiagoCGMixed<std::complex<float>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoCGMixed<std::complex<double>, base_device::DEVICE_GPU>;
template class DiagoCGMixed<std::complex<float>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
