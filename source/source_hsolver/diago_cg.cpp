#include <ATen/core/tensor_map.h>
#include <ATen/core/tensor_utils.h>
#include <ATen/kernels/lapack.h>
#include <ATen/kernels/memory.h>
#include <ATen/ops/einsum_op.h>
#include <ATen/ops/linalg_op.h>
#include <source_base/constants.h>
#include <source_base/memory.h>
#include <source_base/parallel_reduce.h>
#include <source_base/timer.h>
#include <source_base/tool_title.h>             // ModuleBase::TITLE
#include <source_base/global_function.h>        // ModuleBase::GlobalFunc::NOTE
#include <source_hsolver/diago_cg.h>

using namespace hsolver;

template <typename T, typename Device>
DiagoCG<T, Device>::DiagoCG(const std::string& basis_type, const std::string& calculation)
{
    basis_type_ = basis_type;
    calculation_ = calculation;
    this->one_ = new T(static_cast<T>(1.0));
    this->zero_ = new T(static_cast<T>(0.0));
    this->neg_one_ = new T(static_cast<T>(-1.0));
}

template <typename T, typename Device>
DiagoCG<T, Device>::DiagoCG(const std::string& basis_type,
                            const std::string& calculation,
                            const bool& need_subspace,
                            const SubspaceFunc& subspace_func,
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

template <typename T, typename Device>
DiagoCG<T, Device>::~DiagoCG()
{
    delete this->one_;
    delete this->zero_;
    delete this->neg_one_;
}

// P1: 确保 workspace 大小正确（复用内存）
template <typename T, typename Device>
void DiagoCG<T, Device>::ensure_workspace(int n_basis, int n_band)
{
    auto ensure_size = [&](ct::Tensor& t, const std::vector<int64_t>& shape) {
        if (t.shape().dim_size(0) != shape[0])
        {
            t = ct::Tensor(ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, shape);
        }
    };
    ensure_size(ws_phi_m_,     {n_basis});
    ensure_size(ws_hphi_,      {n_basis});
    ensure_size(ws_sphi_,      {n_basis});
    ensure_size(ws_pphi_,      {n_basis});
    ensure_size(ws_cg_,        {n_basis});
    ensure_size(ws_scg_,       {n_basis});
    ensure_size(ws_grad_,      {n_basis});
    ensure_size(ws_g0_,        {n_basis});
    ensure_size(ws_temp_,      {n_basis});
    ensure_size(ws_lagrange_,  {n_band});
}

template <typename T, typename Device>
void DiagoCG<T, Device>::diag_once(const ct::Tensor& prec_in,
                                   ct::Tensor& psi,
                                   ct::Tensor& eigen,
                                   const std::vector<double>& ethr_band)
{
    ModuleBase::TITLE("DiagoCG", "diag_once");
    // ModuleBase::timer::tick("DiagoCG", "diag_once");  // disabled for compilation

    this->notconv_ = 0;
    this->n_band_ = psi.shape().dim_size(0);
    this->n_basis_ = psi.shape().dim_size(1);

    int avg = 0;

    // P1: 使用成员 workspace，一次性 resize
    // ModuleBase::timer::tick("DiagoCG", "setup_workspace");
    ensure_workspace(this->n_basis_, this->n_band_);
    // ModuleBase::timer::tick("DiagoCG", "setup_workspace");

    auto& phi_m = ws_phi_m_;
    auto& hphi  = ws_hphi_;
    auto& sphi  = ws_sphi_;
    auto& pphi  = ws_pphi_;
    auto& cg    = ws_cg_;
    auto& scg   = ws_scg_;
    auto& grad  = ws_grad_;
    auto& g0    = ws_g0_;
    auto& lagrange = ws_lagrange_;

    // ModuleBase::timer::tick("DiagoCG", "setup_prec");
    auto prec = prec_in;
    if (prec.NumElements() == 0)
    {
        prec = ct::Tensor(ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, {this->n_basis_});
        prec.set_value(static_cast<Real>(1.0));
    }
    // ModuleBase::timer::tick("DiagoCG", "setup_prec");

    ModuleBase::Memory::record("DiagoCG", this->n_basis_ * 10);

    eigen.zero();
    auto eigen_pack = eigen.accessor<Real, 1>();

    for (int m = 0; m < this->n_band_; m++)
    {
        // ModuleBase::timer::tick("DiagoCG", "band_copy");
        phi_m.sync(psi[m]);
        // ModuleBase::timer::tick("DiagoCG", "band_copy");

        // ModuleBase::timer::tick("DiagoCG", "spsi_phi_before_orth");
        this->spsi_func_(phi_m.data<T>(), sphi.data<T>(), this->n_basis_, 1);
        // ModuleBase::timer::tick("DiagoCG", "spsi_phi_before_orth");
        this->schmit_orth(m, psi, sphi, phi_m);

        // ModuleBase::timer::tick("DiagoCG", "spsi_phi_after_orth");
        this->spsi_func_(phi_m.data<T>(), sphi.data<T>(), this->n_basis_, 1);
        // ModuleBase::timer::tick("DiagoCG", "spsi_phi_after_orth");
        // ModuleBase::timer::tick("DiagoCG", "hpsi_phi");
        this->hpsi_func_(phi_m.data<T>(), hphi.data<T>(), this->n_basis_, 1);
        // ModuleBase::timer::tick("DiagoCG", "hpsi_phi");

        // ModuleBase::timer::tick("DiagoCG", "eigen_dot");
        eigen_pack[m] = dot_real_op()(this->n_basis_, phi_m.data<T>(), hphi.data<T>());
        // ModuleBase::timer::tick("DiagoCG", "eigen_dot");

        int iter = 0;
        Real gg_last = 0.0;
        Real cg_norm = 0.0;
        Real theta = 0.0;
        bool converged = false;

        do
        {
            this->calc_grad(prec, grad, hphi, sphi, pphi);
            this->orth_grad(psi, m, grad, scg, lagrange);
            this->calc_gamma_cg(iter, cg_norm, theta, prec, scg, grad, phi_m, gg_last, g0, cg);

            // ModuleBase::timer::tick("DiagoCG", "hpsi_cg");
            this->hpsi_func_(cg.data<T>(), pphi.data<T>(), this->n_basis_, 1);
            // ModuleBase::timer::tick("DiagoCG", "hpsi_cg");
            // ModuleBase::timer::tick("DiagoCG", "spsi_cg");
            this->spsi_func_(cg.data<T>(), scg.data<T>(), this->n_basis_, 1);
            // ModuleBase::timer::tick("DiagoCG", "spsi_cg");

            converged = this->update_psi(pphi, cg, scg, ethr_band[m], cg_norm, theta, eigen_pack[m], phi_m, sphi, hphi);

        } while (!converged && ++iter < pw_diag_nmax_);

        // ModuleBase::timer::tick("DiagoCG", "band_save");
        psi[m].sync(phi_m);
        // ModuleBase::timer::tick("DiagoCG", "band_save");

        if (!converged) ++this->notconv_;
        iter_band.push_back(iter);
        avg += static_cast<Real>(iter) + 1.00;

        if (m > 0)
        {
            // ModuleBase::timer::tick("DiagoCG", "reorder_check");
            ModuleBase::GlobalFunc::NOTE("reorder bands!");
            if (eigen_pack[m] - eigen_pack[m - 1] < -2.0 * pw_diag_thr_)
            {
                int ii = 0;
                for (ii = m - 2; ii >= 0; ii--)
                    if (eigen_pack[m] - eigen_pack[ii] > 2.0 * pw_diag_thr_) break;
                ii++;
                Real e0 = eigen_pack[m];
                pphi.sync(psi[m]);
                for (int jj = m; jj >= ii + 1; jj--)
                {
                    eigen_pack[jj] = eigen_pack[jj - 1];
                    psi[jj].sync(psi[jj - 1]);
                }
                eigen_pack[ii] = e0;
                psi[ii].sync(pphi);
            }
            // ModuleBase::timer::tick("DiagoCG", "reorder_check");
        }
    }

    avg /= this->n_band_;
    avg_iter_ += avg;

    // ModuleBase::timer::tick("DiagoCG", "diag_once");
}

// P1: 融合除法与点积，减少 kernel 调用
template <typename T, typename Device>
void DiagoCG<T, Device>::calc_grad(const ct::Tensor& prec,
                                   ct::Tensor& grad,
                                   ct::Tensor& hphi,
                                   ct::Tensor& sphi,
                                   ct::Tensor& pphi)
{
    // ModuleBase::timer::tick("DiagoCG", "calc_grad");

    const auto* prec_data = prec.data<Real>();
    auto* grad_data  = grad.data<T>();
    auto* hphi_data  = hphi.data<T>();
    auto* sphi_data  = sphi.data<T>();
    auto* pphi_data  = pphi.data<T>();

    Real eh = 0.0, es = 0.0;

    for (int i = 0; i < this->n_basis_; ++i)
    {
        T h_val = hphi_data[i] / prec_data[i];
        T s_val = sphi_data[i] / prec_data[i];
        grad_data[i] = h_val;
        pphi_data[i] = s_val;
        eh += std::real(std::conj(sphi_data[i]) * h_val);
        es += std::real(std::conj(sphi_data[i]) * s_val);
    }

    const Real lambda = eh / es;

    for (int i = 0; i < this->n_basis_; ++i)
    {
        grad_data[i] -= lambda * pphi_data[i];
    }

    // ModuleBase::timer::tick("DiagoCG", "calc_grad");
}

// P1: 用临时向量投影代替两次 gemv 更新
template <typename T, typename Device>
void DiagoCG<T, Device>::orth_grad(const ct::Tensor& psi,
                                   const int& m,
                                   ct::Tensor& grad,
                                   ct::Tensor& scg,
                                   ct::Tensor& lagrange)
{
    // ModuleBase::timer::tick("DiagoCG", "orth_grad");

    // ModuleBase::timer::tick("DiagoCG", "orth_grad_spsi");
    this->spsi_func_(grad.data<T>(), scg.data<T>(), this->n_basis_, 1);
    // ModuleBase::timer::tick("DiagoCG", "orth_grad_spsi");

    // ModuleBase::timer::tick("DiagoCG", "orth_grad_lagrange");
    ModuleBase::gemv_op<T, Device>()('C', this->n_basis_, m,
                                     this->one_, psi.data<T>(), this->n_basis_,
                                     scg.data<T>(), 1,
                                     this->zero_, lagrange.data<T>(), 1);
    // ModuleBase::timer::tick("DiagoCG", "orth_grad_lagrange");

    // ModuleBase::timer::tick("DiagoCG", "orth_grad_reduce");
    Parallel_Reduce::reduce_pool(lagrange.data<T>(), m);
    // ModuleBase::timer::tick("DiagoCG", "orth_grad_reduce");

    // 计算投影向量 proj = psi * lagrange
    // ModuleBase::timer::tick("DiagoCG", "orth_grad_project");
    ModuleBase::gemv_op<T, Device>()('N', this->n_basis_, m,
                                     this->one_, psi.data<T>(), this->n_basis_,
                                     lagrange.data<T>(), 1,
                                     this->zero_, ws_temp_.data<T>(), 1);
    // ModuleBase::timer::tick("DiagoCG", "orth_grad_project");

    const T neg_one = static_cast<T>(-1.0);
    // ModuleBase::timer::tick("DiagoCG", "orth_grad_update_grad");
    ModuleBase::axpy_op<T, Device>()(this->n_basis_, &neg_one, ws_temp_.data<T>(), 1, grad.data<T>(), 1);
    // ModuleBase::timer::tick("DiagoCG", "orth_grad_update_grad");

    // ModuleBase::timer::tick("DiagoCG", "orth_grad_update_scg");
    ModuleBase::axpy_op<T, Device>()(this->n_basis_, &neg_one, ws_temp_.data<T>(), 1, scg.data<T>(), 1);
    // ModuleBase::timer::tick("DiagoCG", "orth_grad_update_scg");

    // ModuleBase::timer::tick("DiagoCG", "orth_grad");
}

template <typename T, typename Device>
void DiagoCG<T, Device>::calc_gamma_cg(const int& iter,
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
    // ModuleBase::timer::tick("DiagoCG", "calc_gamma_cg");
    Real gg_inter;
    if (iter > 0)
    {
        // ModuleBase::timer::tick("DiagoCG", "calc_gamma_inter_dot");
        gg_inter = ModuleBase::dot_real_op<T, Device>()(this->n_basis_, grad.data<T>(), g0.data<T>());
        // ModuleBase::timer::tick("DiagoCG", "calc_gamma_inter_dot");
    }

    // ModuleBase::timer::tick("DiagoCG", "calc_gamma_precond");
    ModuleBase::vector_mul_vector_op<T, Device>()(this->n_basis_, g0.data<T>(), scg.data<T>(), prec.data<Real>());
    // ModuleBase::timer::tick("DiagoCG", "calc_gamma_precond");

    // ModuleBase::timer::tick("DiagoCG", "calc_gamma_now_dot");
    const Real gg_now = ModuleBase::dot_real_op<T, Device>()(this->n_basis_, grad.data<T>(), g0.data<T>());
    // ModuleBase::timer::tick("DiagoCG", "calc_gamma_now_dot");

    if (iter == 0)
    {
        // ModuleBase::timer::tick("DiagoCG", "calc_gamma_init_dir");
        gg_last = gg_now;
        cg.sync(grad);
        // ModuleBase::timer::tick("DiagoCG", "calc_gamma_init_dir");
    }
    else
    {
        // ModuleBase::timer::tick("DiagoCG", "calc_gamma_update_dir");
        REQUIRES_OK(gg_last != 0.0, "DiagoCG_New::calc_gamma_cg: gg_last is zero, which is not allowed!");
        const Real gamma = (gg_now - gg_inter) / gg_last;
        gg_last = gg_now;
        ModuleBase::vector_add_vector_op<T, Device>()(this->n_basis_, cg.data<T>(), cg.data<T>(), gamma, grad.data<T>(), 1.0);

        const Real norma = gamma * cg_norm * sin(theta);
        T znorma = static_cast<T>(norma * -1);
        ModuleBase::axpy_op<T, Device>()(this->n_basis_, &znorma, phi_m.data<T>(), 1, cg.data<T>(), 1);
        // ModuleBase::timer::tick("DiagoCG", "calc_gamma_update_dir");
    }
    // ModuleBase::timer::tick("DiagoCG", "calc_gamma_cg");
}

// P1: 融合 sphi 和 hphi 的更新循环
template <typename T, typename Device>
bool DiagoCG<T, Device>::update_psi(const ct::Tensor& pphi,
                                    const ct::Tensor& cg,
                                    const ct::Tensor& scg,
                                    const double& ethreshold,
                                    Real& cg_norm,
                                    Real& theta,
                                    Real& eigen,
                                    ct::Tensor& phi_m,
                                    ct::Tensor& sphi,
                                    ct::Tensor& hphi)
{
    // ModuleBase::timer::tick("DiagoCG", "update_psi");
    // ModuleBase::timer::tick("DiagoCG", "update_norm");
    cg_norm = sqrt(ModuleBase::dot_real_op<T, Device>()(this->n_basis_, cg.data<T>(), scg.data<T>()));
    // ModuleBase::timer::tick("DiagoCG", "update_norm");

    if (cg_norm < 1.0e-10)
    {
        // ModuleBase::timer::tick("DiagoCG", "update_psi");
        return true;
    }

    // ModuleBase::timer::tick("DiagoCG", "update_theta");
    const Real a0 = ModuleBase::dot_real_op<T, Device>()(this->n_basis_, phi_m.data<T>(), pphi.data<T>()) * 2.0 / cg_norm;
    const Real b0 = ModuleBase::dot_real_op<T, Device>()(this->n_basis_, cg.data<T>(), pphi.data<T>()) / (cg_norm * cg_norm);
    const Real e0 = eigen;
    theta = atan(a0 / (e0 - b0)) / 2.0;
    const Real new_e = (e0 - b0) * cos(2.0 * theta) + a0 * sin(2.0 * theta);
    const Real e1 = (e0 + b0 + new_e) / 2.0;
    const Real e2 = (e0 + b0 - new_e) / 2.0;
    if (e1 > e2) theta += ModuleBase::PI_HALF;
    eigen = std::min(e1, e2);
    // ModuleBase::timer::tick("DiagoCG", "update_theta");

    const Real cost = cos(theta);
    const Real sint_norm = sin(theta) / cg_norm;

    // ModuleBase::timer::tick("DiagoCG", "update_phi");
    ModuleBase::vector_add_vector_op<T, Device>()(this->n_basis_, phi_m.data<T>(), phi_m.data<T>(), cost, cg.data<T>(), sint_norm);
    // ModuleBase::timer::tick("DiagoCG", "update_phi");

    if (std::abs(eigen - e0) < ethreshold)
    {
        // ModuleBase::timer::tick("DiagoCG", "update_psi");
        return true;
    }
    else
    {
        // 一次循环同时更新 sphi 和 hphi
        // ModuleBase::timer::tick("DiagoCG", "update_sh");
        auto* sphi_data = sphi.data<T>();
        auto* scg_data  = scg.data<T>();
        auto* hphi_data = hphi.data<T>();
        auto* pphi_data = pphi.data<T>();
        for (int i = 0; i < this->n_basis_; ++i)
        {
            sphi_data[i] = cost * sphi_data[i] + sint_norm * scg_data[i];
            hphi_data[i] = cost * hphi_data[i] + sint_norm * pphi_data[i];
        }
        // ModuleBase::timer::tick("DiagoCG", "update_sh");
        // ModuleBase::timer::tick("DiagoCG", "update_psi");
        return false;
    }
}

// P1: 使用成员 ws_lagrange_，避免 schmit_orth 内部重复分配
template <typename T, typename Device>
void DiagoCG<T, Device>::schmit_orth(const int& m, const ct::Tensor& psi, const ct::Tensor& sphi, ct::Tensor& phi_m)
{
    // ModuleBase::timer::tick("DiagoCG","schmit_orth");
    REQUIRES_OK(m >= 0, "DiagoCG_New::schmit_orth: m < 0");
    REQUIRES_OK(this->n_band_ >= m, "DiagoCG_New::schmit_orth: n_band < m");

    T* lagrange_so = ws_lagrange_.data<T>();
    int inc = 1;

    // ModuleBase::timer::tick("DiagoCG","schmit_lagrange");
    ModuleBase::gemv_op<T, Device>()('C', this->n_basis_, m + 1,
                                     this->one_, psi.data<T>(), this->n_basis_,
                                     sphi.data<T>(), inc,
                                     this->zero_, lagrange_so, inc);
    // ModuleBase::timer::tick("DiagoCG","schmit_lagrange");

    // ModuleBase::timer::tick("DiagoCG","schmit_reduce");
    Parallel_Reduce::reduce_pool(lagrange_so, m + 1);
    // ModuleBase::timer::tick("DiagoCG","schmit_reduce");

    // ModuleBase::timer::tick("DiagoCG","schmit_project");
    ModuleBase::gemv_op<T, Device>()('N', this->n_basis_, m,
                                     this->neg_one_, psi.data<T>(), this->n_basis_,
                                     lagrange_so, inc,
                                     this->one_, phi_m.data<T>(), inc);
    // ModuleBase::timer::tick("DiagoCG","schmit_project");

    Real psi_norm = ct::extract<Real>(lagrange_so[m])
                    - dot_real_op()(m, lagrange_so, lagrange_so, false);

    if (psi_norm <= 0.0)
    {
        std::cout << " m = " << m << std::endl;
        for (int j = 0; j <= m; ++j)
            std::cout << "j = " << j << " lagrange norm = " << ct::extract<Real>(lagrange_so[j] * lagrange_so[j]) << std::endl;
        std::cout << " in DiagoCG, psi norm = " << psi_norm << std::endl;
        // ModuleBase::timer::tick("DiagoCG","schmit_orth");
        ModuleBase::WARNING_QUIT("schmit_orth", "psi_norm <= 0.0");
    }

    psi_norm = sqrt(psi_norm);
    // ModuleBase::timer::tick("DiagoCG","schmit_normalize");
    ModuleBase::vector_mul_real_op<T, Device>()(this->n_basis_, phi_m.data<T>(), phi_m.data<T>(), Real(1.0 / psi_norm));
    // ModuleBase::timer::tick("DiagoCG","schmit_normalize");

    // ModuleBase::timer::tick("DiagoCG","schmit_orth");
}

template <typename T, typename Device>
bool DiagoCG<T, Device>::test_exit_cond(const int& ntry, const int& notconv) const
{
    const bool scf = calculation_ != "nscf";
    const bool f1 = ntry <= 5;
    const bool f2 = !scf && notconv > 0;
    const bool f3 = scf && notconv > 5;
    return f1 && (f2 || f3);
}

template <typename T, typename Device>
double DiagoCG<T, Device>::diag(const HPsiFunc& hpsi_func,
                                const SPsiFunc& spsi_func,
                                const int ld_psi,
                                const int nband,
                                const int dim,
                                T* psi_in,
                                Real* eigenvalue_in,
                                const std::vector<double>& ethr_band,
                                const Real* prec)
{
    REQUIRES_OK(ld_psi >= dim, "DiagoCG::diag: ld_psi must be >= dim");
    REQUIRES_OK(static_cast<int>(ethr_band.size()) >= nband, "DiagoCG::diag: ethr_band size must be >= nband");

    auto psi = ct::TensorMap(psi_in, ct::DataTypeToEnum<T>::value, ct::DeviceTypeToEnum<ct_Device>::value, ct::TensorShape({nband, ld_psi}));
    auto eigen = ct::TensorMap(eigenvalue_in, ct::DataTypeToEnum<Real>::value, ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value, ct::TensorShape({nband}));

    ct::Tensor prec_tensor;
    if (prec != nullptr)
    {
        prec_tensor = ct::TensorMap(const_cast<Real*>(prec), ct::DataTypeToEnum<Real>::value, ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value, ct::TensorShape({dim}))
                          .template to_device<ct_Device>();
    }

    int ntry = 0;
    this->notconv_ = 0;
    hpsi_func_ = hpsi_func;
    spsi_func_ = spsi_func;

    ct::Tensor psi_temp = psi.slice({0, 0}, {nband, dim});
    do
    {
        if (ntry > 0)
        {
            ct::TensorMap psi_map = ct::TensorMap(psi.data(), psi_temp);
            const bool assume_S_orthogonal = true;
            this->subspace_func_(psi_temp.data<T>(), psi_map.data<T>(), dim, nband, assume_S_orthogonal);
            psi_temp.sync(psi_map);
        }
        else if (need_subspace_)
        {
            ct::TensorMap psi_map = ct::TensorMap(psi.data(), psi_temp);
            const bool assume_S_orthogonal = false;
            this->subspace_func_(psi_temp.data<T>(), psi_map.data<T>(), dim, nband, assume_S_orthogonal);
            psi_temp.sync(psi_map);
        }

        ++ntry;
        avg_iter_ += 1.0;
        this->diag_once(prec_tensor, psi_temp, eigen, ethr_band);
    } while (this->test_exit_cond(ntry, this->notconv_));

    if (this->notconv_ > std::max(5, this->n_band_ / 4))
    {
        std::cout << "\n notconv = " << this->notconv_;
        std::cout << "\n DiagoCG::diag', too many bands are not converged! \n";
    }
    psi.zero();
    psi.sync(psi_temp);

    return avg_iter_;
}

// 显式实例化
namespace hsolver
{
template class DiagoCG<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoCG<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoCG<std::complex<float>, base_device::DEVICE_GPU>;
template class DiagoCG<std::complex<double>, base_device::DEVICE_GPU>;
#endif
#ifdef __LCAO
template class DiagoCG<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoCG<double, base_device::DEVICE_GPU>;
#endif
#endif
}