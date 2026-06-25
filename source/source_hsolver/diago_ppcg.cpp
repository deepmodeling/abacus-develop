#include "source_hsolver/diago_ppcg.h"

#include "diago_iter_assist.h"
#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/memory_recorder.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"

#include <ATen/core/tensor_map.h>
#include <ATen/core/tensor_utils.h>
#include <ATen/kernels/lapack.h>
#include <ATen/kernels/memory.h>
#include <cmath>
#include <limits>

namespace hsolver {

template <typename T, typename Device>
DiagoPPCG<T, Device>::DiagoPPCG(const std::string& basis_type, const std::string& calculation)
{
    basis_type_ = basis_type;
    calculation_ = calculation;
    this->one_ = new T(static_cast<T>(1.0));
    this->zero_ = new T(static_cast<T>(0.0));
    this->neg_one_ = new T(static_cast<T>(-1.0));
}

template <typename T, typename Device>
DiagoPPCG<T, Device>::DiagoPPCG(const std::string& basis_type,
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
DiagoPPCG<T, Device>::~DiagoPPCG()
{
    delete this->one_;
    delete this->zero_;
    delete this->neg_one_;
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::compute_gradient(const ct::Tensor& prec,
                                            const ct::Tensor& hpsi,
                                            const ct::Tensor& spsi,
                                            const Real& eigenvalue,
                                            ct::Tensor& grad,
                                            ct::Tensor& ppsi)
{
    // Following ABACUS DiagoCG::calc_grad pattern:
    // grad = P * H|ψ>
    ModuleBase::vector_div_vector_op<T, Device>()(this->n_basis_, grad.data<T>(), hpsi.data<T>(), prec.data<Real>());
    // ppsi = P * S|ψ>
    ModuleBase::vector_div_vector_op<T, Device>()(this->n_basis_, ppsi.data<T>(), spsi.data<T>(), prec.data<Real>());

    // λ_c = <Sψ|P*Hψ> / <Sψ|P*Sψ>
    Real eh = ModuleBase::dot_real_op<T, Device>()(this->n_basis_, spsi.data<T>(), grad.data<T>());
    Real es = ModuleBase::dot_real_op<T, Device>()(this->n_basis_, spsi.data<T>(), ppsi.data<T>());

    // g = P*H|ψ> - λ_c * P*S|ψ>  (gradient orthogonal to S|ψ>)
    Real lambda_c = (es > std::numeric_limits<Real>::epsilon()) ? eh / es : 0.0;
    ModuleBase::vector_add_vector_op<T, Device>()(this->n_basis_,
                                                  grad.data<T>(),
                                                  grad.data<T>(),
                                                  1.0,
                                                  ppsi.data<T>(),
                                                  (-lambda_c));
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::orth_gradient(const ct::Tensor& psi,
                                         const int& m,
                                         ct::Tensor& grad,
                                         ct::Tensor& s_grad,
                                         ct::Tensor& lagrange)
{
    // s_grad = S|grad>
    this->spsi_func_(grad.data<T>(), s_grad.data<T>(), this->n_basis_, 1);

    // lagrange[i] = <ψ_i|S|grad> for i = 0..m-1
    ModuleBase::gemv_op<T, Device>()('C',
                                     this->n_basis_,
                                     m,
                                     this->one_,
                                     psi.data<T>(),
                                     this->ld_psi_,
                                     s_grad.data<T>(),
                                     1,
                                     this->zero_,
                                     lagrange.data<T>(),
                                     1);

    // All-reduce over MPI pools
    Parallel_Reduce::reduce_pool(lagrange.data<T>(), m);

    // |grad> -= Σ ψ_i * lagrange_i
    ModuleBase::gemv_op<T, Device>()('N',
                                     this->n_basis_,
                                     m,
                                     this->neg_one_,
                                     psi.data<T>(),
                                     this->ld_psi_,
                                     lagrange.data<T>(),
                                     1,
                                     this->one_,
                                     grad.data<T>(),
                                     1);

    // |s_grad> -= Σ ψ_i * lagrange_i
    ModuleBase::gemv_op<T, Device>()('N',
                                     this->n_basis_,
                                     m,
                                     this->neg_one_,
                                     psi.data<T>(),
                                     this->ld_psi_,
                                     lagrange.data<T>(),
                                     1,
                                     this->one_,
                                     s_grad.data<T>(),
                                     1);
}

template <typename T, typename Device>
typename DiagoPPCG<T, Device>::Real DiagoPPCG<T, Device>::update_cg_direction(
    ct::Tensor& cg,
    ct::Tensor& scg,
    const ct::Tensor& grad,
    const ct::Tensor& s_grad,
    const ct::Tensor& prec,
    ct::Tensor& g0,
    Real& gg_last,
    const int& iter)
{
    // Following ABACUS DiagoCG::calc_gamma_cg pattern:
    // g0 is persistent between iterations
    // First compute gg_inter using g0 from previous iteration
    Real gg_inter = 0.0;
    if (iter > 0)
    {
        gg_inter = ModuleBase::dot_real_op<T, Device>()(this->n_basis_, grad.data<T>(), g0.data<T>());
    }

    // g0 = P * S|grad> = prec * scg (element-wise multiply)
    ModuleBase::vector_mul_vector_op<T, Device>()(this->n_basis_,
                                                  g0.data<T>(),
                                                  s_grad.data<T>(),
                                                  prec.data<Real>());

    const Real gg_now = ModuleBase::dot_real_op<T, Device>()(this->n_basis_, grad.data<T>(), g0.data<T>());

    Real gamma = 0.0;
    if (iter == 0)
    {
        gg_last = gg_now;
        // cg = grad
        cg.sync(grad);
        scg.sync(s_grad);
    }
    else
    {
        // Polak-Ribiere: gamma = (gg_now - gg_inter) / gg_last
        if (gg_last > std::numeric_limits<Real>::epsilon())
        {
            gamma = (gg_now - gg_inter) / gg_last;
            if (gamma < 0.0) gamma = 0.0; // Safety: restart if negative
        }
        gg_last = gg_now;

        // cg = gamma * cg + grad
        ModuleBase::vector_add_vector_op<T, Device>()(this->n_basis_,
                                                      cg.data<T>(),
                                                      cg.data<T>(),
                                                      gamma,
                                                      grad.data<T>(),
                                                      1.0);
        // scg = gamma * scg + s_grad
        ModuleBase::vector_add_vector_op<T, Device>()(this->n_basis_,
                                                      scg.data<T>(),
                                                      scg.data<T>(),
                                                      gamma,
                                                      s_grad.data<T>(),
                                                      1.0);
    }
    return gamma;
}

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::line_minimization(const ct::Tensor& ppsi,
                                             const ct::Tensor& cg,
                                             const ct::Tensor& scg,
                                             const double& ethreshold,
                                             Real& cg_norm,
                                             Real& theta,
                                             Real& eigenvalue,
                                             ct::Tensor& psi_m,
                                             ct::Tensor& spsi_m,
                                             ct::Tensor& hpsi_m)
{
    // Compute <cg|S|cg> for normalization
    cg_norm = std::sqrt(ModuleBase::dot_real_op<T, Device>()(this->n_basis_, cg.data<T>(), scg.data<T>()));

    if (cg_norm < 1.0e-12)
    {
        return true;
    }

    // Compute Rayleigh-Ritz coefficients in the 2D subspace [ψ, cg]
    // a0 = 2 * <ψ|H|cg> / ||cg||_S
    // b0 = <cg|H|cg> / ||cg||_S^2
    // e0 = <ψ|H|ψ> (current eigenvalue)
    const Real a0 = ModuleBase::dot_real_op<T, Device>()(this->n_basis_, psi_m.data<T>(), ppsi.data<T>()) * 2.0 / cg_norm;
    const Real b0 = ModuleBase::dot_real_op<T, Device>()(this->n_basis_, cg.data<T>(), ppsi.data<T>()) / (cg_norm * cg_norm);

    const Real e0 = eigenvalue;
    // Use atan matching ABACUS CG exactly: theta = atan(a0 / (e0 - b0)) / 2
    const Real denom = e0 - b0;
    theta = (std::abs(denom) > std::numeric_limits<Real>::epsilon()) 
            ? std::atan(a0 / denom) / 2.0 
            : 0.0;

    const Real new_e = (e0 - b0) * std::cos(2.0 * theta) + a0 * std::sin(2.0 * theta);
    const Real e1 = (e0 + b0 + new_e) / 2.0;
    const Real e2 = (e0 + b0 - new_e) / 2.0;

    // Select the smaller eigenvalue
    if (e1 > e2)
    {
        theta += ModuleBase::PI_HALF;
    }
    eigenvalue = std::min(e1, e2);

    const Real cost = std::cos(theta);
    const Real sint_norm = std::sin(theta) / cg_norm;

    // ψ_new = ψ * cos(θ) + cg * sin(θ)/||cg||
    ModuleBase::vector_add_vector_op<T, Device>()(this->n_basis_,
                                                  psi_m.data<T>(),
                                                  psi_m.data<T>(),
                                                  cost,
                                                  cg.data<T>(),
                                                  sint_norm);

    if (std::abs(eigenvalue - e0) < ethreshold)
    {
        return true;
    }
    else
    {
        // Update S|ψ> and H|ψ> for the next iteration
        ModuleBase::vector_add_vector_op<T, Device>()(this->n_basis_,
                                                      spsi_m.data<T>(),
                                                      spsi_m.data<T>(),
                                                      cost,
                                                      scg.data<T>(),
                                                      sint_norm);
        ModuleBase::vector_add_vector_op<T, Device>()(this->n_basis_,
                                                      hpsi_m.data<T>(),
                                                      hpsi_m.data<T>(),
                                                      cost,
                                                      ppsi.data<T>(),
                                                      sint_norm);
        return false;
    }
}

template <typename T, typename Device>
void DiagoPPCG<T, Device>::schmidt_orth(const int& m,
                                        const ct::Tensor& psi,
                                        const ct::Tensor& spsi_m,
                                        ct::Tensor& psi_m)
{
    REQUIRES_OK(m >= 0, "DiagoPPCG::schmidt_orth: m < 0");
    REQUIRES_OK(this->n_band_ >= m, "DiagoPPCG::schmidt_orth: n_band < m");

    ct::Tensor lagrange(ct::DataTypeToEnum<T>::value,
                        ct::DeviceTypeToEnum<ct_Device>::value,
                        {m + 1});

    // lagrange[i] = <ψ_i|S|ψ_m> for i = 0..m
    ModuleBase::gemv_op<T, Device>()('C',
                                     this->n_basis_,
                                     m + 1,
                                     this->one_,
                                     psi.data<T>(),
                                     this->ld_psi_,
                                     spsi_m.data<T>(),
                                     1,
                                     this->zero_,
                                     lagrange.data<T>(),
                                     1);

    Parallel_Reduce::reduce_pool(lagrange.data<T>(), m + 1);

    // Gram-Schmidt: ψ_m -= Σ ψ_i * lagrange_i  (i = 0..m-1)
    ModuleBase::gemv_op<T, Device>()('N',
                                     this->n_basis_,
                                     m,
                                     this->neg_one_,
                                     psi.data<T>(),
                                     this->ld_psi_,
                                     lagrange.data<T>(),
                                     1,
                                     this->one_,
                                     psi_m.data<T>(),
                                     1);

    // Compute new norm and renormalize
    Real psi_norm = ModuleBase::dot_real_op<T, Device>()(this->n_basis_, psi_m.data<T>(), psi_m.data<T>());
    if (psi_norm <= 0.0)
    {
        std::cout << " DiagoPPCG: psi_norm <= 0.0, m = " << m << std::endl;
        ModuleBase::WARNING_QUIT("schmidt_orth", "psi_norm <= 0.0");
    }
    psi_norm = std::sqrt(psi_norm);

    // Normalize
    ModuleBase::vector_mul_real_op<T, Device>()(this->n_basis_, psi_m.data<T>(), psi_m.data<T>(), Real(1.0 / psi_norm));
}

template <typename T, typename Device>
bool DiagoPPCG<T, Device>::test_exit_cond(const int& ntry, const int& notconv) const
{
    const bool scf = calculation_ != "nscf";
    const bool f1 = ntry <= 5;
    const bool f2 = !scf && notconv > 0;
    const bool f3 = scf && notconv > 5;
    return f1 && (f2 || f3);
}

template <typename T, typename Device>
double DiagoPPCG<T, Device>::diag(const HPsiFunc& hpsi_func,
                                  const SPsiFunc& spsi_func,
                                  const int ld_psi,
                                  const int nband,
                                  const int dim,
                                  T* psi_in,
                                  Real* eigenvalue_in,
                                  const std::vector<double>& ethr_band,
                                  const Real* prec)
{
    ModuleBase::TITLE("DiagoPPCG", "diag");
    ModuleBase::timer::start("DiagoPPCG", "diag");

    REQUIRES_OK(ld_psi >= dim, "DiagoPPCG::diag: ld_psi must be >= dim");
    REQUIRES_OK(static_cast<int>(ethr_band.size()) >= nband,
                "DiagoPPCG::diag: ethr_band size must be >= nband");

    // Set up tensor views
    auto psi = ct::TensorMap(psi_in,
                             ct::DataTypeToEnum<T>::value,
                             ct::DeviceTypeToEnum<ct_Device>::value,
                             ct::TensorShape({nband, ld_psi}));
    auto eigen = ct::TensorMap(eigenvalue_in,
                               ct::DataTypeToEnum<Real>::value,
                               ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value,
                               ct::TensorShape({nband}));

    // Store parameters
    this->n_band_ = nband;
    this->n_basis_ = dim;
    this->ld_psi_ = ld_psi;
    this->notconv_ = 0;
    this->avg_iter_ = 0.0;

    this->hpsi_func_ = hpsi_func;
    this->spsi_func_ = spsi_func;

    // Use provided convergence parameters (from constructor) or fall back to DiagoIterAssist
    if (this->pw_diag_nmax_ == 0)
    {
        this->pw_diag_nmax_ = DiagoIterAssist<T, Device>::PW_DIAG_NMAX;
    }
    if (this->pw_diag_thr_ == 1e-5)
    {
        this->pw_diag_thr_ = DiagoIterAssist<T, Device>::PW_DIAG_THR;
    }

    // Allocate working arrays
    // psi_m: current band vector
    auto psi_m = ct::Tensor(ct::DataTypeToEnum<T>::value,
                            ct::DeviceTypeToEnum<ct_Device>::value,
                            {this->n_basis_});
    // hpsi_m: H|psi_m>
    auto hpsi_m = ct::Tensor(ct::DataTypeToEnum<T>::value,
                             ct::DeviceTypeToEnum<ct_Device>::value,
                             {this->n_basis_});
    // spsi_m: S|psi_m>
    auto spsi_m = ct::Tensor(ct::DataTypeToEnum<T>::value,
                             ct::DeviceTypeToEnum<ct_Device>::value,
                             {this->n_basis_});
    // ppsi: preconditioned workspace / H|cg>
    auto ppsi = ct::Tensor(ct::DataTypeToEnum<T>::value,
                           ct::DeviceTypeToEnum<ct_Device>::value,
                           {this->n_basis_});
    // gradient
    auto grad = ct::Tensor(ct::DataTypeToEnum<T>::value,
                           ct::DeviceTypeToEnum<ct_Device>::value,
                           {this->n_basis_});
    // s_grad: S|grad>
    auto s_grad = ct::Tensor(ct::DataTypeToEnum<T>::value,
                             ct::DeviceTypeToEnum<ct_Device>::value,
                             {this->n_basis_});
    // cg: CG direction
    auto cg = ct::Tensor(ct::DataTypeToEnum<T>::value,
                         ct::DeviceTypeToEnum<ct_Device>::value,
                         {this->n_basis_});
    // scg: S|cg>
    auto scg = ct::Tensor(ct::DataTypeToEnum<T>::value,
                          ct::DeviceTypeToEnum<ct_Device>::value,
                          {this->n_basis_});
    // lagrange: for orthogonalization
    auto lagrange = ct::Tensor(ct::DataTypeToEnum<T>::value,
                               ct::DeviceTypeToEnum<ct_Device>::value,
                               {this->n_band_});

    // Set up preconditioner tensor
    ct::Tensor prec_tensor;
    if (prec != nullptr)
    {
        prec_tensor = ct::TensorMap(const_cast<Real*>(prec),
                                    ct::DataTypeToEnum<Real>::value,
                                    ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value,
                                    ct::TensorShape({dim}))
                          .template to_device<ct_Device>();
    }
    else
    {
        prec_tensor = ct::Tensor(ct::DataTypeToEnum<Real>::value,
                                 ct::DeviceTypeToEnum<ct_Device>::value,
                                 {this->n_basis_});
        prec_tensor.set_value(static_cast<Real>(1.0));
    }

    ModuleBase::Memory::record("DiagoPPCG", this->n_basis_ * 12 * sizeof(T));

    // Outer retry loop matching ABACUS CG pattern
    int ntry = 0;
    this->notconv_ = 0;
    int avg = 0;
    eigen.zero();
    auto eigen_pack = eigen.accessor<Real, 1>();

    // Slice of psi for subspace function
    auto psi_temp = psi.slice({0, 0}, {nband, dim});

    do
    {
        // Reset for this try
        this->notconv_ = 0;
        avg = 0;

        // Subspace diagonalization to improve initial guess (matching ABACUS CG)
        if (ntry > 0 && subspace_func_)
        {
            ct::TensorMap psi_map(psi.data<T>(),
                                  ct::DataTypeToEnum<T>::value,
                                  ct::DeviceTypeToEnum<ct_Device>::value,
                                  ct::TensorShape({nband, dim}));
            const bool assume_S_orthogonal = true;
            this->subspace_func_(psi_temp.data<T>(),
                                 psi_map.data<T>(),
                                 dim,
                                 nband,
                                 assume_S_orthogonal);
            psi_temp.sync(psi_map);
        }
        else if (ntry == 0 && need_subspace_ && subspace_func_)
        {
            ct::TensorMap psi_map(psi.data<T>(),
                                  ct::DataTypeToEnum<T>::value,
                                  ct::DeviceTypeToEnum<ct_Device>::value,
                                  ct::TensorShape({nband, dim}));
            const bool assume_S_orthogonal = false;
            this->subspace_func_(psi_temp.data<T>(),
                                 psi_map.data<T>(),
                                 dim,
                                 nband,
                                 assume_S_orthogonal);
            psi_temp.sync(psi_map);
        }

        // Main band-by-band loop
        for (int m = 0; m < this->n_band_; m++)
        {
            // Copy initial guess
            psi_m.sync(psi[m]);

            // Compute S|psi_m>
            this->spsi_func_(psi_m.data<T>(), spsi_m.data<T>(), this->n_basis_, 1);

            // Schmidt orthogonalize against lower bands
            this->schmidt_orth(m, psi, spsi_m, psi_m);

            // Recompute S|psi_m> after orthogonalization
            this->spsi_func_(psi_m.data<T>(), spsi_m.data<T>(), this->n_basis_, 1);

            // Compute H|psi_m>
            this->hpsi_func_(psi_m.data<T>(), hpsi_m.data<T>(), this->n_basis_, 1);

            // Initial eigenvalue estimate: λ = <ψ|H|ψ>
            eigen_pack[m] = ModuleBase::dot_real_op<T, Device>()(this->n_basis_, psi_m.data<T>(), hpsi_m.data<T>());

            int iter = 0;
            Real gg_last = 0.0;
            Real cg_norm = 0.0;
            Real theta = 0.0;
            bool converged = false;
            // g0 = P * S|grad> for CG update (persistent between iterations)
            auto g0_tensor = ct::Tensor(ct::DataTypeToEnum<T>::value,
                                        ct::DeviceTypeToEnum<ct_Device>::value,
                                        {this->n_basis_});

            do
            {
                // Compute preconditioned gradient (ABACUS calc_grad)
                this->compute_gradient(prec_tensor, hpsi_m, spsi_m, eigen_pack[m], grad, ppsi);

                // Orthogonalize gradient against lower bands (ABACUS orth_grad)
                this->orth_gradient(psi, m, grad, s_grad, lagrange);

                // Update CG direction with persistent g0 and Polak-Ribiere (ABACUS calc_gamma_cg)
                Real gamma = this->update_cg_direction(cg, scg, grad, s_grad, prec_tensor, g0_tensor, gg_last, iter);

                // ABACUS CG correction: cg -= gamma * cg_norm * sin(theta) * psi_m
                if (iter > 0 && std::abs(gamma) > std::numeric_limits<Real>::epsilon())
                {
                    const Real norma = gamma * cg_norm * std::sin(theta);
                    if (std::abs(norma) > std::numeric_limits<Real>::epsilon())
                    {
                        T znorma = static_cast<T>(-norma);
                        ModuleBase::axpy_op<T, Device>()(this->n_basis_, &znorma, psi_m.data<T>(), 1, cg.data<T>(), 1);
                        ModuleBase::axpy_op<T, Device>()(this->n_basis_, &znorma, spsi_m.data<T>(), 1, scg.data<T>(), 1);
                    }
                }

                // Compute H|cg> and S|cg>
                this->hpsi_func_(cg.data<T>(), ppsi.data<T>(), this->n_basis_, 1);
                this->spsi_func_(cg.data<T>(), scg.data<T>(), this->n_basis_, 1);

                // Line minimization in [ψ, cg] subspace (ABACUS update_psi)
                converged = this->line_minimization(ppsi,
                                                    cg,
                                                    scg,
                                                    ethr_band[m],
                                                    cg_norm,
                                                    theta,
                                                    eigen_pack[m],
                                                    psi_m,
                                                    spsi_m,
                                                    hpsi_m);

                iter++;
            } while (!converged && iter < this->pw_diag_nmax_);

            // Store converged eigenvector
            psi[m].sync(psi_m);

            if (!converged)
            {
                ++this->notconv_;
            }

            avg += iter;

            // Reorder eigenvalues if necessary
            if (m > 0)
            {
                ModuleBase::GlobalFunc::NOTE("reorder bands!");
                if (eigen_pack[m] - eigen_pack[m - 1] < -2.0 * this->pw_diag_thr_)
                {
                    int ii = 0;
                    for (ii = m - 2; ii >= 0; ii--)
                    {
                        if (eigen_pack[m] - eigen_pack[ii] > 2.0 * this->pw_diag_thr_)
                        {
                            break;
                        }
                    }
                    ii++;
                    Real e0 = eigen_pack[m];
                    ppsi.sync(psi[m]);

                    for (int jj = m; jj >= ii + 1; jj--)
                    {
                        eigen_pack[jj] = eigen_pack[jj - 1];
                        psi[jj].sync(psi[jj - 1]);
                    }
                    eigen_pack[ii] = e0;
                    psi[ii].sync(ppsi);
                }
            }
        }

        avg = (this->n_band_ > 0) ? avg / this->n_band_ : 0;
        avg_iter_ += avg;
        ntry++;

    } while (this->test_exit_cond(ntry, this->notconv_));

    ModuleBase::timer::end("DiagoPPCG", "diag");
    return this->avg_iter_;
}

// ============================================================
// Block PPCG: true simultaneous block algorithm
// (Vecharynski–Yang–Pask, 2015)
// ============================================================

template <typename T, typename Device>
double DiagoPPCG<T, Device>::diag_block(const HPsiFunc& hpsi_func,
                                        const SPsiFunc& spsi_func,
                                        const int ld_psi,
                                        const int nband,
                                        const int dim,
                                        T* psi_in,
                                        Real* eigenvalue_in,
                                        const std::vector<double>& ethr_band,
                                        const Real* prec)
{
    ModuleBase::TITLE("DiagoPPCG", "diag_block");
    ModuleBase::timer::start("DiagoPPCG", "diag_block");

    // ---- parameters ----
    const int k = nband;          // block size = number of bands
    const int two_k = 2 * k;      // subspace dimension for Rayleigh-Ritz
    const int ld = ld_psi;        // leading dimension

    this->n_band_ = k;
    this->n_basis_ = dim;
    this->ld_psi_ = ld;
    this->hpsi_func_ = hpsi_func;
    this->spsi_func_ = spsi_func;
    this->notconv_ = 0;
    this->avg_iter_ = 0.0;

    const int max_iter = (this->pw_diag_nmax_ > 0)
                             ? this->pw_diag_nmax_
                             : DiagoIterAssist<T, Device>::PW_DIAG_NMAX;
    const Real thr = this->pw_diag_thr_;

    // ---- tensor views for I/O ----
    auto psi_view = ct::TensorMap(psi_in,
                                  ct::DataTypeToEnum<T>::value,
                                  ct::DeviceTypeToEnum<ct_Device>::value,
                                  ct::TensorShape({k, ld}));
    auto eigen_view = ct::TensorMap(eigenvalue_in,
                                    ct::DataTypeToEnum<Real>::value,
                                    ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value,
                                    ct::TensorShape({k}));

    // ---- working memory helpers ----
    auto alloc = [&](const std::string& tag, int n) {
        auto t = ct::Tensor(ct::DataTypeToEnum<T>::value,
                            ct::DeviceTypeToEnum<ct_Device>::value,
                            {n});
        ModuleBase::Memory::record("DiagoPPCG::block", n * sizeof(T));
        return t;
    };

    // ---- allocate block arrays (each size = k * ld) ----
    ct::Tensor X_t  = alloc("X",  k * ld);
    ct::Tensor HX_t = alloc("HX", k * ld);
    ct::Tensor SX_t = alloc("SX", k * ld);
    ct::Tensor P_t  = alloc("P",  k * ld);
    ct::Tensor HP_t = alloc("HP", k * ld);
    ct::Tensor SP_t = alloc("SP", k * ld);
    ct::Tensor Z_t  = alloc("Z",  k * ld);
    ct::Tensor SZ_t = alloc("SZ", k * ld);
    ct::Tensor Zo_t = alloc("Zo", k * ld);   // Z_old
    ct::Tensor SZo_t= alloc("SZo",k * ld);   // SZ_old

    T* const X   = X_t.data<T>();
    T* const HX  = HX_t.data<T>();
    T* const SX  = SX_t.data<T>();
    T* const P   = P_t.data<T>();
    T* const HP  = HP_t.data<T>();
    T* const SP  = SP_t.data<T>();
    T* const Z   = Z_t.data<T>();
    T* const SZ  = SZ_t.data<T>();
    T* const Zo_ = Zo_t.data<T>();
    T* const SZo = SZo_t.data<T>();

    // ---- misc working arrays ----
    ct::Tensor lag_t = alloc("lag", k * k);   // Lagrange multipliers
    T* const lag = lag_t.data<T>();

    // Subspace matrices on CPU (for LAPACK)
    ct::Tensor hcc_t(ct::DataTypeToEnum<T>::value,
                     ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value,
                     {two_k * two_k});
    ct::Tensor scc_t(ct::DataTypeToEnum<T>::value,
                     ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value,
                     {two_k * two_k});
    ct::Tensor vcc_t(ct::DataTypeToEnum<T>::value,
                     ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value,
                     {two_k * two_k});
    ct::Tensor eig_t(ct::DataTypeToEnum<Real>::value,
                     ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value,
                     {two_k});
    T* const hcc = hcc_t.data<T>();
    T* const scc = scc_t.data<T>();
    T* const vcc = vcc_t.data<T>();
    Real* const eig_all = eig_t.data<Real>();

    // ---- preconditioner ----
    ct::Tensor prec_t;
    const Real* prec_ptr = nullptr;
    if (prec != nullptr) {
        prec_t = ct::TensorMap(const_cast<Real*>(prec),
                               ct::DataTypeToEnum<Real>::value,
                               ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value,
                               ct::TensorShape({dim}))
                     .template to_device<ct_Device>();
        prec_ptr = prec_t.data<Real>();
    } else {
        prec_t = ct::Tensor(ct::DataTypeToEnum<Real>::value,
                            ct::DeviceTypeToEnum<ct_Device>::value,
                            {dim});
        prec_t.set_value(static_cast<Real>(1.0));
        prec_ptr = prec_t.data<Real>();
    }

    // Own constants
    const T cone  = static_cast<T>(1.0);
    const T czero = static_cast<T>(0.0);
    const T cneg  = static_cast<T>(-1.0);
    const T* const pcone  = &cone;
    const T* const pczero = &czero;
    const T* const pcneg  = &cneg;

    // =========================================================
    // STEP 1:  Copy initial guess and S-orthonormalize
    // =========================================================
    base_device::memory::synchronize_memory_op<T, Device, Device>()(
        X, psi_view.data<T>(), k * ld);

    if (subspace_func_) {
        subspace_func_(X, psi_view.data<T>(), dim, k, false);
        base_device::memory::synchronize_memory_op<T, Device, Device>()(
            X, psi_view.data<T>(), k * ld);
    }

    // Compute HX = H*X, SX = S*X
    hpsi_func(X, HX, ld, k);
    spsi_func(X, SX, ld, k);

    // Initial eigenvalues: λ_i = Re( X[:,i]^H * HX[:,i] )
    // Since X is S-orthonormal, X^H * SX = I
    for (int i = 0; i < k; ++i) {
        Real re = 0.0;
        for (int j = 0; j < dim; ++j) {
            re += (std::conj(X[j + i * ld]) * HX[j + i * ld]).real();
        }
        eigenvalue_in[i] = re;
    }

    // =========================================================
    // STEP 2:  Compute initial preconditioned residual Z
    // =========================================================
    // Z = HX - X * diag(λ)  (residual per band)
    for (int i = 0; i < k; ++i) {
        const Real lam = eigenvalue_in[i];
        for (int j = 0; j < dim; ++j) {
            Z[j + i * ld] = HX[j + i * ld] - static_cast<T>(lam) * X[j + i * ld];
        }
    }

    // Precondition: Z[i] = Z[i] / prec[:]  (element-wise)
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < dim; ++j) {
            Z[j + i * ld] = Z[j + i * ld] / static_cast<T>(prec_ptr[j]);
        }
    }

    // Orthogonalize Z against X:  lag = X^H * S * Z  →  Z -= X * lag
    // First compute SZ = S*Z
    spsi_func(Z, SZ, ld, k);

    // lag[i][j] = X[:,i]^H * SZ[:,j] = (X^H * SZ)[i][j]
    ModuleBase::gemm_op<T, Device>()('C', 'N', k, k, dim,
                                     pcone, X, ld, SZ, ld,
                                     pczero, lag, k);
    // Z -= X * lag
    ModuleBase::gemm_op<T, Device>()('N', 'N', dim, k, k,
                                     pcneg, X, ld, lag, k,
                                     pcone, Z, ld);
    // Recompute SZ for the orthogonalized Z
    spsi_func(Z, SZ, ld, k);

    // Copy to Z_old, SZ_old
    base_device::memory::synchronize_memory_op<T, Device, Device>()(Zo_, Z, k * ld);
    base_device::memory::synchronize_memory_op<T, Device, Device>()(SZo, SZ, k * ld);

    // Search direction P = Z
    base_device::memory::synchronize_memory_op<T, Device, Device>()(P, Z, k * ld);
    // HP = H*P, SP = S*P
    hpsi_func(P, HP, ld, k);
    spsi_func(P, SP, ld, k);

    // =========================================================
    // STEP 3:  Main PPCG iteration
    // =========================================================
    Real old_eigenvalues[64];   // k <= 64 typical for PW
    for (int i = 0; i < k; ++i) old_eigenvalues[i] = eigenvalue_in[i];

    int iter = 0;
    for (; iter < max_iter; ++iter)
    {
        // ---- 3a. Form 2k × 2k subspace matrices hcc, scc ----
        // hcc = [X, P]^H * H * [X, P] = [X^H*HX , X^H*HP;
        //                                P^H*HX , P^H*HP]
        // scc = [X, P]^H * S * [X, P] = [X^H*SX , X^H*SP;
        //                                P^H*SX , P^H*SP]

        // Work matrix tmp_k2[k * k] for intermediate gemm results
        T tmp_k2[4096]; // enough for k <= 64 → k*k <= 4096

        // Block (0,0): X^H * HX
        ModuleBase::gemm_op<T, Device>()('C', 'N', k, k, dim,
                                         pcone, X, ld, HX, ld,
                                         pczero, tmp_k2, k);
        // Copy to hcc[0:k, 0:k] (column-major, ldh = two_k)
        for (int col = 0; col < k; ++col)
            for (int row = 0; row < k; ++row)
                hcc[row + col * two_k] = tmp_k2[row + col * k];

        // Block (0,1): X^H * HP  → hcc[0:k, k:2k]
        ModuleBase::gemm_op<T, Device>()('C', 'N', k, k, dim,
                                         pcone, X, ld, HP, ld,
                                         pczero, tmp_k2, k);
        for (int col = 0; col < k; ++col)
            for (int row = 0; row < k; ++row)
                hcc[row + (col + k) * two_k] = tmp_k2[row + col * k];

        // Block (1,1): P^H * HP  → hcc[k:2k, k:2k]
        ModuleBase::gemm_op<T, Device>()('C', 'N', k, k, dim,
                                         pcone, P, ld, HP, ld,
                                         pczero, tmp_k2, k);
        for (int col = 0; col < k; ++col)
            for (int row = 0; row < k; ++row)
                hcc[(row + k) + (col + k) * two_k] = tmp_k2[row + col * k];

        // Block (1,0) = Hermitian of (0,1)
        for (int col = 0; col < k; ++col)
            for (int row = 0; row < k; ++row)
                hcc[(row + k) + col * two_k] = std::conj(hcc[row + (col + k) * two_k]);

        // ---- Same for scc ----
        // Block (0,0): X^H * SX
        ModuleBase::gemm_op<T, Device>()('C', 'N', k, k, dim,
                                         pcone, X, ld, SX, ld,
                                         pczero, tmp_k2, k);
        for (int col = 0; col < k; ++col)
            for (int row = 0; row < k; ++row)
                scc[row + col * two_k] = tmp_k2[row + col * k];

        // Block (0,1): X^H * SP
        ModuleBase::gemm_op<T, Device>()('C', 'N', k, k, dim,
                                         pcone, X, ld, SP, ld,
                                         pczero, tmp_k2, k);
        for (int col = 0; col < k; ++col)
            for (int row = 0; row < k; ++row)
                scc[row + (col + k) * two_k] = tmp_k2[row + col * k];

        // Block (1,1): P^H * SP
        ModuleBase::gemm_op<T, Device>()('C', 'N', k, k, dim,
                                         pcone, P, ld, SP, ld,
                                         pczero, tmp_k2, k);
        for (int col = 0; col < k; ++col)
            for (int row = 0; row < k; ++row)
                scc[(row + k) + (col + k) * two_k] = tmp_k2[row + col * k];

        // Block (1,0) = Hermitian of (0,1)
        for (int col = 0; col < k; ++col)
            for (int row = 0; row < k; ++row)
                scc[(row + k) + col * two_k] = std::conj(scc[row + (col + k) * two_k]);

        // ---- 3b. Solve generalized EVP: hcc * Y = scc * Y * diag(Λ) ----
        // Copy hcc → vcc (diag_hegvd overwrites hcc, stores eigenvectors in vcc)
        base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, base_device::DEVICE_CPU>()(
            vcc, hcc, two_k * two_k);

        DiagoIterAssist<T, Device>::diag_hegvd(two_k, k, hcc, scc, two_k, eig_all, vcc);

        // eig_all[0..k-1] are the k smallest eigenvalues (sorted ascending)
        // vcc[:, 0..k-1] are the corresponding eigenvectors (2k × k)

        // ---- 3c. Rotate X, HX, SX ----
        // X_new  = X  * V1 + P  * V2   where V1=vcc[0:k,0:k], V2=vcc[k:2k,0:k]
        // HX_new = HX * V1 + HP * V2
        // SX_new = SX * V1 + SP * V2

        // Temporary storage: temp[dim * k]
        ct::Tensor temp_t = alloc("temp", dim * k);
        T* const temp = temp_t.data<T>();

        // X_new = X * V1  (gemm: dim×k = dim×k * k×k)
        ModuleBase::gemm_op<T, Device>()('N', 'N', dim, k, k,
                                         pcone, X, ld, vcc, two_k,
                                         pczero, temp, dim);
        // X_new += P * V2
        ModuleBase::gemm_op<T, Device>()('N', 'N', dim, k, k,
                                         pcone, P, ld, &vcc[k], two_k,
                                         pcone, temp, dim);
        base_device::memory::synchronize_memory_op<T, Device, Device>()(X, temp, dim * k);

        // HX_new = HX * V1 + HP * V2
        ModuleBase::gemm_op<T, Device>()('N', 'N', dim, k, k,
                                         pcone, HX, ld, vcc, two_k,
                                         pczero, temp, dim);
        ModuleBase::gemm_op<T, Device>()('N', 'N', dim, k, k,
                                         pcone, HP, ld, &vcc[k], two_k,
                                         pcone, temp, dim);
        base_device::memory::synchronize_memory_op<T, Device, Device>()(HX, temp, dim * k);

        // SX_new = SX * V1 + SP * V2
        ModuleBase::gemm_op<T, Device>()('N', 'N', dim, k, k,
                                         pcone, SX, ld, vcc, two_k,
                                         pczero, temp, dim);
        ModuleBase::gemm_op<T, Device>()('N', 'N', dim, k, k,
                                         pcone, SP, ld, &vcc[k], two_k,
                                         pcone, temp, dim);
        base_device::memory::synchronize_memory_op<T, Device, Device>()(SX, temp, dim * k);

        // ---- 3d. Convergence check ----
        int nconv = 0;
        for (int i = 0; i < k; ++i) {
            const Real de = std::abs(eig_all[i] - old_eigenvalues[i]);
            old_eigenvalues[i] = eig_all[i];
            eigenvalue_in[i] = eig_all[i];
            if (de < static_cast<Real>(ethr_band[i])) ++nconv;
        }
        this->avg_iter_ += 1.0;
        if (nconv == k) break;  // all bands converged

        // ---- 3e. Compute new preconditioned residual ----
        // Z_new = HX - X * diag(eigenvalues)
        for (int i = 0; i < k; ++i) {
            const Real lam = eig_all[i];
            for (int j = 0; j < dim; ++j) {
                Z[j + i * ld] = HX[j + i * ld] - static_cast<T>(lam) * X[j + i * ld];
            }
        }

        // Precondition (element-wise division)
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < dim; ++j) {
                Z[j + i * ld] = Z[j + i * ld] / static_cast<T>(prec_ptr[j]);
            }
        }

        // Orthogonalize Z against X: lag = X^H * (S*Z)
        spsi_func(Z, SZ, ld, k);
        ModuleBase::gemm_op<T, Device>()('C', 'N', k, k, dim,
                                         pcone, X, ld, SZ, ld,
                                         pczero, lag, k);
        ModuleBase::gemm_op<T, Device>()('N', 'N', dim, k, k,
                                         pcneg, X, ld, lag, k,
                                         pcone, Z, ld);
        spsi_func(Z, SZ, ld, k);

        // ---- 3f. Polak–Ribiere beta (block, single scalar) ----
        // beta = max(0, trace(Z_new^H * SZ_new - Z_new^H * SZ_old)
        //               / trace(Z_old^H * SZ_old) )
        Real tr_num = 0.0, tr_den = 0.0;
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < dim; ++j) {
                tr_num += (std::conj(Z[j + i * ld]) * SZ[j + i * ld]
                          - std::conj(Zo_[j + i * ld]) * SZo[j + i * ld]).real();
                tr_den += (std::conj(Zo_[j + i * ld]) * SZo[j + i * ld]).real();
            }
        }
        const Real beta = (tr_den > std::numeric_limits<Real>::epsilon())
                              ? std::max(Real(0.0), tr_num / tr_den) : Real(0.0);

        // ---- 3g. Update search direction ----
        // P_new = Z_new + beta * P
        // HP_new = H*P_new, SP_new = S*P_new
        // But first: compute HP_new = H*Z_new + beta * HP
        // (since HP_new = H*(Z_new + beta*P) = H*Z_new + beta*HP)

        // Save old Z → Z_old, SZ → SZ_old
        base_device::memory::synchronize_memory_op<T, Device, Device>()(Zo_, Z, k * ld);
        base_device::memory::synchronize_memory_op<T, Device, Device>()(SZo, SZ, k * ld);

        // Store H*Z_new in HP_temp, then combine
        ct::Tensor HZ_t = alloc("HZ", dim * k);
        hpsi_func(Z, HZ_t.data<T>(), ld, k);

        // P_new = Z_new + beta * P
        // Z currently holds Z_new
        for (int i = 0; i < k * ld; ++i) {
            P[i] = Z[i] + static_cast<T>(beta) * P[i];
        }

        // HP_new = H*Z_new + beta * HP
        for (int i = 0; i < k * ld; ++i) {
            HP[i] = HZ_t.data<T>()[i] + static_cast<T>(beta) * HP[i];
        }
        // SP_new = S*Z_new + beta * SP
        for (int i = 0; i < k * ld; ++i) {
            SP[i] = SZ[i] + static_cast<T>(beta) * SP[i];
        }
    }

    // ---- convergence summary ----
    if (iter >= max_iter) {
        this->notconv_ = k;
    }
    this->avg_iter_ = (k > 0) ? this->avg_iter_ / k : 0.0;

    // ---- write back ----
    base_device::memory::synchronize_memory_op<T, Device, Device>()(
        psi_view.data<T>(), X, k * ld);

    ModuleBase::timer::end("DiagoPPCG", "diag_block");
    return this->avg_iter_;
}

// Template instantiations
template class DiagoPPCG<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoPPCG<std::complex<float>, base_device::DEVICE_GPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hsolver
