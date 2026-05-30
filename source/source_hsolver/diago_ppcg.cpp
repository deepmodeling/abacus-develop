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
    for (int i = 0; i < this->n_basis_; i++)
    {
        g0.data<T>()[i] = s_grad.data<T>()[i] * prec.data<Real>()[i];
    }

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
        std::cout << " DiagoPPGC: psi_norm <= 0.0, m = " << m << std::endl;
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

    this->hpsi_func_ = hpsi_func;
    this->spsi_func_ = spsi_func;

    // Get convergence parameters from DiagoIterAssist
    this->pw_diag_nmax_ = DiagoIterAssist<T, Device>::PW_DIAG_NMAX;
    this->pw_diag_thr_ = DiagoIterAssist<T, Device>::PW_DIAG_THR;

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
    return avg;
}

// Template instantiations
template class DiagoPPCG<std::complex<float>, base_device::DEVICE_CPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoPPCG<std::complex<float>, base_device::DEVICE_GPU>;
template class DiagoPPCG<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hsolver
