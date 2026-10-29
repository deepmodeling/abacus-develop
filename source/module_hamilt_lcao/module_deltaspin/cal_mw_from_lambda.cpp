#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_parameter/parameter.h"
#include "spin_constrain.h"
#include "module_hamilt_pw/hamilt_pwdft/onsite_projector.h"
#include "module_hamilt_pw/hamilt_pwdft/operator_pw/nonlocal_pw.h"
#include "module_base/parallel_reduce.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_hsolver/hsolver_pw.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/onsite_op.h"
#include <algorithm>
#if defined(__CUDA) || defined(__ROCM)
#include <ATen/core/tensor.h>
#endif
#if defined(__CUDA)
#include <base/macros/macros.h>
#include <cuda_runtime_api.h>
#endif

#ifdef __LCAO
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/dspin_lcao.h"
#endif

namespace
{

#if defined(__CUDA)
using Complex = std::complex<double>;
using GPU = base_device::DEVICE_GPU;
using NonlocalPW = hamilt::Nonlocal<hamilt::OperatorPW<Complex, GPU>>;

NonlocalPW* find_nonlocal_operator(hamilt::Hamilt<Complex, GPU>* hamiltonian)
{
    for (auto* op = hamiltonian->ops; op != nullptr; op = op->next_op)
    {
        if (auto* nonlocal = dynamic_cast<NonlocalPW*>(op)) { return nonlocal; }
    }
    return nullptr;
}

hamilt::Operator<Complex, GPU>::Action projected_nonlocal_action(
    NonlocalPW* nonlocal, const std::vector<Complex>& projections, int bands, int polarizations)
{
    const size_t count = size_t(nonlocal->get_ppcell()->nkb) * bands * polarizations;
    if (projections.size() != count) { return {}; }
    auto* projected = nonlocal->prepare_becp(bands * polarizations);
    base_device::memory::synchronize_memory_op<Complex, GPU, base_device::DEVICE_CPU>()(
        nullptr, nullptr, projected, projections.data(), count);
    return [nonlocal, projected](const hamilt::Operator<Complex, GPU>* op,
        int columns, int leading, int npol, const Complex*, Complex* output, int plane_waves) {
        if (op != nonlocal) { return false; }
        nonlocal->act_projected(columns, leading, npol, projected, output, plane_waves);
        return true;
    };
}

class NonlocalBasisObserver
{
    using Solver = hsolver::HSolverPW<Complex, GPU>;
    using Cache = std::vector<std::vector<Complex>>;

    NonlocalPW* nonlocal_;
    const int bands_;
    const int rows_;
    Cache& cache_;
    ct::Tensor basis_;
    ct::Tensor rotated_;
    Solver::BasisObserver callbacks_;

    void capture(int ik, int first, int count, bool cached_action)
    {
        auto* target = basis_.data<Complex>() + size_t(first) * rows_;
        const size_t size = size_t(count) * rows_;
        if (cached_action)
        {
            assert(first == 0 && count == bands_ && cache_[ik].size() == size);
            base_device::memory::synchronize_memory_op<Complex, GPU, base_device::DEVICE_CPU>()(
                nullptr, nullptr, target, cache_[ik].data(), size);
        }
        else
        {
            base_device::memory::synchronize_memory_op<Complex, GPU, GPU>()(
                nullptr, nullptr, target, nonlocal_->get_becp(), size);
        }
    }

    void rotate(int ik, int count, const Complex* coefficients, int leading, bool final)
    {
        auto* basis = basis_.data<Complex>();
        auto* rotated = rotated_.data<Complex>();
        hsolver::gemm_op<Complex, GPU>()(nullptr, 'N', 'N', rows_, bands_, count,
            &ModuleBase::ONE, basis, rows_, coefficients, leading, &ModuleBase::ZERO, rotated, rows_);
        const size_t size = size_t(rows_) * bands_;
        if (final)
        {
            cache_[ik].resize(size);
            base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_CPU, GPU>()(
                nullptr, nullptr, cache_[ik].data(), rotated, size);
        }
        else
        {
            base_device::memory::synchronize_memory_op<Complex, GPU, GPU>()(
                nullptr, nullptr, basis, rotated, size);
        }
    }

  public:
    NonlocalBasisObserver(NonlocalPW* nonlocal, int bands, int polarizations, int capacity, Cache& cache)
        : nonlocal_(nonlocal), bands_(bands), rows_(nonlocal->get_ppcell()->nkb * polarizations), cache_(cache),
          basis_(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice, {capacity, rows_}),
          rotated_(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice, {bands_, rows_})
    {
        callbacks_.capture = [this](int ik, int first, int count, bool cached) { capture(ik, first, count, cached); };
        callbacks_.rotate = [this](int ik, int count, const Complex* coefficients, int leading, bool final) {
            rotate(ik, count, coefficients, leading, final);
        };
        callbacks_.invalidate = [this](int ik) { cache_[ik].clear(); };
    }

    const Solver::BasisObserver* callbacks() const { return &callbacks_; }
};

std::vector<double> pack_constraint_field(int atoms, int nkb, const int* projectors,
                                         const ModuleBase::Vector3<double>* field)
{
    std::vector<double> coefficients(3 * nkb);
    int offset = 0;
    for (int atom = 0; atom < atoms; ++atom)
    {
        for (int projector = 0; projector < projectors[atom]; ++projector)
        {
            for (int axis = 0; axis < 3; ++axis)
            {
                coefficients[3 * (offset + projector) + axis] = field[atom][axis];
            }
        }
        offset += projectors[atom];
    }
    return coefficients;
}
#endif

void apply_spin_constraint_cpu(int atoms, int bands, int npol, int nkb, int sign,
                               const int* projectors, const ModuleBase::Vector3<double>* lambda,
                               const std::complex<double>* becp, std::complex<double>* ps)
{
    int offset = 0;
    for (int atom = 0; atom < atoms; ++atom)
    {
        const std::complex<double> z(lambda[atom].z, 0.0);
        const std::complex<double> plus(lambda[atom].x, lambda[atom].y);
        const std::complex<double> minus(lambda[atom].x, -lambda[atom].y);
        for (int band = 0; band < bands; ++band)
        {
            for (int ip = 0; ip < projectors[atom]; ++ip)
            {
                const int index = band * npol * nkb + offset + ip;
                if (npol == 2)
                {
                    ps[index] = z * becp[index] + minus * becp[index + nkb];
                    ps[index + nkb] = plus * becp[index] - z * becp[index + nkb];
                }
                else { ps[index] = lambda[atom].z * sign * becp[index]; }
            }
        }
        offset += projectors[atom];
    }
}

} // namespace

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::calculate_delta_hcc(
    std::complex<double>* h_tmp, const std::complex<double>* becp_k,
    const ModuleBase::Vector3<double>* delta_lambda, const int nbands,
    const int nkb, const int* nh_iat, const int sign)
{
    using Complex = std::complex<double>;
    const int rows = nkb * this->npol_;
    const int size = rows * nbands;
    const int atoms = this->Mi_.size();
#if defined(__CUDA) || defined(__ROCM)
    if (PARAM.inp.device == "gpu")
    {
        const base_device::DEVICE_GPU* gpu = nullptr;
        const base_device::DEVICE_CPU* cpu = nullptr;
        ct::Tensor ps(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice, {nbands, rows});
#if defined(__CUDA)
        // Each band shares the same three constraint coefficients per projector.
        // Only these coefficients cross the host/device boundary.
        const auto coefficients = pack_constraint_field(atoms, nkb, nh_iat, delta_lambda);
        ct::Tensor lambda(ct::DataType::DT_DOUBLE, ct::DeviceType::GpuDevice, {nkb, 3});
        base_device::memory::synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            gpu, cpu, lambda.data<double>(), coefficients.data(), coefficients.size());
        hamilt::spin_constraint_projectors_gpu(nbands, this->npol_, nkb, sign, lambda.data<double>(),
                                               becp_k, ps.data<Complex>());
#else
        std::vector<Complex> input(size), output(size);
        base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(
            cpu, gpu, input.data(), becp_k, size);
        apply_spin_constraint_cpu(atoms, nbands, this->npol_, nkb, sign, nh_iat, delta_lambda,
                                  input.data(), output.data());
        base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
            gpu, cpu, ps.data<Complex>(), output.data(), size);
#endif
        hsolver::gemm_op<Complex, base_device::DEVICE_GPU>()(
            gpu, 'C', 'N', nbands, nbands, rows, &ModuleBase::ONE, becp_k, rows,
            ps.data<Complex>(), rows, &ModuleBase::ONE, h_tmp, nbands);
        return;
    }
#endif
    std::vector<Complex> ps(size);
    apply_spin_constraint_cpu(atoms, nbands, this->npol_, nkb, sign, nh_iat, delta_lambda, becp_k, ps.data());
    hsolver::gemm_op<Complex, base_device::DEVICE_CPU>()(
        nullptr, 'C', 'N', nbands, nbands, rows, &ModuleBase::ONE, becp_k, rows,
        ps.data(), rows, &ModuleBase::ONE, h_tmp, nbands);
}

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::cal_mw_from_lambda(int i_step, const ModuleBase::Vector3<double>* delta_lambda)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "cal_mw_from_lambda");
    ModuleBase::timer::tick("spinconstrain::SpinConstrain", "cal_mw_from_lambda");
    // lambda has been updated in the lambda loop
#ifdef __LCAO
    if (PARAM.inp.basis_type == "lcao")
    {
        psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
        hamilt::Hamilt<std::complex<double>>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt);
        hsolver::HSolverLCAO<std::complex<double>> hsolver_t(this->ParaV, PARAM.inp.ks_solver);
        if (PARAM.inp.nspin == 2)
        {
            dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator)
                ->update_lambda();
        }
        else if (PARAM.inp.nspin == 4)
        {
            dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>*>(
                this->p_operator)
                ->update_lambda();
        }
        // diagonalization without update charge
        hsolver_t.solve(hamilt_t, psi_t[0], this->pelec, true);
        this->pelec->calculate_weights();
        this->pelec->calEBand();
        elecstate::ElecStateLCAO<std::complex<double>>* pelec_lcao
            = dynamic_cast<elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec);
        elecstate::cal_dm_psi(this->ParaV, pelec_lcao->wg, *psi_t, *(pelec_lcao->get_DM()));
        pelec_lcao->get_DM()->cal_DMR();
        this->cal_mi_lcao(i_step);
    }
    else
#endif
    {
        /*if (i_step == -1 && this->higher_mag_prec)
        {
            // std::cout<<__FILE__<<__LINE__<<"istep == 0"<<std::endl;
            if (PARAM.inp.device == "cpu")
            {
                psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
                hamilt::Hamilt<std::complex<double>>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>>*>(this->p_hamilt);
                hsolver::HSolver<std::complex<double>, base_device::DEVICE_CPU>* hsolver_t = static_cast<hsolver::HSolver<std::complex<double>, base_device::DEVICE_CPU>*>(this->phsol);
                hsolver_t->solve(hamilt_t, psi_t[0], this->pelec, PARAM.inp.ks_solver, true);
            }
            else
            {
                psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*>(this->psi);
                hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>*>(this->p_hamilt);
                hsolver::HSolver<std::complex<double>, base_device::DEVICE_GPU>* hsolver_t = static_cast<hsolver::HSolver<std::complex<double>, base_device::DEVICE_GPU>*>(this->phsol);
                hsolver_t->solve(hamilt_t, psi_t[0], this->pelec, PARAM.inp.ks_solver, true);
            }
            this->pelec->calculate_weights();
            this->cal_Mi_pw();
        }
        else*/
        {
            this->zero_Mi();
            int size_becp = 0;
            std::vector<std::complex<double>> becp_tmp;
            int nk = 0;
            int nkb = 0;
            int nbands = 0;
            int npol = 0;
            const int* nh_iat = nullptr;
            if (PARAM.inp.device == "cpu")
            {
                psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
                hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>*>(this->p_hamilt);
                auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::get_instance();
                nbands = psi_t->get_nbands();
                npol = psi_t->npol;
                nkb = onsite_p->get_tot_nproj();
                nk = psi_t->get_nk();
                nh_iat = &onsite_p->get_nh(0);
                size_becp = nbands * nkb * npol;
                becp_tmp.resize(size_becp * nk);
                std::vector<std::complex<double>> h_tmp(nbands * nbands), s_tmp(nbands * nbands);
                int initial_hs = 0;
                if(this->sub_h_save == nullptr)
                {
                    initial_hs = 1;
                    this->sub_h_save = new std::complex<double>[nbands * nbands * nk];
                    this->sub_s_save = new std::complex<double>[nbands * nbands * nk];
                    this->becp_save = new std::complex<double>[size_becp * nk];
                }
                for (int ik = 0; ik < nk; ++ik)
                {

                    psi_t->fix_k(ik);

                    std::complex<double>* h_k = this->sub_h_save + ik * nbands * nbands;
                    std::complex<double>* s_k = this->sub_s_save + ik * nbands * nbands;
                    std::complex<double>* becp_k = this->becp_save + ik * size_becp;
                    if(initial_hs)
                    {
                        /// update H(k) for each k point
                        hamilt_t->updateHk(ik);
                        hsolver::DiagoIterAssist<std::complex<double>>::cal_hs_subspace(hamilt_t, psi_t[0], h_k, s_k);
                        memcpy(becp_k, onsite_p->get_becp(), sizeof(std::complex<double>) * size_becp);
                    }
                    memcpy(h_tmp.data(), h_k, sizeof(std::complex<double>) * nbands * nbands);
                    memcpy(s_tmp.data(), s_k, sizeof(std::complex<double>) * nbands * nbands);
                    // update h_tmp by delta_lambda
                    const int sign = this->pelec->klist->isk[ik] == 0 ? 1 : -1;
                    if (i_step != -1) this->calculate_delta_hcc(h_tmp.data(), becp_k, delta_lambda, nbands, nkb, nh_iat, sign);

                    hsolver::DiagoIterAssist<std::complex<double>>::diag_responce(h_tmp.data(),
                                                                                  s_tmp.data(),
                                                                                  nbands,
                                                                                  becp_k,
                                                                                  &becp_tmp[ik * size_becp],
                                                                                  nkb * 2,
                                                                                  &this->pelec->ekb(ik, 0));
                }
            }
#if ((defined __CUDA) || (defined __ROCM))
            else
            {
                base_device::DEVICE_GPU* ctx = {};
                base_device::DEVICE_CPU* cpu_ctx = {};
                psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*>(this->psi);
                hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>*>(this->p_hamilt);
                auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();
                nbands = psi_t->get_nbands();
                npol = psi_t->npol;
                nkb = onsite_p->get_tot_nproj();
                nk = psi_t->get_nk();
                nh_iat = &onsite_p->get_nh(0);
                size_becp = nbands * nkb * npol;
                becp_tmp.resize(size_becp * nk);
                std::complex<double>* h_tmp = nullptr;
                std::complex<double>* s_tmp = nullptr;
                base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, h_tmp, nbands * nbands);
                base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, s_tmp, nbands * nbands);
                int initial_hs = 0;
                if(this->sub_h_save == nullptr)
                {
                    initial_hs = 1;
                    this->subspace_lambda = this->lambda_;
                    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, this->sub_h_save, nbands * nbands * nk);
                    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, this->sub_s_save, nbands * nbands * nk);
                    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, this->becp_save, size_becp * nk);
                }
                std::complex<double>* hpsi_host = nullptr;
#if defined(__CUDA)
                if (initial_hs && PARAM.inp.ks_solver == "dav_subspace" && !PARAM.inp.use_paw && !PARAM.globalv.use_uspp)
                {
                    const size_t count = size_t(nbands) * psi_t->get_nbasis() * nk;
                    if (count > this->hpsi_cache_capacity)
                    {
                        this->hpsi_cache.reset();
                        this->hpsi_cache_capacity = 0;
                        std::complex<double>* storage = nullptr;
                        const auto status = cudaHostAlloc(reinterpret_cast<void**>(&storage),
                                                          count * sizeof(*storage), cudaHostAllocDefault);
                        if (status == cudaSuccess)
                        {
                            this->hpsi_cache = {storage, [](std::complex<double>* pointer) {
                                cudaErrcheck(cudaFreeHost(pointer));
                            }};
                            this->hpsi_cache_capacity = count;
                        }
                        else
                        {
                            // Host caching is optional when the pinned-memory allocation cannot be satisfied.
                            GlobalV::ofs_running << "Hamiltonian action cache unavailable: "
                                                 << cudaGetErrorString(status) << "; recomputing H psi.\n";
                            cudaGetLastError();
                        }
                    }
                    hpsi_host = this->hpsi_cache.get();
                }
#endif
#if defined(__CUDA)
                NonlocalPW* nonlocal = nullptr;
                if (initial_hs && PARAM.inp.ks_solver == "dav_subspace" && !PARAM.inp.use_paw
                    && !PARAM.globalv.use_uspp && GlobalV::NPROC_IN_POOL == 1)
                {
                    nonlocal = find_nonlocal_operator(hamilt_t);
                    if (nonlocal && nonlocal->get_ppcell()->nkb > 0)
                    {
                        this->nonlocal_projections.resize(nk);
                    }
                    else
                    {
                        nonlocal = nullptr;
                    }
                }
#endif
                std::complex<double>* becp_pointer = nullptr;
                // allocate memory for becp_pointer in GPU device
                base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, becp_pointer, size_becp);
                for (int ik = 0; ik < nk; ++ik)
                {
                    psi_t->fix_k(ik);

                    std::complex<double>* h_k = this->sub_h_save + ik * nbands * nbands;
                    std::complex<double>* s_k = this->sub_s_save + ik * nbands * nbands;
                    std::complex<double>* becp_k = this->becp_save + ik * size_becp;
                    if(initial_hs)
                    {
                        /// update H(k) for each k point
                        hamilt_t->updateHk(ik);
                        typename hamilt::Operator<std::complex<double>, base_device::DEVICE_GPU>::Action projected_action;
#if defined(__CUDA)
                        size_t projection_count = 0;
                        if (nonlocal)
                        {
                            projection_count = size_t(nonlocal->get_ppcell()->nkb) * nbands * npol;
                            projected_action = projected_nonlocal_action(nonlocal, this->nonlocal_projections[ik], nbands, npol);
                        }
#endif
                        const size_t offset = size_t(ik) * nbands * psi_t->get_nbasis();
                        // Host action caching is restricted to norm-conserving inputs with S = I.
                        hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::cal_hs_subspace(
                            hamilt_t, psi_t[0], h_k, s_k, hpsi_host == nullptr ? nullptr : hpsi_host + offset,
                            hpsi_host != nullptr, projected_action);
#if defined(__CUDA)
                        if (nonlocal && !projected_action)
                        {
                            this->nonlocal_projections[ik].resize(projection_count);
                            base_device::memory::synchronize_memory_op<Complex, base_device::DEVICE_CPU, GPU>()(
                                cpu_ctx, ctx, this->nonlocal_projections[ik].data(), nonlocal->get_becp(), projection_count);
                        }
#endif
                        base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(ctx, ctx, becp_k, onsite_p->get_becp(), size_becp);
                    }
                    base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(ctx, ctx, h_tmp, h_k, nbands * nbands);
                    base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(ctx, ctx, s_tmp, s_k, nbands * nbands);
                    // update h_tmp by delta_lambda
                    const int sign = this->pelec->klist->isk[ik] == 0 ? 1 : -1;
                    if (i_step != -1) this->calculate_delta_hcc(h_tmp, becp_k, delta_lambda, nbands, nkb, nh_iat, sign);

                    hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::diag_responce(h_tmp,
                                                                                  s_tmp,
                                                                                  nbands,
                                                                                  becp_k,
                                                                                  becp_pointer,
                                                                                  nkb * npol,
                                                                                  &this->pelec->ekb(ik, 0));
                    // copy becp_pointer from GPU to CPU
                    base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(cpu_ctx, ctx, &becp_tmp[ik * size_becp], becp_pointer, size_becp);   
                }

                // free memory for becp_pointer in GPU device
                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, becp_pointer);
                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, h_tmp);
                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, s_tmp);
            }
#endif
            // calculate weights from ekb to update wg
            this->pelec->calculate_weights();
            // calculate Mi from existed becp
            for (int ik = 0; ik < nk; ik++)
            {
                const std::complex<double>* becp = &becp_tmp[ik * size_becp];
                // becp(nbands*npol , nkb)
                // mag = wg * \sum_{nh}becp * becp
                for (int ib = 0; ib < nbands; ib++)
                {
                    const double weight = this->pelec->wg(ik, ib);
                    int begin_ih = 0;
                    for (int iat = 0; iat < this->Mi_.size(); iat++)
                    {
                        const int nh = nh_iat[iat];
                        std::complex<double> occ[4]
                            = {ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO};
                        for (int ih = 0; ih < nh; ih++)
                        {
                            const int index = ib * npol * nkb + begin_ih + ih;
                            occ[0] += conj(becp[index]) * becp[index];
                            occ[1] += conj(becp[index]) * becp[index + nkb];
                            occ[2] += conj(becp[index + nkb]) * becp[index];
                            occ[3] += conj(becp[index + nkb]) * becp[index + nkb];
                        }
                        // occ has been reduced and calculate mag
                        this->Mi_[iat].x += weight * (occ[1] + occ[2]).real();
                        this->Mi_[iat].y += weight * (occ[1] - occ[2]).imag();
                        this->Mi_[iat].z += weight * (occ[0] - occ[3]).real();
                        begin_ih += nh;
                    }
                }
            }
            Parallel_Reduce::reduce_double_allpool(GlobalV::KPAR,
                                                   GlobalV::NPROC_IN_POOL,
                                                   &(this->Mi_[0][0]),
                                                   3 * this->Mi_.size());
            // for(int i = 0; i < this->Mi_.size(); i++)
            //{
            //     std::cout<<"atom"<<i<<": "<<" mag: "<<this->Mi_[i].x<<" "<<this->Mi_[i].y<<" "<<this->Mi_[i].z<<"
            //     "<<this->lambda_[i].x<<" "<<this->lambda_[i].y<<" "<<this->lambda_[i].z<<std::endl;
            // }
        }
    }
    ModuleBase::timer::tick("SpinConstrain", "cal_mw_from_lambda");
}

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::update_psi_charge(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve)
{
    ModuleBase::TITLE("spinconstrain::SpinConstrain", "update_psi_charge");
    ModuleBase::timer::tick("spinconstrain::SpinConstrain", "update_psi_charge");
#ifdef __LCAO
    if (PARAM.inp.basis_type == "lcao")
    {
        psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
        this->pelec->psiToRho(*psi_t);
    }
    else
#endif
    {
        int size_becp = 0;
        std::vector<std::complex<double>> becp_tmp;
        int nk = 0;
        int nkb = 0;
        int nbands = 0;
        int npol = 0;
        const int* nh_iat = nullptr;
        if (PARAM.inp.device == "cpu")
        {
            psi::Psi<std::complex<double>>* psi_t = static_cast<psi::Psi<std::complex<double>>*>(this->psi);
            hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>*>(this->p_hamilt);
            auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::get_instance();
            nbands = psi_t->get_nbands();
            npol = psi_t->npol;
            nkb = onsite_p->get_tot_nproj();
            nk = psi_t->get_nk();
            nh_iat = &onsite_p->get_nh(0);
            size_becp = nbands * nkb * npol;
            becp_tmp.resize(size_becp * nk);
            std::vector<std::complex<double>> h_tmp(nbands * nbands), s_tmp(nbands * nbands);
            assert(this->sub_h_save != nullptr);
            assert(this->sub_s_save != nullptr);
            assert(this->becp_save != nullptr);
            for (int ik = 0; ik < nk; ++ik)
            {
                std::complex<double>* h_k = this->sub_h_save + ik * nbands * nbands;
                std::complex<double>* s_k = this->sub_s_save + ik * nbands * nbands;
                std::complex<double>* becp_k = this->becp_save + ik * size_becp;

                psi_t->fix_k(ik);
                memcpy(h_tmp.data(), h_k, sizeof(std::complex<double>) * nbands * nbands);
                memcpy(s_tmp.data(), s_k, sizeof(std::complex<double>) * nbands * nbands);
                const int sign = this->pelec->klist->isk[ik] == 0 ? 1 : -1;
                this->calculate_delta_hcc(h_tmp.data(), becp_k, delta_lambda, nbands, nkb, nh_iat, sign);
                hsolver::DiagoIterAssist<std::complex<double>>::diag_subspace_psi(h_tmp.data(),
                                                                                s_tmp.data(),
                                                                                nbands,
                                                                                psi_t[0],
                                                                                &this->pelec->ekb(ik, 0));
            }

            delete[] this->sub_h_save;
            delete[] this->sub_s_save;
            delete[] this->becp_save;
            this->sub_h_save = nullptr;
            this->sub_s_save = nullptr;
            this->becp_save = nullptr;

            if(pw_solve)
            {
                hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_CPU> hsolver_pw_obj(this->pw_wfc_,
                                                 PARAM.inp.calculation,
                                                 PARAM.inp.basis_type,
                                                 PARAM.inp.ks_solver,
                                                 PARAM.inp.use_paw,
                                                 PARAM.globalv.use_uspp,
                                                 PARAM.inp.nspin,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::SCF_ITER,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::PW_DIAG_NMAX,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::PW_DIAG_THR,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::need_subspace);

                hsolver_pw_obj.solve(hamilt_t,
                         psi_t[0],
                         this->pelec,
                         this->pelec->ekb.c,
                         GlobalV::RANK_IN_POOL,
                         GlobalV::NPROC_IN_POOL,
                         false,
                         this->tpiba,
                         this->get_nat());
            }
            else
            {// update charge density only
                this->pelec->psiToRho(*psi_t);
            }
        }
#if ((defined __CUDA) || (defined __ROCM))
        else
        {
            base_device::DEVICE_GPU* ctx = {};
            base_device::DEVICE_CPU* cpu_ctx = {};
            psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*>(this->psi);
            hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>* hamilt_t = static_cast<hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>*>(this->p_hamilt);
            assert(this->sub_h_save != nullptr);
            assert(this->sub_s_save != nullptr);
            assert(this->becp_save != nullptr);
            auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();
            nbands = psi_t->get_nbands();
            npol = psi_t->npol;
            nkb = onsite_p->get_tot_nproj();
            nk = psi_t->get_nk();
            nh_iat = &onsite_p->get_nh(0);
            size_becp = nbands * nkb * npol;
            const bool davidson_solve = pw_solve && PARAM.inp.ks_solver == "dav_subspace";
            // Davidson uses the Euclidean overlap. The saved generalized
            // overlap has the same meaning only for norm-conserving inputs.
            const bool reuse_subspace = davidson_solve && !PARAM.inp.use_paw && !PARAM.globalv.use_uspp;
            using Solver = hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_GPU>;
            const Solver::BasisObserver* observer = nullptr;
#if defined(__CUDA)
            std::unique_ptr<NonlocalBasisObserver> projection_observer;
            auto* nonlocal = reuse_subspace && GlobalV::NPROC_IN_POOL == 1
                ? find_nonlocal_operator(hamilt_t) : nullptr;
            const int projection_rows = nonlocal ? nonlocal->get_ppcell()->nkb * npol : 0;
            const size_t projection_count = size_t(projection_rows) * nbands;
            const bool initial_projections_valid = projection_rows > 0
                && this->nonlocal_projections.size() == size_t(nk)
                && std::all_of(this->nonlocal_projections.begin(), this->nonlocal_projections.end(),
                               [projection_count](const std::vector<Complex>& entry) { return entry.size() == projection_count; });
            if (initial_projections_valid)
            {
                projection_observer.reset(new NonlocalBasisObserver(nonlocal, nbands, npol,
                    nbands * PARAM.inp.pw_diag_ndim, this->nonlocal_projections));
                observer = projection_observer->callbacks();
            }
#endif
            if (observer == nullptr) { this->invalidate_nonlocal_projections(); }
            hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_GPU>::InitialHPsiFunc initial_hpsi;
            if (reuse_subspace)
            {
                // The local potential and initial wavefunctions stay fixed
                // throughout the constraint loop. Only the low-rank onsite
                // field changes.
                std::vector<ModuleBase::Vector3<double>> difference(this->lambda_.size());
                for (size_t atom = 0; atom < difference.size(); ++atom)
                {
                    difference[atom] = this->lambda_[atom] - this->subspace_lambda[atom];
                }
                for (int ik = 0; ik < nk; ++ik)
                {
                    const int sign = this->pelec->klist->isk[ik] == 0 ? 1 : -1;
                    this->calculate_delta_hcc(this->sub_h_save + ik * nbands * nbands,
                                             this->becp_save + ik * size_becp, difference.data(),
                                             nbands, nkb, nh_iat, sign);
                }
#if defined(__CUDA)
                if (this->hpsi_cache)
                {
                    const auto coefficients = pack_constraint_field(difference.size(), nkb, nh_iat, difference.data());
                    const int source_leading = psi_t->get_nbasis();
                    // The cached action uses the same entry basis and local potential.
                    // Only the onsite constraint field needs a low-rank correction.
                    initial_hpsi = [=](int ik, std::complex<double>* output, int leading) {
                        using Complex = std::complex<double>;
                        using GPU = base_device::DEVICE_GPU;
                        using CPU = base_device::DEVICE_CPU;
                        const size_t offset = size_t(ik) * nbands * source_leading;
                        cudaErrcheck(cudaMemcpy2D(output, size_t(leading) * sizeof(Complex),
                            this->hpsi_cache.get() + offset, size_t(source_leading) * sizeof(Complex),
                            size_t(leading) * sizeof(Complex), nbands, cudaMemcpyHostToDevice));
                        ct::Tensor field(ct::DataType::DT_DOUBLE, ct::DeviceType::GpuDevice, {nkb, 3});
                        ct::Tensor weighted(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice,
                                            {nbands * npol, nkb});
                        base_device::memory::synchronize_memory_op<double, GPU, CPU>()(
                            ctx, nullptr, field.data<double>(), coefficients.data(), coefficients.size());
                        const int sign = this->pelec->klist->isk[ik] == 0 ? 1 : -1;
                        hamilt::spin_constraint_projectors_gpu(nbands, npol, nkb, sign, field.data<double>(),
                            this->becp_save + size_t(ik) * size_becp, weighted.data<Complex>());
                        onsite_p->get_fs_tools()->prepare_projectors(ik);
                        const int npw = onsite_p->get_npw();
                        hsolver::gemm_op<Complex, GPU>()(ctx, 'N', 'N', npw, nbands * npol, nkb,
                            &ModuleBase::ONE, onsite_p->get_tab_atomic(), npw, weighted.data<Complex>(), nkb,
                            &ModuleBase::ONE, output, leading / npol);
                    };
                }
#endif
            }
            else if (!davidson_solve)
            {
                std::complex<double>* h_tmp = nullptr;
                std::complex<double>* s_tmp = nullptr;
                base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, h_tmp, nbands * nbands);
                base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, s_tmp, nbands * nbands);
                for (int ik = 0; ik < nk; ++ik)
                {
                    std::complex<double>* h_k = this->sub_h_save + ik * nbands * nbands;
                    std::complex<double>* s_k = this->sub_s_save + ik * nbands * nbands;
                    std::complex<double>* becp_k = this->becp_save + ik * size_becp;

                    psi_t->fix_k(ik);
                    const int sign = this->pelec->klist->isk[ik] == 0 ? 1 : -1;
                    base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(ctx, ctx, h_tmp, h_k, nbands * nbands);
                    base_device::memory::synchronize_memory_op<std::complex<double>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(ctx, ctx, s_tmp, s_k, nbands * nbands);
                    this->calculate_delta_hcc(h_tmp, becp_k, delta_lambda, nbands, nkb, nh_iat, sign);
                    hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::diag_subspace_psi(h_tmp,
                                                                                    s_tmp,
                                                                                    nbands,
                                                                                    psi_t[0],
                                                                                    &this->pelec->ekb(ik, 0));
                }

                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, h_tmp);
                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, s_tmp);
            }
            if (!reuse_subspace)
            {
                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, sub_h_save);
                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, sub_s_save);
                this->sub_h_save = nullptr;
                this->sub_s_save = nullptr;
            }
            if (!initial_hpsi)
            {
                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, becp_save);
                this->becp_save = nullptr;
            }

            if(pw_solve)
            {
                hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_GPU> hsolver_pw_obj(this->pw_wfc_,
                                                 PARAM.inp.calculation,
                                                 PARAM.inp.basis_type,
                                                 PARAM.inp.ks_solver,
                                                 PARAM.inp.use_paw,
                                                 PARAM.globalv.use_uspp,
                                                 PARAM.inp.nspin,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::SCF_ITER,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::PW_DIAG_NMAX,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::PW_DIAG_THR,
                                                 hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_GPU>::need_subspace);

                hsolver_pw_obj.solve(hamilt_t,
                         psi_t[0],
                         this->pelec,
                         this->pelec->ekb.c,
                         GlobalV::RANK_IN_POOL,
                         GlobalV::NPROC_IN_POOL,
                         false,
                         this->tpiba,
                         this->get_nat(),
                         this->sub_h_save,
                         this->sub_s_save,
                         initial_hpsi,
                         observer);
            }
            else
            {// update charge density only
                reinterpret_cast<elecstate::ElecStatePW<std::complex<double>, base_device::DEVICE_GPU>*>(this->pelec)->psiToRho(*psi_t);
            }
            if (reuse_subspace)
            {
                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, sub_h_save);
                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, sub_s_save);
                this->sub_h_save = nullptr;
                this->sub_s_save = nullptr;
            }
            if (initial_hpsi)
            {
                base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(ctx, becp_save);
                this->becp_save = nullptr;
            }
            this->subspace_lambda.clear();
            
        }
#endif       
    }
    ModuleBase::timer::tick("spinconstrain::SpinConstrain", "update_psi_charge");
}
