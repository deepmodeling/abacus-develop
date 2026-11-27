#include "source_pw/module_rt/td_orbital_ops.h"

#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_hsolver/module_linalg/orb_buffer.h"

#include <iomanip>
#include <iostream>
#include <sstream>

namespace hsolver
{
template <typename T, typename Device>
TDOrbOps<T, Device>::TDOrbOps(int npol,
                              OrthMethod orth_method,
                              bool out_stat,
                              const std::vector<int>& global_k,
                              const diag_comm_info& comm,
                              std::ostream& log)
    : comm_(comm), npol_(npol), orth_method_(orth_method), out_stat_(out_stat), global_k_(global_k), algebra_(comm),
      orthonormal_(comm, algebra_), log_(log), band_products_(comm, 9)
{
    log_ << " PW RT-TDDFT orthonormalization: " << orth_method_name(orth_method_) << '\n';
}

template <typename T, typename Device>
void TDOrbOps<T, Device>::prepare_buffers(int nbands, int nbasis)
{
    const int64_t size = static_cast<int64_t>(nbands) * nbasis;
    linear_buffer<T, Device>(&clean_, size);
    linear_buffer<T, Device>(&hpsi_, size);
}

template <typename T, typename Device>
void TDOrbOps<T, Device>::correct(psi::Psi<T, Device>* current, int step, int iter, bool detailed_output)
{
    if (out_stat_)
    {
        orth_norms_.resize(current->get_nk());
    }
    zero_padding(current);
    for (int ik = 0; ik < current->get_nk(); ++ik)
    {
        current->fix_k(ik);
        const int ld = current->get_nbasis();
        const int dim = npol_ == 2 ? ld : current->get_ngk(ik);
        const int bands = current->get_nbands();
        T* data = current->get_pointer();
        correct_orbitals(data, ld, dim, bands, ik, step, iter, detailed_output);
    }
}

template <typename T, typename Device>
void TDOrbOps<T, Device>::cal_band_energy(HSOperator<T, Device>& op, const psi::Psi<T, Device>& current, ModuleBase::matrix* energies)
{
    ModuleBase::timer::start("TDOrbOps", "cal_band_energy");
    const int nband = current.get_nbands();
    const int ld = current.get_nbasis();
    prepare_buffers(nband, ld);
    std::vector<T> expectations(nband);
    T* hpsi = hpsi_.template data<T>();
    for (int ik = 0; ik < current.get_nk(); ++ik)
    {
        op.update_k(ik);
        current.fix_k(ik);
        const int npw = current.get_ngk(ik);
        const int dim = npol_ == 2 ? ld : npw;
        const OrbLayout layout{ld, npw, npol_, nband};
        const T* data = clean_orbs<T, Device>(current.get_pointer(), layout, &clean_);
        op.hpsi(data, hpsi, ld, nband);
        if (npol_ == 2)
        {
            zero_pad<T, Device>(hpsi, layout);
        }
        band_products_.dot(ld, dim, nband, data, hpsi, expectations.data());
        for (int band = 0; band < nband; ++band)
        {
            (*energies)(ik, band) = std::real(expectations[band]);
        }
    }
    ModuleBase::timer::end("TDOrbOps", "cal_band_energy");
}

template <typename T, typename Device>
std::vector<std::complex<double>> TDOrbOps<T, Device>::cal_eij(HSOperator<T, Device>& op, const psi::Psi<T, Device>& current, const int ik)
{
    ModuleBase::timer::start("TDOrbOps", "cal_eij");
    const int nband = current.get_nbands();
    const int ld = current.get_nbasis();
    const int npw = current.get_ngk(ik);
    const int dim = npol_ == 2 ? ld : npw;
    prepare_buffers(nband, ld);
    op.update_k(ik);
    current.fix_k(ik);
    const OrbLayout layout{ld, npw, npol_, nband};
    const T* psi = clean_orbs<T, Device>(current.get_pointer(), layout, &clean_);
    T* hpsi = hpsi_.template data<T>();
    op.hpsi(psi, hpsi, ld, nband);
    if (npol_ == 2)
    {
        zero_pad<T, Device>(hpsi, layout);
    }
    std::vector<std::complex<double>> result = algebra_.cross(ld, dim, nband, nband, psi, hpsi);
    ModuleBase::timer::end("TDOrbOps", "cal_eij");
    return result;
}

template <typename T, typename Device>
void TDOrbOps<T, Device>::reset_orth_stats()
{
    if (out_stat_)
    {
        orth_stats_ = TDOrthStats();
    }
    warned_events_ = 0;
}

template <typename T, typename Device>
void TDOrbOps<T, Device>::correct_orbitals(T* current, int ld, int dim, int bands, int ik, int istep, int iter, bool detailed_output)
{
    ModuleBase::timer::start("TDOrbOps", "correct_orbitals");
    const OrthResult result = orthonormal_.apply(current, ld, dim, bands, orth_method_, out_stat_);
    record_orth(result, ik, istep, iter);
    if (detailed_output && result.gram_checked && log_.good())
    {
        const int evolution_step = istep + 1;
        std::ostringstream record;
        record << std::scientific << std::setprecision(6) << " PW RT-TDDFT orth: evolution_step=" << evolution_step << " scf_iter=" << iter
               << " global_k=" << global_k_.at(ik) + 1 << " requested=" << orth_method_name(orth_method_)
               << " orth_before=" << result.before << " orth_after=" << result.after << '\n';
        log_ << record.str();
    }
    ModuleBase::timer::end("TDOrbOps", "correct_orbitals");
}

template <typename T, typename Device>
void TDOrbOps<T, Device>::record_orth(const OrthResult& result, int ik, int istep, int iter)
{
    if (result.status == OrthStatus::failed)
    {
        const char* stage = istep == 0 ? "initialization" : "propagation";
        std::ostringstream message;
        const int evolution_step = istep + 1;
        message << std::setprecision(16) << " PW RT-TDDFT orbital validation failed (" << stage << "): evolution_step=" << evolution_step
                << " scf_iter=" << iter << " global_k=" << global_k_.at(ik) + 1 << " requested=" << orth_method_name(orth_method_) << '\n'
                << "   last_attempted=" << orth_method_name(result.actual) << " passes=" << result.passes
                << " fallbacks=" << result.fallbacks << '\n'
                << "   " << orth_failure_name(result.failure) << "; " << result.reason;
        if (result.gram_checked)
        {
            message << '\n' << "   before=" << result.before << " after=" << result.after << " tolerance=" << orth_tolerance<T>();
        }
        // stdout is disabled on non-world-root ranks; a failing pool must still report its error.
        if (comm_.rank == 0)
        {
            std::cerr << message.str() << std::endl;
        }
#ifdef __MPI
        // Let the pool root flush its diagnostic before another rank's exit stops the MPI job.
        double diagnostic_written = 1.0;
        Parallel_Common::bcast_data(&diagnostic_written, 1, comm_.comm, 0);
#endif
        ModuleBase::WARNING_QUIT("TDOrbOps", message.str());
    }
    report_orth_warning(result, ik, istep, iter);
    if (out_stat_)
    {
        orth_norms_[ik] = result.norms;
        if (result.gram_checked)
        {
            orth_stats_.before = std::max(orth_stats_.before, result.before);
            orth_stats_.after = std::max(orth_stats_.after, result.after);
        }
    }
}

template <typename T, typename Device>
void TDOrbOps<T, Device>::report_orth_warning(const OrthResult& result, int ik, int istep, int iter)
{
    if (comm_.rank != 0)
    {
        return;
    }
    unsigned int events = 0;
    if (result.fallbacks > 0)
    {
        events |= 1;
    }
    if (result.rejected > 0)
    {
        events |= 2;
    }
    const unsigned int fresh_events = events & ~warned_events_;
    if (fresh_events == 0)
    {
        return;
    }
    warned_events_ |= events;
    const int evolution_step = istep + 1;
    std::ostringstream message;
    message << std::setprecision(16) << " PW RT-TDDFT orbital warning: evolution_step=" << evolution_step << " scf_iter=" << iter;
    message << " global_k=" << global_k_.at(ik) + 1;
    message << " requested=" << orth_method_name(orth_method_) << '\n'
            << "   before=" << result.before << " after=" << result.after << " tolerance=" << orth_tolerance<T>() << '\n'
            << "   ";
    message << "last_attempted=" << orth_method_name(result.actual) << ". ";
    if (fresh_events & 1)
    {
        message << "A fallback was attempted. ";
    }
    if (fresh_events & 2)
    {
        message << "A candidate was rejected. ";
    }
    message << "The retained state satisfies the tolerance." << '\n'
            << "   " << result.reason << '\n'
            << "   Further events of these types in this pool are suppressed for this electronic evolution step.";
    const std::string text = message.str();
    std::cerr << text << std::endl;
    if (log_.good())
    {
        log_ << text << std::endl;
    }
}

template <typename T, typename Device>
void TDOrbOps<T, Device>::zero_padding(psi::Psi<T, Device>* current)
{
    ModuleBase::timer::start("TDOrbOps", "zero_padding");
    if (npol_ != 2)
    {
        ModuleBase::timer::end("TDOrbOps", "zero_padding");
        return;
    }
    const int old_k = current->get_current_k();
    for (int ik = 0; ik < current->get_nk(); ++ik)
    {
        current->fix_k(ik);
        const int npw = current->get_ngk(ik);
        const int ld = current->get_nbasis();
        const int bands = current->get_nbands();
        const OrbLayout layout{ld, npw, npol_, bands};
        zero_pad<T, Device>(current->get_pointer(), layout);
    }
    current->fix_k(old_k);
    ModuleBase::timer::end("TDOrbOps", "zero_padding");
}

template <typename T, typename Device>
bool TDOrbOps<T, Device>::orth_initial(psi::Psi<T, Device>* current, int iter)
{
    ModuleBase::timer::start("TDOrbOps", "orth_initial");
    bool changed = false;
    if (orth_method_ != OrthMethod::none)
    {
        if (out_stat_)
        {
            // Initial diagnostics describe this SCF state, not earlier discarded iterates.
            orth_stats_ = TDOrthStats();
            orth_norms_.resize(current->get_nk());
        }
        for (int ik = 0; ik < current->get_nk(); ++ik)
        {
            current->fix_k(ik);
            const int dim = npol_ == 2 ? current->get_nbasis() : current->get_ngk(ik);
            const OrthResult result
                = orthonormal_.apply(current->get_pointer(), current->get_nbasis(), dim, current->get_nbands(), orth_method_, out_stat_);
            record_orth(result, ik, 0, iter);
            changed = changed || result.status == OrthStatus::accepted;
        }
    }
    ModuleBase::timer::end("TDOrbOps", "orth_initial");
    return changed;
}

template <typename T, typename Device>
void TDOrbOps<T, Device>::check_initial(const psi::Psi<T, Device>& current, int iter)
{
    ModuleBase::timer::start("TDOrbOps", "check_initial");
    const bool full_gram = orth_method_ != OrthMethod::none;
    const int ld = current.get_nbasis();
    const int bands = current.get_nbands();
    prepare_buffers(bands, ld);
    if (out_stat_)
    {
        orth_norms_.resize(current.get_nk());
    }
    for (int ik = 0; ik < current.get_nk(); ++ik)
    {
        current.fix_k(ik);
        const int npw = current.get_ngk(ik);
        const int dim = npol_ == 2 ? ld : npw;
        const OrbLayout layout{ld, npw, npol_, bands};
        const T* data = clean_orbs<T, Device>(current.get_pointer(), layout, &clean_);
        const OrthResult result = orthonormal_.inspect(data, current.get_nbasis(), dim, current.get_nbands(), full_gram, out_stat_);
        record_orth(result, ik, 0, iter);
    }
    ModuleBase::timer::end("TDOrbOps", "check_initial");
}

template <typename T, typename Device>
double TDOrbOps<T, Device>::wave_electrons(const ModuleBase::matrix& occupations) const
{
    double electrons = 0.0;
    for (std::size_t ik = 0; ik < orth_norms_.size(); ++ik)
    {
        for (std::size_t band = 0; band < orth_norms_[ik].size(); ++band)
        {
            electrons += occupations(ik, band) * orth_norms_[ik][band];
        }
    }
    return electrons;
}

template class TDOrbOps<std::complex<float>, base_device::DEVICE_CPU>;
template class TDOrbOps<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class TDOrbOps<std::complex<float>, base_device::DEVICE_GPU>;
template class TDOrbOps<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace hsolver
