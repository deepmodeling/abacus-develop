#include "sto_dos.h"

#include "source_base/parallel_comm.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_io/module_parameter/parameter.h"
#include "sto_dos_quality.h"
#include "sto_tool.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace
{
/// Single entry point for the DOS diagnostics of this file.  Keeping the
/// running-log reference here (instead of at every call site) keeps the global
/// dependency of this module to one place; none of these messages aborts a run.
void dos_log(const std::string& message, const bool also_stdout)
{
    GlobalV::ofs_running << message << std::endl;
    if (also_stdout)
    {
        std::cout << message << std::endl;
    }
}
} // namespace

template <typename FPTYPE, typename Device>
Sto_DOS<FPTYPE, Device>::~Sto_DOS()
{
}

template <typename FPTYPE, typename Device>
Sto_DOS<FPTYPE, Device>::Sto_DOS(ModulePW::PW_Basis_K* p_wfcpw_in,
                                 K_Vectors* p_kv_in,
                                 elecstate::ElecState* p_elec_in,
                                 psi::Psi<std::complex<double>>* p_psi_in,
                                 hamilt::Hamilt<std::complex<double>>* p_hamilt_in,
                                 StoChe<FPTYPE, Device>& stoche,
                                 Stochastic_WF<std::complex<double>, base_device::DEVICE_CPU>* p_stowf_in)
{
    this->p_wfcpw = p_wfcpw_in;
    this->p_kv = p_kv_in;
    this->p_elec = p_elec_in;
    this->p_psi = p_psi_in;
    this->p_hamilt = p_hamilt_in;
    this->p_hamilt_sto = static_cast<StoHamiltPW<std::complex<double>>*>(p_hamilt_in);
    this->p_stowf = p_stowf_in;
    this->nbands_ks = p_psi_in->get_nbands();
    // Stochastic_WF::nchi is now initialized by allocate_chi0() (MDFT) and
    // init_com_orbitals() (complete basis).  Keep a safe fallback for old
    // restart paths where only nchip is available.
    this->nbands_sto = (p_stowf_in->nchi > 0) ? p_stowf_in->nchi : p_stowf_in->nchip_max;
    this->method_sto = stoche.method_sto;
    this->stofunc.set_E_range((double*)&stoche.emin_sto, (double*)&stoche.emax_sto);
}

template <typename FPTYPE, typename Device>
void Sto_DOS<FPTYPE, Device>::decide_param(const int& dos_nche,
                                           const double& emin_sto,
                                           const double& emax_sto,
                                           const bool& dos_setemin,
                                           const bool& dos_setemax,
                                           const double& dos_emin_ev,
                                           const double& dos_emax_ev,
                                           const double& dos_scale)
{
    this->dos_nche = dos_nche;
    if (this->dos_nche < 2)
    {
        ModuleBase::WARNING_QUIT("Sto_DOS", "dos_nche should be at least 2.");
    }

    // If emin_sto/emax_sto are left at zero, keep the bounds already found by
    // the SCF range check instead of resetting them to the bare input defaults.
    // The repeated reset discarded a valid converged range and made the DOS
    // range check depend on the random test vector alone.
    double try_emin = emin_sto;
    double try_emax = emax_sto;
    if (try_emin == 0.0)
    {
        try_emin = *p_hamilt_sto->emin;
    }
    if (try_emax == 0.0)
    {
        try_emax = *p_hamilt_sto->emax;
    }

    check_che_op<FPTYPE, Device>()(
        this->dos_nche,
        try_emin,
        try_emax,
        this->nbands_sto,
        this->p_kv,
        reinterpret_cast<Stochastic_WF<std::complex<FPTYPE>, Device>*>(this->p_stowf),
        reinterpret_cast<StoHamiltPW<std::complex<FPTYPE>, Device>*>(this->p_hamilt_sto));

    // The range determined by the SCF check only learns the high-energy end
    // reliably.  Run the pure, unit-tested sanitiser (module
    // sto_dos_quality): it covers the known deterministic eigenvalues and the
    // plane-wave spectrum top, so an occupied state outside the interval (T_n
    // grows exponentially there) and an underestimated upper bound are both
    // removed.  Healthy ranges are returned bit for bit unchanged.
    {
        std::vector<double> det_eigs;
        if (this->p_elec != nullptr)
        {
            const auto& ekb = this->p_elec->ekb;
            for (int ik = 0; ik < ekb.nr; ++ik)
            {
                for (int ib = 0; ib < ekb.nc; ++ib)
                {
                    const double e = ekb(ik, ib);
                    if (std::isfinite(e))
                    {
                        det_eigs.push_back(e);
                    }
                }
            }
        }
        const ModuleIO::CheRange fixed
            = ModuleIO::sanitize_che_range(*p_hamilt_sto->emin, *p_hamilt_sto->emax, PARAM.inp.ecutwfc, det_eigs);
        if (fixed.emin_repaired || fixed.emax_repaired || fixed.degenerate)
        {
            std::ostringstream range_msg;
            range_msg << " DOS: Chebyshev range repaired from [" << *p_hamilt_sto->emin << ", " << *p_hamilt_sto->emax
                      << "] to [" << fixed.emin << ", " << fixed.emax << "] Ry";
            dos_log(range_msg.str(), false);
        }
        *p_hamilt_sto->emin = fixed.emin;
        *p_hamilt_sto->emax = fixed.emax;
    }

    if (dos_setemax)
    {
        this->emax = dos_emax_ev;
    }
    else
    {
        this->emax = *p_hamilt_sto->emax * ModuleBase::Ry_to_eV;
    }
    if (dos_setemin)
    {
        this->emin = dos_emin_ev;
    }
    else
    {
        this->emin = *p_hamilt_sto->emin * ModuleBase::Ry_to_eV;
    }

    if (!dos_setemax && !dos_setemin)
    {
        double delta = (emax - emin) * dos_scale;
        this->emax = emax + delta / 2.0;
        this->emin = emin - delta / 2.0;
    }

    if (!std::isfinite(this->emin) || !std::isfinite(this->emax) || this->emax <= this->emin)
    {
        ModuleBase::WARNING_QUIT("Sto_DOS", "Invalid DOS energy range: dos_emax_ev must be larger than dos_emin_ev.");
    }

    // The Chebyshev interval is a property of the spectrum and is kept as it is:
    // widening it to fit the requested output window would reduce the resolution
    // of a fixed-order expansion and could hide a peak that a narrower interval
    // resolves.  A requested window that sticks out of the interval is therefore
    // clipped to it (the `dos_scale` widening above makes that routine) and the
    // change is reported, because energies outside the interval cannot be
    // represented by the expansion at all.
    const double interval_emin_ev = *p_hamilt_sto->emin * ModuleBase::Ry_to_eV;
    const double interval_emax_ev = *p_hamilt_sto->emax * ModuleBase::Ry_to_eV;
    const double requested_emin_ev = this->emin;
    const double requested_emax_ev = this->emax;
    const double clipped_emin_ev = std::max(requested_emin_ev, interval_emin_ev);
    const double clipped_emax_ev = std::min(requested_emax_ev, interval_emax_ev);
    // The `dos_scale` widening above intentionally makes the window a little wider
    // than the spectrum, so a sub-percent clip is expected and is not reported;
    // only a window that really asks for energies outside the interval is.
    const double interval_span_ev = interval_emax_ev - interval_emin_ev;
    const double dropped_ev
        = std::max(requested_emin_ev - clipped_emin_ev, clipped_emax_ev - requested_emax_ev);
    const bool report_clip = dropped_ev > 0.01 * interval_span_ev;
    if (clipped_emax_ev > clipped_emin_ev)
    {
        // Always clip, so that |t| <= 1 holds by construction for every output
        // energy.  Only a clip that removes a visible part of the window is
        // reported, because `dos_scale` intentionally widens the window by a
        // fraction of a percent.
        if (report_clip)
        {
            std::ostringstream clip_msg;
            clip_msg << " DOS: output window clipped from [" << requested_emin_ev << ", " << requested_emax_ev << "] to ["
                     << clipped_emin_ev << ", " << clipped_emax_ev << "] eV to stay inside the Chebyshev interval ["
                     << interval_emin_ev << ", " << interval_emax_ev << "] eV";
            dos_log(clip_msg.str(), false);
        }
        this->emin = clipped_emin_ev;
        this->emax = clipped_emax_ev;
    }

    std::ostringstream range_msg;
    range_msg << " DOS: Chebyshev range = [" << *p_hamilt_sto->emin << ", " << *p_hamilt_sto->emax
              << "] Ry; output window = [" << this->emin << ", " << this->emax << "] eV";
    dos_log(range_msg.str(), false);
}

template <typename FPTYPE, typename Device>
void Sto_DOS<FPTYPE, Device>::caldos(const double sigmain, const double de, const int npart)
{
    ModuleBase::TITLE("Sto_DOS", "caldos");
    ModuleBase::timer::start("Sto_DOS", "caldos");
    std::cout << "=========================" << std::endl;
    std::cout << "###Calculating Dos....###" << std::endl;
    std::cout << "=========================" << std::endl;
    ModuleBase::Chebyshev<double> che(dos_nche);
    const int nk = p_kv->get_nks();
    const int npwx = p_wfcpw->npwk_max;

    std::vector<double> spolyv;
    std::vector<std::complex<double>> allorderchi;
    if (this->method_sto == 1)
    {
        spolyv.resize(dos_nche, 0);
    }
    else
    {
        spolyv.resize(dos_nche * dos_nche, 0);
        int nchip_new = ceil((double)this->p_stowf->nchip_max / npart);
        allorderchi.resize(nchip_new * npwx * dos_nche);
    }
    ModuleBase::timer::start("Sto_DOS", "Tracepoly");
    std::cout << "1. TracepolyA:" << std::endl;
    for (int ik = 0; ik < nk; ik++)
    {
        std::cout << "ik: " << ik + 1 << std::endl;
        if (nk > 1)
        {
            this->p_hamilt->updateHk(ik);
        }
        const int npw = p_kv->ngk[ik];
        const int nchipk = this->p_stowf->nchip[ik];

        std::complex<double>* pchi;
        if (PARAM.inp.nbands > 0)
        {
            p_stowf->chiortho->fix_k(ik);
            pchi = p_stowf->chiortho->get_pointer();
        }
        else
        {
            p_stowf->chi0->fix_k(ik);
            pchi = p_stowf->chi0->get_pointer();
        }
        auto hchi_norm = std::bind(&StoHamiltPW<std::complex<double>>::hPsi_norm,
                                   p_hamilt_sto,
                                   std::placeholders::_1,
                                   std::placeholders::_2,
                                   std::placeholders::_3);
        if (this->method_sto == 1)
        {
            che.tracepolyA(hchi_norm, pchi, npw, npwx, nchipk);
            for (int i = 0; i < dos_nche; ++i)
            {
                spolyv[i] += che.polytrace[i] * p_kv->wk[ik] / 2;
            }
        }
        else
        {
            int N = dos_nche;
            double kweight = p_kv->wk[ik] / 2;
            char trans = 'T';
            char normal = 'N';
            double one = 1;
            for (int ipart = 0; ipart < npart; ++ipart)
            {
                int nchipk_new = nchipk / npart;
                int start_nchipk = ipart * nchipk_new + nchipk % npart;
                if (ipart < nchipk % npart)
                {
                    nchipk_new++;
                    start_nchipk = ipart * nchipk_new;
                }
                ModuleBase::GlobalFunc::ZEROS(allorderchi.data(), nchipk_new * npwx * dos_nche);
                std::complex<double>* tmpchi = pchi + start_nchipk * npwx;
                che.calpolyvec_complex(hchi_norm, tmpchi, allorderchi.data(), npw, npwx, nchipk_new);
                double* vec_all = (double*)allorderchi.data();
                int LDA = npwx * nchipk_new * 2;
                int M = npwx * nchipk_new * 2;
                BlasConnector::gemm(normal, trans, N, N, M, kweight, vec_all, LDA, vec_all, LDA, one, spolyv.data(), N);
            }
        }
    }

    allorderchi.clear();
    allorderchi.shrink_to_fit();

    std::ofstream ofsdos;
    if (de <= 0.0)
    {
        ModuleBase::WARNING_QUIT("Sto_DOS", "dos_edelta_ev must be positive.");
    }
    int ndos = int((emax - emin) / de) + 1;
    if (ndos < 1)
    {
        ModuleBase::WARNING_QUIT("Sto_DOS", "The DOS output window is empty; check dos_emin_ev/dos_emax_ev.");
    }
    this->stofunc.sigma = sigmain / ModuleBase::Ry_to_eV;

    // Report DOS settings that cannot resolve the curve: grid spacing versus
    // Gaussian broadening, and Chebyshev order versus range/sigma.  The checks
    // are pure functions (unittests/test_sto_dos_quality.cpp) and only report;
    // a long run is never aborted by a quality suggestion.
    {
        const double range_ev = (*this->stofunc.Emax - *this->stofunc.Emin) * ModuleBase::Ry_to_eV;
        const double sigma_ev = std::abs(sigmain);
        for (const std::string& msg: ModuleIO::check_dos_parameters(range_ev, sigma_ev, de, this->dos_nche))
        {
            dos_log(" WARNING: DOS: " + msg, true);
        }
    }
    ModuleBase::timer::end("Sto_DOS", "Tracepoly");

    std::cout << "2. Dos:" << std::endl;
    ModuleBase::timer::start("Sto_DOS", "DOS Loop");
    // For very small output windows ndos/10 is zero; the old modulus then
    // divided by zero and terminated the run with SIGFPE.
    int n10 = std::max(1, ndos / 10);
    int percent = 10;
    std::vector<double> sto_dos(ndos);
    std::vector<double> ks_dos(ndos);
    std::vector<double> error(ndos);
    for (int ie = 0; ie < ndos; ++ie)
    {
        double tmpks = 0;
        double tmpsto = 0;
        this->stofunc.targ_e = (emin + ie * de) / ModuleBase::Ry_to_eV;
        if (this->method_sto == 1)
        {
            auto ngauss = std::bind(&Sto_Func<double>::ngauss, &this->stofunc, std::placeholders::_1);
            che.calcoef_real(ngauss);
            tmpsto = BlasConnector::dot(dos_nche, che.coef_real, 1, spolyv.data(), 1);
        }
        else
        {
            auto nroot_gauss = std::bind(&Sto_Func<double>::nroot_gauss, &this->stofunc, std::placeholders::_1);
            che.calcoef_real(nroot_gauss);
            tmpsto = vTMv<double, base_device::DEVICE_CPU>(che.coef_real, spolyv.data(), dos_nche);
        }
        if (PARAM.inp.nbands > 0)
        {
            for (int ik = 0; ik < nk; ++ik)
            {
                double* en = &(this->p_elec->ekb(ik, 0));
                for (int ib = 0; ib < PARAM.inp.nbands; ++ib)
                {
                    tmpks += this->stofunc.gauss(en[ib]) * p_kv->wk[ik] / 2;
                }
            }
        }
        tmpks /= GlobalV::NPROC_IN_POOL;

        double tmperror = 0;
        if (this->method_sto == 1)
        {
            tmperror = che.coef_real[dos_nche - 1] * spolyv[dos_nche - 1];
        }
        else
        {
            const int norder = dos_nche;
            double last_coef = che.coef_real[norder - 1];
            double last_spolyv = spolyv[norder * norder - 1];
            tmperror = last_coef
                       * (BlasConnector::dot(norder, che.coef_real, 1, spolyv.data() + norder * (norder - 1), 1)
                          + BlasConnector::dot(norder, che.coef_real, 1, spolyv.data() + norder - 1, norder)
                          - last_coef * last_spolyv);
        }

        if (ie % n10 == n10 - 1)
        {
            std::cout << percent << "%"
                      << " ";
            percent += 10;
        }
        sto_dos[ie] = tmpsto;
        ks_dos[ie] = tmpks;
        error[ie] = tmperror;
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, ks_dos.data(), ndos, MPI_DOUBLE, MPI_SUM, INT_BGROUP);
    Parallel_Reduce::reduce_all(sto_dos.data(), ndos);
    Parallel_Reduce::reduce_all(error.data(), ndos);
    #endif
    if (GlobalV::MY_RANK == 0)
    {
        std::string dosfile = PARAM.globalv.global_out_dir + "dos_sdft.txt";
        ofsdos.open(dosfile.c_str());
        double maxerror = 0;
        double sum = 0;
        double min_dos = std::numeric_limits<double>::max();

        ofsdos << ndos << " # number of points" << std::endl;
        ofsdos << "#" << std::setw(19) << "energy(eV)"
            << std::setw(20) << "dos(eV^-1)" << std::setw(20) << "sum"
            << std::setw(20) << "error(eV^-1)" << std::endl;
        for (int ie = 0; ie < ndos; ++ie)
        {
            double tmperror = 2.0 * std::abs(error[ie]);
            if (maxerror < tmperror)
            {
                maxerror = tmperror;
            }
            double dos = 2.0 * (ks_dos[ie] + sto_dos[ie]) / ModuleBase::Ry_to_eV;
            min_dos = std::min(min_dos, dos);
            sum += dos;
            ofsdos << std::setw(20) << emin + ie * de << std::setw(20) << dos << std::setw(20) << sum * de
                   << std::setw(20) << tmperror << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Finish DOS" << std::endl;
        std::cout << std::scientific << "DOS max absolute Chebyshev Error: " << maxerror << std::endl;
        ofsdos.close();

        // Result-level checks: a curve that goes visibly negative, or an output
        // energy outside the Chebyshev interval, means the printed numbers must
        // not be used without changing dos_sigma / dos_nche / the window.
        const double dE = *this->stofunc.Emax - *this->stofunc.Emin;
        const double t_lo = (2.0 * this->emin / ModuleBase::Ry_to_eV - *this->stofunc.Emax - *this->stofunc.Emin) / dE;
        const double t_hi = (2.0 * this->emax / ModuleBase::Ry_to_eV - *this->stofunc.Emax - *this->stofunc.Emin) / dE;
        const double max_abs_t = std::max(std::fabs(t_lo), std::fabs(t_hi));
        for (const std::string& msg: ModuleIO::check_dos_result(min_dos, max_abs_t))
        {
            dos_log(" WARNING: DOS: " + msg, true);
        }
    }
    ModuleBase::timer::end("Sto_DOS", "DOS Loop");
    ModuleBase::timer::end("Sto_DOS", "caldos");
    return;
}

template class Sto_DOS<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Sto_DOS<double, base_device::DEVICE_GPU>;
#endif
