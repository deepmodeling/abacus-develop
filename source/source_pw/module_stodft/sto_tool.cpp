#include "sto_tool.h"

#include <sstream>

#include "source_base/math_chebyshev.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_device.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include <algorithm>
#include <cmath>
#include <vector>

template <typename FPTYPE, typename Device>
void check_che_op<FPTYPE, Device>::operator()(const int& nche_in,
                                              const double& try_emin,
                                              const double& try_emax,
                                              const int& nbands_sto,
                                              K_Vectors* p_kv,
                                              Stochastic_WF<std::complex<FPTYPE>, Device>* p_stowf,
                                              StoHamiltPW<std::complex<FPTYPE>, Device>* p_hamilt_sto)
{
    //------------------------------
    //      Convergence test
    //------------------------------
    bool change = false;
    const int nk = p_kv->get_nks();
    ModuleBase::Chebyshev<FPTYPE, Device> chetest(nche_in);
    int ntest0 = 5;
    *p_hamilt_sto->emax = try_emax;
    *p_hamilt_sto->emin = try_emin;
    // if (PARAM.inp.nbands > 0)
    // {
    //     double tmpemin = 1e10;
    //     for (int ik = 0; ik < nk; ++ik)
    //     {
    //         tmpemin = std::min(tmpemin, this->pelec->ekb(ik, PARAM.inp.nbands - 1));
    //     }
    //     *p_hamilt_sto->emin = tmpemin;
    // }
    // else
    // {
    //     *p_hamilt_sto->emin = 0;
    // }
    for (int ik = 0; ik < nk; ++ik)
    {
        p_hamilt_sto->updateHk(ik);
        const int npw = p_kv->ngk[ik];
        std::complex<FPTYPE>* pchi = nullptr;
        psi::Psi<std::complex<FPTYPE>, Device> randchi_d;
        randchi_d.resize(1, 1, npw); // holds the independent random probe of this k point
        // Use the real stochastic (or complete-basis) orbitals whenever they
        // exist.  Previously this decision was based only on "nbands_sto==0",
        // while Stochastic_WF::nchi was never initialized; as a result the DOS
        // range check always used random vectors, even for MDFT.  Random
        // vectors have a Rayleigh quotient close to the mean spectral energy
        // and cannot locate a low-lying eigenvalue, so Emin could be set far
        // above the occupied states.
        const bool use_random = (p_stowf->nchi <= 0) || (nbands_sto == 0 && p_stowf->nchip_max <= 0);
        const int ntest_real = use_random ? 0 : std::min(ntest0, p_stowf->nchip[ik]);
        // Always add one independent random probe on top of the real orbitals.
        // With `nbands_sto all` the stochastic orbitals are coordinate vectors,
        // so a handful of them can miss a diagonal state of another invariant
        // subspace; a dense random vector does not depend on that choice.
        const int ntest = use_random ? 1 : ntest_real + 1;
        for (int i = 0; i < ntest; ++i)
        {
            const bool use_this_random = use_random || i == ntest_real;
            if (use_this_random)
            {
                std::vector<std::complex<FPTYPE>> randchi(npw);   
                for (int ig = 0; ig < npw; ++ig)
                {
                    FPTYPE rr = std::rand() / FPTYPE(RAND_MAX);
                    FPTYPE arg = std::rand() / FPTYPE(RAND_MAX);
                    randchi[ig] = std::complex<FPTYPE>(rr * cos(arg), rr * sin(arg));
                }
                syncmem_complex_h2d_op()(randchi_d.get_pointer(), randchi.data(), npw);
                pchi = randchi_d.get_pointer();
            }
            else if (PARAM.inp.nbands > 0)
            {
                pchi = &p_stowf->chiortho[0](ik, i, 0);
            }
            else
            {
                pchi = &p_stowf->chi0[0](ik, i, 0);
            }
            while (true)
            {
                bool converge;
                auto hchi_norm = std::bind(&StoHamiltPW<std::complex<FPTYPE>, Device>::hPsi_norm,
                                           p_hamilt_sto,
                                           std::placeholders::_1,
                                           std::placeholders::_2,
                                           std::placeholders::_3);
                converge = chetest.checkconverge(hchi_norm,
                                                 pchi,
                                                 npw,
                                                 p_stowf->npwx,
                                                 *p_hamilt_sto->emax,
                                                 *p_hamilt_sto->emin,
                                                 2.0);

                if (!converge)
                {
                    change = true;
                }
                else
                {
                    break;
                }
            }
        }

        if (ik == nk - 1)
        {
#ifdef __MPI
            Parallel_Reduce::reduce_max(*p_hamilt_sto->emax);
            Parallel_Reduce::reduce_min(*p_hamilt_sto->emin);
#endif
            // checkconverge() starts from the user's (0,0) defaults and can
            // leave the automatic lower bound positive when the first trial
            // Rayleigh quotient is positive.  A positive Emin excludes every
            // bound state and is the root of the high-temperature SCF and DOS
            // failures.  Only correct that pathological outcome here so that
            // cases whose final bound is already negative keep the original
            // numerical result bit for bit.
            // Only an automatically determined bound may be repaired: a positive
            // emin_sto given by the user (or a genuinely positive spectrum) is a
            // valid interval and must be preserved.
            const bool emin_is_automatic = (PARAM.inp.emin_sto == 0.0);
            if (emin_is_automatic && (!std::isfinite(*p_hamilt_sto->emin) || *p_hamilt_sto->emin > 0))
            {
                const double positive_emin = *p_hamilt_sto->emin;
                const double emax_abs = std::abs(*p_hamilt_sto->emax);
                const double fallback_span = std::max(5.0, 0.01 * emax_abs);
                const double safe_emin = -fallback_span;
                std::ostringstream range_msg;
                range_msg << " WARNING: positive automatic Chebyshev Emin " << positive_emin << " Ry is reset to "
                          << safe_emin << " Ry";
                *p_hamilt_sto->emin = safe_emin;
                // The pathological path that produced a positive Emin also
                // tends to underestimate Emax.  The highest plane-wave
                // eigenvalue is bounded by roughly ecutwfc; use it as a floor.
                const double ecutwfc = PARAM.inp.ecutwfc;
                const double ecut_floor = 1.05 * ecutwfc;
                if (ecutwfc > 0 && *p_hamilt_sto->emax < ecut_floor)
                {
                    range_msg << "; Chebyshev Emax " << *p_hamilt_sto->emax
                              << " Ry is raised to the ecutwfc-based floor " << ecut_floor << " Ry";
                    *p_hamilt_sto->emax = ecut_floor;
                }
                GlobalV::ofs_running << range_msg.str() << std::endl;
            }
            GlobalV::ofs_running << "New Emax " << *p_hamilt_sto->emax << " Ry; new Emin " << *p_hamilt_sto->emin
                                 << " Ry" << std::endl;
            change = false;
        }
    }
}

template <typename FPTYPE, typename Device>
psi::Psi<std::complex<FPTYPE>, Device>* gatherchi_op<FPTYPE, Device>::operator()(
    psi::Psi<std::complex<FPTYPE>, Device>& chi,
    psi::Psi<std::complex<FPTYPE>, Device>& chi_all,
    const int& npwx,
    int* nrecv_sto,
    int* displs_sto,
    const int perbands_sto)
{
    psi::Psi<std::complex<FPTYPE>, Device>* p_chi;
    p_chi = &chi;
#ifdef __MPI
    if (PARAM.inp.bndpar > 1)
    {
        p_chi = &chi_all;
        ModuleBase::timer::start("sKG", "bands_gather");
        Parallel_Common::gatherv_dev<std::complex<FPTYPE>, Device>(chi.get_pointer(),
                                                                   perbands_sto * npwx,
                                                                   chi_all.get_pointer(),
                                                                   nrecv_sto,
                                                                   displs_sto,
                                                                   BP_WORLD);
        ModuleBase::timer::end("sKG", "bands_gather");
    }
#endif
    return p_chi;
}

template struct check_che_op<double, base_device::DEVICE_CPU>;
#ifdef __FLOAT_FFTW
template struct check_che_op<float, base_device::DEVICE_CPU>;
#endif
template struct gatherchi_op<double, base_device::DEVICE_CPU>;
template struct gatherchi_op<float, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template struct check_che_op<double, base_device::DEVICE_GPU>;
#ifdef __FLOAT_FFTW
template struct check_che_op<float, base_device::DEVICE_GPU>;
#endif
template struct gatherchi_op<double, base_device::DEVICE_GPU>;
template struct gatherchi_op<float, base_device::DEVICE_GPU>;
#endif