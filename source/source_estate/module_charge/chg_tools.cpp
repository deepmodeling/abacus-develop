#include "chg_tools.h"

#include <algorithm>
#include <functional>

#include "source_base/complexmatrix.h"
#include "source_base/constants.h"
#include "source_base/math_integral.h"
#include "source_base/math_sphbes.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_threading.h"
#include "source_base/tool_title.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"

#include <cassert>
#include <cmath>
#include <new>

namespace module_charge
{

double sum_rho(double* const* rho,
               const int nspin0,
               const int nrxx,
               const double omega,
               const int nxyz)
{
    ModuleBase::TITLE("module_charge", "sum_rho");

    double sum_rho = 0.0;

    for (int is = 0; is < nspin0; is++)
    {
        for (int ir = 0; ir < nrxx; ir++)
        {
            sum_rho += rho[is][ir];
        }
    }

    // multiply the sum of charge density by a factor
    sum_rho *= omega / static_cast<double>(nxyz);

#ifdef __MPI
    Parallel_Reduce::reduce_pool(sum_rho);
#endif

    // mohan fixed bug 2010-01-18,
    // sum_rho may be smaller than 1, like Na bcc.
    if (sum_rho <= 0.1)
    {
        ModuleBase::WARNING_QUIT("module_charge::sum_rho", "Can't find even an electron!");
    }

    return sum_rho;
}

double cal_rho2ne(const double* rho_in,
                  const int nrxx,
                  const double omega,
                  const int nxyz)
{
    assert(nxyz > 0); // mohan add 2025-12-02

    double ne = 0.0;
    for (int ir = 0; ir < nrxx; ir++)
    {
        ne += rho_in[ir];
    }
#ifdef __MPI
    Parallel_Reduce::reduce_pool(ne);
#endif
    ne = ne * omega / static_cast<double>(nxyz);

    return ne;
}

void check_rho(double* const* rho,
               const int nspin,
               const int nrxx,
               const double omega,
               const int nxyz,
               const double nelec)
{
    ModuleBase::TITLE("module_charge", "check_rho");

    assert(rho != nullptr);
    assert(nelec > 0.0);

    if (nspin == 1 || nspin == 4)
    {
        const double ne = cal_rho2ne(rho[0], nrxx, omega, nxyz);
        if (std::abs(ne - nelec) > 1.0e-6)
        {
            ModuleBase::WARNING("Charge", "Charge is not equal to the number of electrons!");
        }
    }
    else if (nspin == 2)
    {
        // for spin up
        const double ne_up = cal_rho2ne(rho[0], nrxx, omega, nxyz);
        if (ne_up < 0.0)
        {
            ModuleBase::WARNING_QUIT("Charge",
                "Number of spin-down electrons set in starting magnetization exceeds all available.");
        }
        // for spin down
        const double ne_dn = cal_rho2ne(rho[1], nrxx, omega, nxyz);
        if (ne_dn < 0.0)
        {
            ModuleBase::WARNING_QUIT("Charge",
                "Number of spin-up electrons set in starting magnetization exceeds all available.");
        }
        // for total charge
        if (std::abs(ne_up + ne_dn - nelec) > 1.0e-6)
        {
            ModuleBase::WARNING("Charge", "Charge is not equal to the number of electrons!");
        }
    }
}

void non_linear_core_correction(const NlcCtx& ctx,
                                double* rhocg)
{
    ModuleBase::TITLE("module_charge", "drhoc");

    const bool numeric = ctx.numeric;
    const double omega = ctx.omega;
    const double tpiba2 = ctx.tpiba2;
    const int mesh = ctx.mesh;
    const double* r = ctx.r;
    const double* rab = ctx.rab;
    const double* rhoc = ctx.rhoc;
    const double* gg_uniq = ctx.gg_uniq;
    const int ngg = ctx.ngg;

    // use labmda instead of repeating codes
    const std::function<void(int, int)> kernel = [&](int num_threads, int thread_id)
    {

    double gx = 0.0;
    double rhocg1 = 0.0;
    std::vector<double> aux_vec;

    // here we compute the fourier transform is the charge in numeric form
    if (numeric)
    {
        aux_vec.resize(mesh);
        double* aux = aux_vec.data();
        // G=0 term

        int igl0 = 0;
        if (gg_uniq [0] < 1.0e-8)
        {
            // single thread term
            if (thread_id == 0)
            {
                for (int ir = 0;ir < mesh; ir++)
                {
                    aux [ir] = r [ir] * r [ir] * rhoc [ir];
                }
                ModuleBase::Integral::Simpson_Integral(mesh, aux, rab, rhocg1);
                //rhocg [1] = fpi * rhocg1 / omega;
                rhocg [0] = ModuleBase::FOUR_PI * rhocg1 / omega;//mohan modify 2008-01-19
            }
            igl0 = 1;
        }

        int igl_beg, igl_end;
        // exclude igl0
        ModuleBase::TASK_DIST_1D(num_threads, thread_id, ngg - igl0, igl_beg, igl_end);
        igl_beg += igl0;
        igl_end += igl_beg;

        // G <> 0 term
        for (int igl = igl_beg; igl < igl_end;igl++)
        {
            gx = sqrt(gg_uniq[igl] * tpiba2);
            ModuleBase::Sphbes::Spherical_Bessel(mesh, r, gx, 0, aux);
            for (int ir = 0;ir < mesh; ir++)
            {
                aux [ir] = r[ir] * r[ir] * rhoc [ir] * aux [ir];
            } //  enddo
            ModuleBase::Integral::Simpson_Integral(mesh, aux, rab, rhocg1);
            rhocg [igl] = ModuleBase::FOUR_PI * rhocg1 / omega;
        } //  enddo
    }
    else
    {
        // here the case where the charge is in analytic form,
        // check old version before 2008-12-9
    }

    }; // end kernel

    // do not use omp parallel when this function is already in parallel block
    //
    // it is called in parallel block in Forces::cal_force_cc,
    // but not in other funtcion such as Stress_Func::stress_cc.
    ModuleBase::TRY_OMP_PARALLEL(kernel);

    return;
}

// computes the core charge on the real space 3D mesh.
void set_rho_core(const UnitCell& ucell,
                  const ModuleBase::ComplexMatrix& structure_factor,
                  const bool* numeric,
                  double* rho_core,
                  std::complex<double>* rhog_core,
                  const ModulePW::PW_Basis& rhopw)
{
    ModuleBase::TITLE("module_charge", "set_rho_core");
    ModuleBase::timer::start("module_charge", "set_rho_core");

    bool bl = false;
    for (int it = 0; it < ucell.ntype; it++)
    {
        if (ucell.atoms[it].ncpp.nlcc)
        {
            bl = true;
            break;
        }
    }

    if (!bl)
    {
        std::fill(rho_core, rho_core + rhopw.nrxx, 0.0);
        ModuleBase::timer::end("module_charge", "set_rho_core");
        return;
    }

    std::vector<double> rhocg(rhopw.ngg, 0.0);

    // three dimension.
    std::vector<std::complex<double>> vg(rhopw.npw);

    for (int it = 0; it < ucell.ntype; it++)
    {
        if (ucell.atoms[it].ncpp.nlcc)
        {
            assert(numeric != nullptr);
//----------------------------------------------------------
// EXPLAIN : drhoc compute the radial fourier transform for
// each shell of g vec
//----------------------------------------------------------
            NlcCtx nlc_ctx{
                numeric[it],
                ucell.omega,
                ucell.tpiba2,
                ucell.atoms[it].ncpp.msh,
                ucell.atoms[it].ncpp.r.data(),
                ucell.atoms[it].ncpp.rab.data(),
                ucell.atoms[it].ncpp.rho_atc.data(),
                rhopw.gg_uniq,
                rhopw.ngg
            };
            non_linear_core_correction(nlc_ctx, rhocg.data());
//----------------------------------------------------------
// EXPLAIN : multiply by the structure factor and sum
//----------------------------------------------------------
            for (int ig = 0; ig < rhopw.npw; ig++)
            {
                vg[ig] += structure_factor(it, ig) * rhocg[rhopw.ig2igg[ig]];
            }
        }
    }

    // for tmp use.
    for (int ig = 0; ig < rhopw.npw; ig++)
    {
        rhog_core[ig] = vg[ig];
    }

    rhopw.recip2real(vg.data(), rho_core);

    // test on the charge and computation of the core energy
    double rhoima = 0.0;
    double rhoneg = 0.0;
    for (int ir = 0; ir < rhopw.nrxx; ir++)
    {
        rhoneg += std::min(0.0, rhopw.fft_bundle.get_auxr_data<double>()[ir].real());
        rhoima += std::abs(rhopw.fft_bundle.get_auxr_data<double>()[ir].imag());
        // NOTE: Core charge is computed in reciprocal space and brought to real
        // space by FFT. For non smooth core charges (or insufficient cut-off)
        // this may result in negative values in some grid points.
        // Up to October 1999 the core charge was forced to be positive definite.
        // This induces an error in the force, and probably stress, calculation if
        // the number of grid points where the core charge would be otherwise neg
        // is large. The error disappears for sufficiently high cut-off, but may be
        // rather large and it is better to leave the core charge as it is.
        // If you insist to have it positive definite (with the possible problems
        // mentioned above) uncomment the following lines.  SdG, Oct 15 1999
    }

#ifdef __MPI
    // mohan fix bug 2011-04-03
    Parallel_Reduce::reduce_pool(rhoneg);
    Parallel_Reduce::reduce_pool(rhoima);
#endif

    // mohan changed 2010-2-2, make this same as in atomic_rho.
    // still lack something......
    rhoneg /= rhopw.nxyz * ucell.omega;
    rhoima /= rhopw.nxyz * ucell.omega;

    // calculate core_only exch-corr energy etxcc=E_xc[rho_core] if required
    // The term was present in previous versions of the code but it shouldn't
    ModuleBase::timer::end("module_charge", "set_rho_core");
}

} // namespace module_charge
