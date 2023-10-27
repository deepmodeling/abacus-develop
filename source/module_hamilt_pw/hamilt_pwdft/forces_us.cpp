#include "forces.h"
#include "module_base/libm/libm.h"
#include "module_base/math_ylmreal.h"
#include "module_base/timer.h"
#include "module_elecstate/elecstate_pw.h"

// This routine computes the contribution to atomic forces due
// to the dependence of the Q function on the atomic position.
// \[ F_{j,\text{at}} = \sum_G \sum_{lm} iG_j\ \text{exp}(-iG*R_\text{at})
//    V^*(G)\ Q_{lm}(G)\ \text{becsum}(lm,\text{at}) \]
// where:
// \[ \text{becsum}(lm,\text{at}) = \sum_i \langle \psi_i|\beta_l\rangle
//    w_i\langle \beta_m|\psi_i\rangle \]
// On output: the contribution is added to \(\text{forcenl}\).
template <typename FPTYPE, typename Device>
void Forces<FPTYPE, Device>::cal_force_us(ModuleBase::matrix& forcenl,
                                          ModulePW::PW_Basis* rho_basis,
                                          pseudopot_cell_vnl* ppcell_in,
                                          const elecstate::ElecState& elec,
                                          const UnitCell& ucell)
{
    ModuleBase::TITLE("Forces", "cal_force_us");
    ModuleBase::timer::tick("Forces", "cal_force_us");

    const int npw = rho_basis->npw;
    const int nh_tot = ppcell_in->nhm * (ppcell_in->nhm + 1) / 2;
    const std::complex<double> fac = ModuleBase::NEG_IMAG_UNIT * ucell.tpiba;
    const std::complex<double> ci_tpi = ModuleBase::IMAG_UNIT * ModuleBase::TWO_PI;
    double* becsum = static_cast<const elecstate::ElecStatePW<std::complex<FPTYPE>, Device>&>(elec).becsum;

    ModuleBase::matrix forceq(ucell.nat, 3);

    ModuleBase::matrix veff = elec.pot->get_effective_v();
    ModuleBase::ComplexMatrix vg(GlobalV::NSPIN, npw);
    // fourier transform of the total effective potential
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        rho_basis->real2recip(&veff(is, 0), &vg(is, 0));
    }

    ModuleBase::matrix ylmk0(ppcell_in->lmaxq * ppcell_in->lmaxq, npw);
    ModuleBase::YlmReal::Ylm_Real(ppcell_in->lmaxq * ppcell_in->lmaxq, npw, rho_basis->gcar, ylmk0);

    double* qnorm = new double[npw];
    for (int ig = 0; ig < npw; ig++)
    {
        qnorm[ig] = rho_basis->gcar[ig].norm() * ucell.tpiba;
    }

    for (int it = 0; it < ucell.ntype; it++)
    {
        Atom* atom = &ucell.atoms[it];
        if (atom->ncpp.tvanp)
        {
            // nij = max number of (ih,jh) pairs per atom type nt
            // qgm contains the Q functions in G space
            const int nij = atom->ncpp.nh * (atom->ncpp.nh + 1) / 2;
            ModuleBase::ComplexMatrix qgm(nij, npw);

            // Compute and store Q(G) for this atomic species
            // (without structure factor)
            int ijh = 0;
            for (int ih = 0; ih < atom->ncpp.nh; ih++)
            {
                for (int jh = ih; jh < atom->ncpp.nh; jh++)
                {
                    ppcell_in->radial_fft_q(npw, ih, jh, it, qnorm, ylmk0, &qgm(ijh, 0));
                    ijh++;
                }
            }
            double* qgm_data = reinterpret_cast<double*>(qgm.c);

            ModuleBase::ComplexArray aux1(3, atom->na, npw);
            ModuleBase::realArray ddeeq(GlobalV::NSPIN, 3, atom->na, nij);
            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                for (int ia = 0; ia < atom->na; ia++)
                {
                    // aux1 = product of potential, structure factor and iG
                    for (int ig = 0; ig < npw; ig++)
                    {
                        double arg = rho_basis->gcar[ig] * atom->tau[ia];
                        std::complex<double> cfac = fac * vg(is, ig) * ModuleBase::libm::exp(ci_tpi * arg);
                        for (int ipol = 0; ipol < 3; ipol++)
                        {
                            aux1(ipol, ia, ig) = cfac * rho_basis->gcar[ig][ipol];
                        }
                    }
                }
                double* aux1_data = reinterpret_cast<double*>(aux1.ptr);

                // ddeeq = dot product of aux1 with the Q functions
                // No need for special treatment of the G=0 term (is zero)
                const char transa = 'C';
                const char transb = 'N';
                const int dim = 2 * npw;
                const double zero = 0;
                for (int ipol = 0; ipol < 3; ipol++)
                {
                    dgemm_(&transa,
                           &transb,
                           &nij,
                           &atom->na,
                           &dim,
                           &(ucell.omega),
                           qgm_data,
                           &dim,
                           &aux1_data[ipol * dim * atom->na],
                           &dim,
                           &zero,
                           &ddeeq(is, ipol, 0, 0),
                           &nij);
                }
            }

            for (int is = 0; is < GlobalV::NSPIN; is++)
            {
                for (int ia = 0; ia < atom->na; ia++)
                {
                    const int iat = ucell.itia2iat(it, ia);
                    const int index = is * ucell.nat * nh_tot + iat * nh_tot;
                    for (int ipol = 0; ipol < 3; ipol++)
                    {
                        for (int ijh = 0; ijh < nij; ijh++)
                        {
                            forceq(iat, ipol) += ddeeq(is, ipol, ia, ijh) * becsum[index + ijh];
                        }
                    }
                }
            }
        }
    }

    Parallel_Reduce::reduce_all(forceq.c, forceq.nr * forceq.nc);
    forcenl += forceq;

    delete[] qnorm;

    ModuleBase::timer::tick("Forces", "cal_force_us");
}

template class Forces<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Forces<double, psi::DEVICE_GPU>;
#endif