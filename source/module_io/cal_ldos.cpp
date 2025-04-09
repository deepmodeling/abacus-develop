#include "cal_ldos.h"

#include "cube_io.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_lcao/module_gint/temp_gint/gint_interface.h"

#include <type_traits>

namespace ModuleIO
{
template <typename T>
void Cal_ldos<T>::cal_ldos_pw(const elecstate::ElecStatePW<std::complex<double>>* pelec,
                              const psi::Psi<std::complex<double>>& psi,
                              const Parallel_Grid& pgrid,
                              const UnitCell& ucell)
{
    for (int ie = 0; ie < PARAM.inp.stm_bias[2]; ie++)
    {
        // energy range for ldos (efermi as reference)
        const double en = PARAM.inp.stm_bias[0] + ie * PARAM.inp.stm_bias[1];
        const double emin = en < 0 ? en : 0;
        const double emax = en > 0 ? en : 0;

        std::vector<double> ldos(pelec->charge->nrxx);
        std::vector<std::complex<double>> wfcr(pelec->basis->nrxx);

        for (int ik = 0; ik < pelec->klist->get_nks(); ++ik)
        {
            psi.fix_k(ik);
            const double efermi = pelec->eferm.get_efval(pelec->klist->isk[ik]);
            int nbands = psi.get_nbands();

            for (int ib = 0; ib < nbands; ib++)
            {
                pelec->basis->recip2real(&psi(ib, 0), wfcr.data(), ik);
                const double eigenval = (pelec->ekb(ik, ib) - efermi) * ModuleBase::Ry_to_eV;
                if (eigenval >= emin && eigenval <= emax)
                {
                    for (int ir = 0; ir < pelec->basis->nrxx; ir++)
                    {
                        ldos[ir] += pelec->klist->wk[ik] * norm(wfcr[ir]);
                    }
                }
            }
        }

        std::stringstream fn;
        fn << PARAM.globalv.global_out_dir << "LDOS_" << en << "eV"
           << ".cube";

        const int precision = PARAM.inp.out_ldos[1];
        ModuleIO::write_vdata_palgrid(pgrid, ldos.data(), 0, PARAM.inp.nspin, 0, fn.str(), 0, &ucell, precision, 0);
    }
}

#ifdef __LCAO
template <typename T>
void Cal_ldos<T>::cal_ldos_lcao(const elecstate::ElecStateLCAO<T>* pelec,
                                const psi::Psi<T>& psi,
                                const Parallel_Grid& pgrid,
                                const UnitCell& ucell)
{
    for (int ie = 0; ie < PARAM.inp.stm_bias[2]; ie++)
    {
        // energy range for ldos (efermi as reference)
        const double en = PARAM.inp.stm_bias[0] + ie * PARAM.inp.stm_bias[1];
        const double emin = en < 0 ? en : 0;
        const double emax = en > 0 ? en : 0;

        // calculate weight (for bands not in the range, weight is zero)
        ModuleBase::matrix weight(pelec->ekb.nr, pelec->ekb.nc);
        for (int ik = 0; ik < pelec->ekb.nr; ++ik)
        {
            const double efermi = pelec->eferm.get_efval(pelec->klist->isk[ik]);

            for (int ib = 0; ib < pelec->ekb.nc; ib++)
            {
                const double eigenval = (pelec->ekb(ik, ib) - efermi) * ModuleBase::Ry_to_eV;
                if (eigenval >= emin && eigenval <= emax)
                {
                    weight(ik, ib) = pelec->klist->wk[ik];
                }
            }
        }

        // calculate dm-like for ldos
        const int nspin_dm = PARAM.inp.nspin == 2 ? 2 : 1;
        elecstate::DensityMatrix<T, double> dm_ldos(pelec->DM->get_paraV_pointer(),
                                                    nspin_dm,
                                                    pelec->klist->kvec_d,
                                                    pelec->klist->get_nks() / nspin_dm);

        elecstate::cal_dm_psi(pelec->DM->get_paraV_pointer(), weight, psi, dm_ldos);
        dm_ldos.init_DMR(*(pelec->DM->get_DMR_pointer(1)));
        dm_ldos.cal_DMR();

        // allocate ldos space
        std::vector<double> ldos_space(PARAM.inp.nspin * pelec->charge->nrxx);
        double** ldos = new double*[PARAM.inp.nspin];
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            ldos[is] = &ldos_space[is * pelec->charge->nrxx];
        }

    // calculate ldos
#ifndef __NEW_GINT
        ModuleBase::WARNING_QUIT("Cal_ldos::dm2ldos",
                                 "do not support old grid integral, please recompile with __NEW_GINT");
#else
        ModuleGint::cal_gint_rho(dm_ldos.get_DMR_vector(), PARAM.inp.nspin, ldos);
#endif

        // I'm not sure whether ldos should be output for each spin or not
        // ldos[0] += ldos[1] for nspin_dm == 2
        if (nspin_dm == 2)
        {
            BlasConnector::axpy(pelec->charge->nrxx, 1.0, ldos[1], 1, ldos[0], 1);
        }

        // write ldos to cube file
        std::stringstream fn;
        fn << PARAM.globalv.global_out_dir << "LDOS_" << en << "eV"
           << ".cube";

        const int precision = PARAM.inp.out_ldos[1];
        ModuleIO::write_vdata_palgrid(pgrid,
                                      ldos_space.data(),
                                      0,
                                      PARAM.inp.nspin,
                                      0,
                                      fn.str(),
                                      0,
                                      &ucell,
                                      precision,
                                      0);

        // free memory
        delete[] ldos;
    }
}

#endif

template class Cal_ldos<double>;               // Gamma_only case
template class Cal_ldos<std::complex<double>>; // multi-k case
} // namespace ModuleIO 
