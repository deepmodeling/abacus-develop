#include "cal_ldos.h"

#include "cube_io.h"

namespace ModuleIO
{
void cal_ldos(const elecstate::ElecStatePW<std::complex<double>>* pelec,
              const psi::Psi<std::complex<double>>& psi,
              const Parallel_Grid& pgrid,
              const UnitCell& ucell)
{
    // energy range for ldos (efermi as reference)
    const double emin = PARAM.inp.stm_bias < 0 ? PARAM.inp.stm_bias : 0;
    const double emax = PARAM.inp.stm_bias > 0 ? PARAM.inp.stm_bias : 0;

    std::vector<double> ldos(pelec->charge->nrxx);
    std::vector<std::complex<double>> wfcr(pelec->basis->nrxx);

    for (int ik = 0; ik < pelec->klist->get_nks(); ++ik)
    {
        psi.fix_k(ik);
        double efermi = pelec->eferm.get_efval(pelec->klist->isk[ik]);
        int nbands = psi.get_nbands();

        for (int ib = 0; ib < nbands; ib++)
        {
            pelec->basis->recip2real(&psi(ib, 0), wfcr.data(), ik);
            double eigenval = (pelec->ekb(ik, ib) - efermi) * ModuleBase::Ry_to_eV;
            if (eigenval >= emin && eigenval <= emax)
            {
                for (int ir = 0; ir < pelec->basis->nrxx; ir++)
                    ldos[ir] += pelec->klist->wk[ik] * norm(wfcr[ir]);
            }
        }
    }

    std::stringstream fn;
    fn << PARAM.globalv.global_out_dir << "LDOS_" << PARAM.inp.stm_bias << "eV"
       << ".cube";

    ModuleIO::write_vdata_palgrid(pgrid, ldos.data(), 0, PARAM.inp.nspin, 0, fn.str(), 0, &ucell, 11, 0);
}

#ifdef __LCAO
// lcao multi-k case
// void cal_ldos(elecstate::ElecState* pelec, const psi::Psi<std::complex<double>>& psi, std::vector<double>& ldos)
// {
// }

// // lcao Gamma_only case
// void cal_ldos(elecstate::ElecState* pelec, const psi::Psi<double>& psi, std::vector<double>& ldos)
// {
// }
#endif
} // namespace elecstate
