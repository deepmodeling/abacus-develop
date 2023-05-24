#include "module_io/dos_lcao.h"
#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_io/input.h"

namespace ModuleIO
{
void out_dos_lcao(const psi::Psi<double>* psid,
                  const psi::Psi<std::complex<double>>* psi,
                  LCAO_Hamilt& uhm,
                  const ModuleBase::matrix& ekb,
                  const ModuleBase::matrix& wg,
                  const double& dos_edelta_ev,
                  const double& dos_scale,
                  const double& dos_sigma,
                  const K_Vectors& kv,
                  const Parallel_Kpoints& Pkpoints,
                  const UnitCell& ucell,
                  const elecstate::efermi& eferm,
                  int nks,
                  int nbands)
{
    ModuleBase::TITLE("Driver", "init");
    write_dos_lcao(psid, psi, uhm, ekb, wg, dos_edelta_ev, dos_scale, dos_sigma, kv);

    int nspin0 = (GlobalV::NSPIN == 2) ? 2 : 1;
    if (INPUT.out_dos == 3)
    {
        for (int i = 0; i < nspin0; i++)
        {
            std::stringstream ss3;
            ss3 << GlobalV::global_out_dir << "Fermi_Surface_" << i << ".bxsf";
            nscf_fermi_surface(ss3.str(), nks, nbands, eferm.ef, kv, Pkpoints, ucell, ekb);
        }
    }

    if (nspin0 == 1)
    {
        GlobalV::ofs_running << " Fermi energy is " << eferm.ef << " Rydberg" << std::endl;
    }
    else if (nspin0 == 2)
    {
        GlobalV::ofs_running << " Fermi energy (spin = 1) is " << eferm.ef_up << " Rydberg" << std::endl;
        GlobalV::ofs_running << " Fermi energy (spin = 2) is " << eferm.ef_dw << " Rydberg" << std::endl;
    }
}
} // namespace ModuleIO