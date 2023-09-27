#include "module_io/dos_nao.h"
#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_io/input.h"

namespace ModuleIO
{
    /// @brief manege the output of dos in numerical atomic basis case
/// @param[in] psi
/// @param[in] uhm
/// @param[in] ekb
/// @param[in] wg
/// @param[in] dos_edelta_ev
/// @param[in] dos_scale
/// @param[in] dos_sigma
/// @param[in] kv
/// @param[in] Pkpoints
/// @param[in] ucell
/// @param[in] eferm
/// @param[in] nbands
    template<typename T>
    void out_dos_nao(
        const psi::Psi<T>* psi,
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
        int nbands,
        hamilt::Hamilt<T>* p_ham)
{
    ModuleBase::TITLE("Driver", "init");
    write_dos_lcao(psi, uhm, ekb, wg, dos_edelta_ev, dos_scale, dos_sigma, kv, p_ham);

    int nspin0 = (GlobalV::NSPIN == 2) ? 2 : 1;
    if (INPUT.out_dos == 3)
    {
        for (int i = 0; i < nspin0; i++)
        {
            std::stringstream ss3;
            ss3 << GlobalV::global_out_dir << "Fermi_Surface_" << i << ".bxsf";
            nscf_fermi_surface(ss3.str(), nbands, eferm.ef, kv, Pkpoints, ucell, ekb);
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

    template void out_dos_nao(const psi::Psi<double>* psi,
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
        int nbands,
        hamilt::Hamilt<double>* p_ham);
    template void out_dos_nao(const psi::Psi<std::complex<double>>* psi,
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
        int nbands,
        hamilt::Hamilt<std::complex<double>>* p_ham);
} // namespace ModuleIO