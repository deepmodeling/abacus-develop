#include "module_io/dos_nao.h"
#include "module_io/write_dos_lcao.h"
#include "module_base/global_variable.h"
#include "module_base/tool_title.h"
#include "module_parameter/parameter.h"

namespace ModuleIO
{

template<typename T>
void out_dos_nao(
			const psi::Psi<T>* psi,
			hamilt::Hamilt<T>* p_ham,
			const Parallel_Orbitals &pv,
			const UnitCell& ucell,
			const K_Vectors& kv,
			const int nbands,
			const elecstate::efermi& eferm,
			const ModuleBase::matrix& ekb,
			const ModuleBase::matrix& wg,
			const double& dos_edelta_ev,
			const double& dos_scale,
			const double& dos_sigma)
{
    ModuleBase::TITLE("Module_IO", "out_dos_nao");

	ModuleIO::write_dos_lcao(
			psi, 
            p_ham,
			pv, 
			ucell, 
			kv, 
			PARAM.inp.nbands,
			eferm,
			ekb, 
			wg, 
			dos_edelta_ev, 
			dos_scale, 
			dos_sigma, 
			GlobalV::ofs_running);

    const int nspin0 = (PARAM.inp.nspin == 2) ? 2 : 1;
    if (PARAM.inp.out_dos == 3)
    {
        for (int i = 0; i < nspin0; i++)
        {
            std::stringstream ss3;
            ss3 << PARAM.globalv.global_out_dir << "Fermi_Surface_" << i << ".bxsf";
            nscf_fermi_surface(ss3.str(), nbands, eferm.ef, kv, ucell, ekb);
        }
    }

}


template void out_dos_nao(
			const psi::Psi<double>* psi,
			hamilt::Hamilt<double>* p_ham,
			const Parallel_Orbitals &pv,
			const UnitCell& ucell,
			const K_Vectors& kv,
			const int nbands,
			const elecstate::efermi& eferm,
			const ModuleBase::matrix& ekb,
			const ModuleBase::matrix& wg,
			const double& dos_edelta_ev,
			const double& dos_scale,
			const double& dos_sigma);

template void out_dos_nao(
			const psi::Psi<std::complex<double>>* psi,
			hamilt::Hamilt<std::complex<double>>* p_ham,
			const Parallel_Orbitals &pv,
			const UnitCell& ucell,
			const K_Vectors& kv,
			const int nbands,
			const elecstate::efermi& eferm,
			const ModuleBase::matrix& ekb,
			const ModuleBase::matrix& wg,
			const double& dos_edelta_ev,
			const double& dos_scale,
			const double& dos_sigma);

} // namespace ModuleIO
