#include "stress_func.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_parameter/parameter.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/fs_kin_tools.h"

//calculate the kinetic stress in PW base
template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_kin(ModuleBase::matrix& sigma,
                                             const ModuleBase::matrix& wg,
                                             ModuleSymmetry::Symmetry* p_symm,
                                             K_Vectors* p_kv,
                                             ModulePW::PW_Basis_K* wfc_basis,
                                             const UnitCell& ucell_in,
                                             const psi::Psi<complex<FPTYPE>, Device>* psi_in)
{
    ModuleBase::TITLE("Stress_Func","stress_kin");
	ModuleBase::timer::tick("Stress_Func","stress_kin");
	this->ucell = &ucell_in;
	hamilt::FS_Kin_tools<FPTYPE, Device> kin_tool(*this->ucell, p_kv, wfc_basis, wg);
    for (int ik = 0; ik < wfc_basis->nks; ++ik)
    {
        int nbands_occ = wg.nc;
        while (wg(ik, nbands_occ - 1) == 0.0)
        {
            nbands_occ--;
            if (nbands_occ == 0)
            {
                break;
            }
        }
        const int npm = nbands_occ;
        kin_tool.cal_gk(ik);
        kin_tool.cal_stress_kin(ik, npm, true, &psi_in[0](ik, 0, 0));
    }
    kin_tool.symmetrize_stress(p_symm, sigma);
		
	ModuleBase::timer::tick("Stress_Func","stress_kin");
	return;
}

template class Stress_Func<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, base_device::DEVICE_GPU>;
#endif