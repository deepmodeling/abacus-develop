#include <complex>
#include <string>
// Description: Getters for elecstate module
namespace elecstate
{

// get the value of GlobalC::rhopw->nrxx
const int get_rhopw_nrxx();
// get the value of GlobalC::rhopw->nxyz
const int get_rhopw_nxyz();
// get the value of GlobalC::ucell.omega
const double get_ucell_omega();
// get the value of H_Ewald_pw::ewald_energy
const double get_ewald_energy();
// get the value of H_Hartree_pw::hartree_energy
const double get_hartree_energy();
// get the value of Efield::etotefield
const double get_etot_efield();
// get the value of Gatefield::etotgatefield
const double get_etot_gatefield();
// get the value of XC_Functional::get_func_type()
const int get_xc_functional_type();
// get the value of INPUT::vdw_method
const std::string get_input_vdw_method();
#ifdef __LCAO
// get dftu energy
const double get_dftu_energy();
#endif
#ifdef __DEEPKS
// get lcao deepks E_delta
const double get_lcao_deepks_E_delta();
// get lcao deepks e_delta_band
const double get_lcao_deepks_e_delta_band();
#endif
// get solvent model Ael
const double get_solvent_model_Ael();
// get solvent model Acav
const double get_solvent_model_Acav();
const double get_tot_magnetization();
const double get_abs_magnetization();
const double get_tot_magnetization_nc_x();
const double get_tot_magnetization_nc_y();
const double get_tot_magnetization_nc_z();
#ifdef __EXX
#ifdef __LCAO
// get exx_lip exx_energy
const double get_exx_lip_exx_energy();
// get exx_info.info_ri.real_number
const bool get_exx_info_ri_real_number();
// get exx_lri_double.Eexx
const double get_exx_lri_double_Eexx();
// get exx_lri_complex.Eexx
const std::complex<double> get_exx_lri_complex_Eexx();
// get exx_info.info_global.cal_exx
const bool get_exx_info_global_cal_exx();
// get exx_info.info_global.hybrid_alpha
const double get_exx_info_global_hybrid_alpha();
#endif // __LCAO
#endif // __EXX

} // namespace elecstate
