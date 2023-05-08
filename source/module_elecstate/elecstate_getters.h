// Description: Getters for elecstate module
namespace elecstate
{

// get the value of GlobalC::rhopw->nrxx
const int get_rhopw_nrxx();
// get the value of GlobalC::rhopw->nxyz
const int get_rhopw_nxyz();
// get the value of GlobalC::ucell.omega
const double get_ucell_omega();
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
}
