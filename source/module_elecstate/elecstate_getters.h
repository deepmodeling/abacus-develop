#ifndef ELECSTATE_GETTERS_H
#define ELECSTATE_GETTERS_H

#include <string>

// Description: Getters for elecstate module
namespace elecstate
{

/// @brief get the value of GlobalC::ucell.omega
double get_ucell_omega();
/// @brief get the value of GlobalC::ucell.tpiba
double get_ucell_tpiba();
/// @brief get the value of XC_Functional::func_type
int get_xc_func_type();
/// @brief get the value of INPUT.vdw_method
std::string get_input_vdw_method();
/// @brief get the value of GlobalC::ucell.magnet.tot_magnetization
double get_ucell_tot_magnetization();
/// @brief get the value of GlobalC::ucell.magnet.abs_magnetization
double get_ucell_abs_magnetization();
/// @brief get the value of GlobalC::ucell.magnet.tot_magnetization_nc[0]
double get_ucell_tot_magnetization_nc_x();
/// @brief get the value of GlobalC::ucell.magnet.tot_magnetization_nc[1]
double get_ucell_tot_magnetization_nc_y();
/// @brief get the value of GlobalC::ucell.magnet.tot_magnetization_nc[2]
double get_ucell_tot_magnetization_nc_z();
/// @brief get the type of KS_SOLVER
std::string get_ks_solver_type();

} // namespace elecstate

#endif // ELECSTATE_GETTERS_H
