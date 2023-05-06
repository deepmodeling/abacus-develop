// Date: 2021/4/10
// Description: Getters for elecstate module
namespace elecstate
{

// get the value of GlobalC::rhopw->nrxx, need to be refactored in the future
const int get_rhopw_nrxx();
// get the value of GlobalC::rhopw->nxyz, need to be refactored in the future
const int get_rhopw_nxyz();
// get the value of GlobalC::ucell.omega, need to be refactored in the future
const double get_ucell_omega();

}