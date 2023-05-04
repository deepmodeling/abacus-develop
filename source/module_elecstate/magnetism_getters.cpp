#include "magnetism.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

int Magnetism::get_rhopw_nrxx() const
{
    return GlobalC::rhopw->nrxx;
}
int Magnetism::get_rhopw_nxyz() const
{
    return GlobalC::rhopw->nxyz;
}
double Magnetism::get_ucell_omega() const
{
    return GlobalC::ucell.omega;
}
