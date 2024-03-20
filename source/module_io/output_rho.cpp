#include "module_io/output_rho.h"

#include "module_io/rho_io.h"

namespace ModuleIO
{
Output_Rho::Output_Rho(const ModulePW::PW_Basis_Big* pw_big,
                       const ModulePW::PW_Basis* pw_rho,
                       int is,
                       int nspin,
                       const double* data,
                       int iter,
                       double ef,
                       const UnitCell* ucell,
                       const std::string directory,
                       int precision,
                       const std::string tag,
                       const std::string prefix)
    : _pw_big(pw_big),
      _pw_rho(pw_rho),
      _is(is),
      _nspin(nspin),
      _data(data),
      _iter(iter),
      _ef(ef),
      _ucell(ucell),
      _directory(directory),
      _precision(precision),
      _tag(tag),
      _prefix(prefix)
{
    if (prefix != "None")
    {
        this->_fn = this->_directory + "/" + _prefix + "SPIN" + std::to_string(this->_is + 1) + "_" + _tag + ".cube";
    }
    else
    {
        this->_fn = this->_directory + "/SPIN" + std::to_string(this->_is + 1) + "_" + _tag + ".cube";
    }
}
void Output_Rho::write()
{
    write_rho(
#ifdef __MPI
        _pw_big->bz,
        _pw_big->nbz,
        _pw_rho->nplane,
        _pw_rho->startz_current,
#endif
        _data,
        _is,
        _nspin,
        _iter,
        _fn,
        _pw_rho->nx,
        _pw_rho->ny,
        _pw_rho->nz,
        _ef,
        _ucell,
        _precision);
}
} // namespace ModuleIO