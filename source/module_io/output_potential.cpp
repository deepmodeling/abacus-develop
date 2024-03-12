#include "output_potential.h"
#include "potential_io.h"

namespace ModuleIO
{

Output_Potential::Output_Potential(const ModulePW::PW_Basis_Big* pw_big,
                                   ModulePW::PW_Basis* pw_rho,
                                   int nspin,
                                   int iter,
                                   int out_pot,
                                   const ModuleBase::matrix& v_effective,
                                   const double* v_effective_fixed,
                                   const UnitCell* ucell,
                                   const Charge* const charge,
                                   int precision,
                                   const std::string directory,
                                   const std::string tag,
                                   const std::string prefix)
    : _pw_big(pw_big),
      _pw_rho(pw_rho),
      _nspin(nspin),
      _iter(iter),
      _out_pot(out_pot),
      _v_effective(v_effective),
      _v_effective_fixed(v_effective_fixed),
      _ucell(ucell),
      _charge(charge),
      _precision(precision),
      _directory(directory),
      _prefix(prefix),
      _tag(tag)
{
    _fn_ElecStatic = _directory + "/" + "ElecStaticPot.cube";
}

void Output_Potential::write()
{
    if (_out_pot == 1)
    {
        for (int is = 0; is < _nspin; is++)
        {
            if (_prefix != "None")
            {
                _fn_Pot = _directory + "/" + _prefix + "SPIN" + std::to_string(is + 1) + "_" + _tag + ".cube";
            }
            else
            {
                _fn_Pot = _directory + "/SPIN" + std::to_string(is + 1) + "_" + _tag + ".cube";
            }
            ModuleIO::write_potential(
#ifdef __MPI
                _pw_big->bz,
                _pw_big->nbz,
                _pw_rho->nplane,
                _pw_rho->startz_current,
#endif
                is,
                _iter,
                _fn_Pot,
                _pw_rho->nx,
                _pw_rho->ny,
                _pw_rho->nz,
                _v_effective,
                _precision);
        }
    }
    else if (_out_pot == 2)
    {
        ModuleIO::write_elecstat_pot(
#ifdef __MPI
            _pw_big->bz,
            _pw_big->nbz,
#endif
            _fn_ElecStatic,
            _pw_rho,
            _charge,
            _ucell,
            _v_effective_fixed); // output 'Hartree + local pseudopot'
    }
}

} // namespace ModuleIO