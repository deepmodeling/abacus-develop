#include "module_io/output_dm.h"
#include "module_io/dm_io.h"

namespace ModuleIO
{

Output_DM::Output_DM(const Grid_Technique& GridT,
                     int is,
                     int iter,
                     const int& precision,
                     const int& out_dm,
                     double*** DM,
                     const double& ef,
                     const UnitCell* ucell,
                     const std::string& directory,
                     bool gamma_only_local)
    : _GridT(GridT),
      _is(is),
      _iter(iter),
      _precision(precision),
      _out_dm(out_dm),
      _DM(DM),
      _ef(ef),
      _ucell(ucell),
      _directory(directory),
      _gamma_only_local(gamma_only_local)
{
    if (gamma_only_local)
    {
        this->_fn = this->_directory + "/tmp_SPIN" + std::to_string(this->_is + 1) + "_DM";
    }
    else
    {
        this->_fn = this->_directory + "/tmp_SPIN" + std::to_string(this->_is + 1) + "_DM_R";
    }
}
void Output_DM::write()
{
    ModuleIO::write_dm(
#ifdef __MPI
        _GridT.trace_lo,
#endif
        _is,
        _iter,
        _fn,
        _precision,
        _out_dm,
        _DM,
        _ef,
        _ucell);
}
} // namespace ModuleIO