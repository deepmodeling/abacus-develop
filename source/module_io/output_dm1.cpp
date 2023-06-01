
#include "module_io/output_dm1.h"
#include "module_io/write_dm_sparse.h"

namespace ModuleIO
{

Output_DM1::Output_DM1(int nspin, int istep, Local_Orbital_Charge& LOC, Record_adj& RA, K_Vectors& kv)
    : _nspin(nspin), _istep(istep), _LOC(LOC), _RA(RA), _kv(kv)
{
}

void Output_DM1::write()
{
    double** dm2d;
    dm2d = new double*[_nspin];
    for (int is = 0; is < _nspin; is++)
    {
        dm2d[is] = new double[_LOC.ParaV->nnr];
        ModuleBase::GlobalFunc::ZEROS(dm2d[is], _LOC.ParaV->nnr);
    }
    _LOC.cal_dm_R(_LOC.dm_k, _RA, dm2d, _kv);
    for (int is = 0; is < _nspin; is++)
    {
        write_dm1(is, _istep, dm2d, _LOC.ParaV, _LOC.DMR_sparse);
    }

    for (int is = 0; is < _nspin; is++)
    {
        delete[] dm2d[is];
    }
    delete[] dm2d;
}

} // namespace ModuleIO