
#include "module_io/output_dm1.h"
#include "module_io/write_dm_sparse.h"

namespace ModuleIO
{

Output_DM1::Output_DM1(
			int nspin, 
			int istep, 
			Local_Orbital_Charge& LOC, 
			Record_adj &RA, 
			K_Vectors& kv, 
			const elecstate::DensityMatrix<std::complex<double>,double>* DM)
	: _nspin(nspin), _istep(istep), loc(LOC), _RA(RA), _kv(kv), _DM(DM)
{
}

void Output_DM1::write(void)
{
    double** dm2d;
    dm2d = new double*[_nspin];
    for (int is = 0; is < _nspin; is++)
    {
        dm2d[is] = new double[loc.ParaV->nnr];
        ModuleBase::GlobalFunc::ZEROS(dm2d[is], loc.ParaV->nnr);
    }

    // get DMK from DM
    loc.dm_k.resize(_kv.nks);
    for (int ik = 0; ik < _kv.nks; ++ik)
    {
        loc.set_dm_k(ik,_DM->get_DMK_pointer(ik));
    }

    // cal DMR in LOC
    loc.cal_dm_R(loc.dm_k, _RA, dm2d, _kv);
    for (int is = 0; is < _nspin; is++)
    {
        write_dm1(is, _istep, dm2d, loc.ParaV, loc.DMR_sparse);
    }

    for (int is = 0; is < _nspin; is++)
    {
        delete[] dm2d[is];
    }
    delete[] dm2d;
}

} // namespace ModuleIO
