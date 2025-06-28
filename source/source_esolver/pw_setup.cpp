#include "pw_setup.h"

// To get POOL_WORLD
#include "source_base/parallel_comm.h"
// Print information 
#include "source_io/print_info.h"

void ModuleESolver::pw_setup(const Input_para& inp,
		const UnitCell& ucell, 
		const ModulePW::PW_Basis& pw_rho,
		K_Vectors& kv,
		ModulePW::PW_Basis_K& pw_wfc)
{
    ModuleBase::TITLE("ModuleESolver", "pw_setup");

    //! new plane wave basis, fft grids, etc.
#ifdef __MPI
    pw_wfc.initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif

	pw_wfc.initgrids(inp.ref_cell_factor * ucell.lat0,
			ucell.latvec,
			pw_rho.nx,
			pw_rho.ny,
			pw_rho.nz);

    pw_wfc.initparameters(false, inp.ecutwfc, kv.get_nks(), kv.kvec_d.data());

#ifdef __MPI
    if (inp.pw_seed > 0)
    {
        MPI_Allreduce(MPI_IN_PLACE, &pw_wfc.ggecut, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    // qianrui add 2021-8-13 to make different kpar parameters can get the same
    // results
#endif

    pw_wfc.fft_bundle.initfftmode(inp.fft_mode);
    pw_wfc.setuptransform();

    //! initialize the number of plane waves for each k point
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        kv.ngk[ik] = pw_wfc.npwk[ik];
    }

    pw_wfc.collect_local_pw(inp.erf_ecut, inp.erf_height, inp.erf_sigma);

    ModuleIO::print_wfcfft(inp, pw_wfc, GlobalV::ofs_running);

}

