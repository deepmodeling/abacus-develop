#include "source_pw/module_pwdft/setup_pwwfc.h" // pw_wfc
#include "source_base/parallel_comm.h" // POOL_WORLD
#include "source_base/parallel_reduce.h" // Parallel_Reduce
#include "source_io/module_output/print_info.h" // print information
#include "source_io/module_parameter/parameter.h"

void pw::teardown_pwwfc(ModulePW::PW_Basis_K* &pw_wfc)
{
    delete pw_wfc;
}

void pw::setup_pwwfc(const Input_para& inp,
        const UnitCell& ucell, 
        const ModulePW::PW_Basis& pw_rho,
        K_Vectors& kv,
        ModulePW::PW_Basis_K* &pw_wfc)
{
    ModuleBase::TITLE("pw", "pw_setup");

    std::string fft_device = inp.device;

    //! setup pw_wfc
    // currently LCAO doesn't support GPU acceleration of FFT
    if(inp.basis_type == "lcao")
    {
        fft_device = "cpu";
    }
    std::string fft_precision = inp.precision;
#ifdef __FLOAT_FFTW
    if (inp.cal_cond && inp.esolver_type == "sdft")
    {
        fft_precision = "mixing";
    }
#endif

    pw_wfc = new ModulePW::PW_Basis_K_Big(fft_device, fft_precision);


    // for LCAO calculations, we need to set bx, by, and bz
    ModulePW::PW_Basis_K_Big* tmp = static_cast<ModulePW::PW_Basis_K_Big*>(pw_wfc);
    tmp->setbxyz(inp.bx, inp.by, inp.bz);



    //! new plane wave basis, fft grids, etc.
#ifdef __MPI
    pw_wfc->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif

    // NOTE(liuyu): ref_cell_factor is currently forced to 1.0 in
    // read_input_item_md.cpp because the reference-cell mechanism is
    // disabled for both pw_rho and pw_wfc. The wfc FFT grid shares the
    // same nx/ny/nz as pw_rho, so the same staleness issue applies:
    // when ref_cell_factor > 1, pw_wfc->lat0/tpiba/G/GGT/omega hold
    // reference-cell values, which leaks into wfc IO (read_wfc_pw,
    // write_wfc_pw), cal_energies, and other paths that read these
    // members as physical-cell quantities. If re-enabled in the future,
    // see the comment in setup_pwrho.cpp for the required refactor.
    pw_wfc->initgrids(inp.ref_cell_factor * ucell.lat0,
            ucell.latvec,
            pw_rho.nx,
            pw_rho.ny,
            pw_rho.nz);

    pw_wfc->initparameters(false, inp.ecutwfc, kv.get_nks(), kv.kvec_d.data());
#ifdef __MPI
    if (inp.pw_seed > 0)
    {
        Parallel_Reduce::reduce_max( pw_wfc->ggecut);
    }
    // qianrui add 2021-8-13 to make different kpar parameters can get the same result
#endif

    pw_wfc->fft_bundle.initfftmode(inp.fft_mode);
    pw_wfc->setuptransform();

    //! initialize the number of plane waves for each k point
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        kv.ngk[ik] = pw_wfc->npwk[ik];
    }

    pw_wfc->collect_local_pw(inp.erf_ecut, inp.erf_height, inp.erf_sigma);

    ModuleIO::print_wfcfft(inp, *pw_wfc, GlobalV::ofs_running);

    return;
}

