#include "source_pw/module_pwdft/setup_pwrho.h"
#include "source_io/module_output/print_info.h" // use print_rhofft
#include "source_base/parallel_comm.h" // use POOL_WORLD

void pw::setup_pwrho(
        UnitCell& ucell, // unitcell 
        const bool double_grid, // for USPP
        bool &pw_rho_flag, // flag for allocation of pw_rho
        ModulePW::PW_Basis* &pw_rho, // pw for rhod
        ModulePW::PW_Basis* &pw_rhod, // pw for rhod
        ModulePW::PW_Basis_Big* &pw_big, // pw for rhod
        const std::string &classname,
        const Input_para& inp) // input parameters *
{
    ModuleBase::TITLE("pw", "setup_pwrho");

    std::string fft_device = inp.device;
    std::string fft_precision = inp.precision;

    // LCAO basis doesn't support GPU acceleration on FFT currently
    if(inp.basis_type == "lcao")
    {
        fft_device = "cpu";
    }

    // single, double, or mixing precision calculations
    if ((inp.precision=="single") || (inp.precision=="mixing"))
    {
        fft_precision = "mixing";
    }
    else if (inp.precision=="double")
    {
        fft_precision = "double";
    }

    // for GPU
#if (not defined(__FLOAT_FFTW) and (defined(__CUDA) || defined(__RCOM)))
    if (fft_device == "gpu")
    {
        fft_precision = "double";
    }
#endif

    // initialize pw_rho
    pw_rho = new ModulePW::PW_Basis_Big(fft_device, fft_precision);
    pw_rho_flag = true;

    // initialize pw_rhod
    if (double_grid)
    {
        pw_rhod = new ModulePW::PW_Basis_Big(fft_device, fft_precision);
    }
    else
    {
        pw_rhod = pw_rho;
    }

    // initialize pw_big
    pw_big = static_cast<ModulePW::PW_Basis_Big*>(pw_rhod);
    pw_big->setbxyz(inp.bx, inp.by, inp.bz);

    //! initialie the plane wave basis for rho
#ifdef __MPI
    pw_rho->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif

    //! for OFDFT calculations
    if (classname == "ESolver_OF" || inp.of_ml_gene_data == 1)
    {
        pw_rho->setfullpw(inp.of_full_pw, inp.of_full_pw_dim);
    }

    //! initialize the FFT grid
    // NOTE(liuyu): ref_cell_factor is currently forced to 1.0 in
    // read_input_item_md.cpp because the reference-cell mechanism is
    // disabled. When ref_cell_factor != 1, the PW_Basis lattice members
    // (lat0/tpiba/G/GGT/omega) become stale relative to ucell in NPT,
    // breaking sum_rho/get_local_pp_energy/cal_delta_escf/makov_payne.
    // If the reference-cell feature is re-enabled in the future, this
    // call site (and the equivalent in setup_pwwfc.cpp) MUST be updated
    // to use the proposed initgrids_ref/initgrids_actual split so that
    // only nx/ny/nz come from the reference cell while lat0/tpiba/G/GGT/omega
    // track the physical cell. Until then, ref_cell_factor * ucell.lat0
    // below is effectively just ucell.lat0.
    if (inp.nx * inp.ny * inp.nz == 0)
    {
        pw_rho->initgrids(inp.ref_cell_factor * ucell.lat0, ucell.latvec, 4.0 * inp.ecutwfc);
    }
    else
    {
        pw_rho->initgrids(inp.ref_cell_factor * ucell.lat0, ucell.latvec, inp.nx, inp.ny, inp.nz);
    }

    pw_rho->initparameters(false, 4.0 * inp.ecutwfc);
    pw_rho->fft_bundle.initfftmode(inp.fft_mode);
    pw_rho->setuptransform();
    pw_rho->collect_local_pw();
    pw_rho->collect_uniqgg();

    //! initialize the double grid (for uspp) if necessary
    if (double_grid)
    {
        ModulePW::PW_Basis_Sup* pw_rhod_sup = static_cast<ModulePW::PW_Basis_Sup*>(pw_rhod);
#ifdef __MPI
        pw_rhod->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
        if (classname == "ESolver_OF")
        {
            pw_rhod->setfullpw(inp.of_full_pw, inp.of_full_pw_dim);
        }
        // NOTE(liuyu): same ref_cell_factor warning applies to pw_rhod.
        if (inp.ndx * inp.ndy * inp.ndz == 0)
        {
            pw_rhod->initgrids(inp.ref_cell_factor * ucell.lat0, ucell.latvec, inp.ecutrho);
        }
        else
        {
            pw_rhod->initgrids(inp.ref_cell_factor * ucell.lat0, ucell.latvec, inp.ndx, inp.ndy, inp.ndz);
        }
        pw_rhod->initparameters(false, inp.ecutrho);
        pw_rhod->fft_bundle.initfftmode(inp.fft_mode);
        pw_rhod_sup->setuptransform(pw_rho);
        pw_rhod->collect_local_pw();
        pw_rhod->collect_uniqgg();
    }

    ModuleIO::print_rhofft(pw_rhod, pw_rho, pw_big, GlobalV::ofs_running);

    return;
}


void pw::teardown_pwrho(bool &pw_rho_flag,
        const bool double_grid,
        ModulePW::PW_Basis* &pw_rho, // pw for rhod
        ModulePW::PW_Basis* &pw_rhod) // pw for rhod
{
    if (pw_rho_flag == true)
    {
        delete pw_rho;
        pw_rho = nullptr;
        pw_rho_flag = false;
    }

    if (double_grid == true)
    {
        delete pw_rhod;
        pw_rhod = nullptr;
    }

   return;
}

