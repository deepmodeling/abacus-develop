#include "source_pw/module_pwdft/setup_pwrho.h"

#include "source_estate/module_charge/symmetry_rho.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_pw/module_pwdft/onsite_projector.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_pw/module_pwdft/VSep_in_pw.h"

template <typename T, typename Device>
void pw::setup_pot(const int istep, 
		UnitCell& ucell, // unitcell 
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Parallel_Grid &para_grid, // parallel of FFT grids
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		psi::Psi<T, Device>* kspw_psi, // electronic wave functions
        hamilt::Hamilt<T, Device>* p_hamilt, // hamiltonian
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
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
#if (not defined(__ENABLE_FLOAT_FFTW) and (defined(__CUDA) || defined(__RCOM)))
    if (fft_device == "gpu")
    {
        fft_precision = "double";
    }
#endif

    // initialize pw_rho
    pw_rho = new ModulePW::PW_Basis_Big(fft_device, fft_precision);
    pw_rho_flag = true;

    // initialize pw_rhod
    if (PARAM.globalv.double_grid)
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

    // setup the structure factors
    sf.set(pw_rhod, inp.nbspline);

    //! read pseudopotentials, move somewhere else??? mohan note 20251005
    elecstate::read_pseudo(GlobalV::ofs_running, ucell);

    //! initialie the plane wave basis for rho
#ifdef __MPI
    this->pw_rho->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif

    //! for OFDFT calculations
    if (this->classname == "ESolver_OF" || inp.of_ml_gene_data == 1)
    {
        this->pw_rho->setfullpw(inp.of_full_pw, inp.of_full_pw_dim);
    }

    //! initialize the FFT grid
    if (inp.nx * inp.ny * inp.nz == 0)
    {
        this->pw_rho->initgrids(inp.ref_cell_factor * ucell.lat0, ucell.latvec, 4.0 * inp.ecutwfc);
    }
    else
    {
        this->pw_rho->initgrids(inp.ref_cell_factor * ucell.lat0, ucell.latvec, inp.nx, inp.ny, inp.nz);
    }

    this->pw_rho->initparameters(false, 4.0 * inp.ecutwfc);
    this->pw_rho->fft_bundle.initfftmode(inp.fft_mode);
    this->pw_rho->setuptransform();
    this->pw_rho->collect_local_pw();
    this->pw_rho->collect_uniqgg();

    //! initialize the double grid (for uspp) if necessary
    if ( PARAM.globalv.double_grid)
    {
        ModulePW::PW_Basis_Sup* pw_rhod_sup = static_cast<ModulePW::PW_Basis_Sup*>(pw_rhod);
#ifdef __MPI
        this->pw_rhod->initmpi(GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL, POOL_WORLD);
#endif
        if (this->classname == "ESolver_OF")
        {
            this->pw_rhod->setfullpw(inp.of_full_pw, inp.of_full_pw_dim);
        }
        if (inp.ndx * inp.ndy * inp.ndz == 0)
        {
            this->pw_rhod->initgrids(inp.ref_cell_factor * ucell.lat0, ucell.latvec, inp.ecutrho);
        }
        else
        {
            this->pw_rhod->initgrids(inp.ref_cell_factor * ucell.lat0, ucell.latvec, inp.ndx, inp.ndy, inp.ndz);
        }
        this->pw_rhod->initparameters(false, inp.ecutrho);
        this->pw_rhod->fft_bundle.initfftmode(inp.fft_mode);
        pw_rhod_sup->setuptransform(this->pw_rho);
        this->pw_rhod->collect_local_pw();
        this->pw_rhod->collect_uniqgg();
    }

    return;
}

template void pw::setup_pot<std::complex<float>, base_device::DEVICE_CPU>(
        const int istep,  // ionic step
		UnitCell& ucell, // unitcell 
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Parallel_Grid &para_grid, // parallel of FFT grids
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		psi::Psi<std::complex<float>, base_device::DEVICE_CPU>* kspw_psi, // electronic wave functions
        hamilt::Hamilt<std::complex<float>, base_device::DEVICE_CPU>* p_hamilt, // hamiltonian
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
		const Input_para& inp); // input parameters


template void pw::setup_pot<std::complex<double>, base_device::DEVICE_CPU>(
        const int istep,  // ionic step
		UnitCell& ucell, // unitcell 
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Parallel_Grid &para_grid, // parallel of FFT grids
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		psi::Psi<std::complex<double>, base_device::DEVICE_CPU>* kspw_psi, // electronic wave functions
        hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>* p_hamilt, // hamiltonian
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
		const Input_para& inp); // input parameters

#if ((defined __CUDA) || (defined __ROCM))

template void pw::setup_pot<std::complex<float>, base_device::DEVICE_GPU>(
        const int istep,  // ionic step
		UnitCell& ucell, // unitcell 
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Parallel_Grid &para_grid, // parallel of FFT grids
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		psi::Psi<std::complex<float>, base_device::DEVICE_GPU>* kspw_psi, // electronic wave functions
        hamilt::Hamilt<std::complex<float>, base_device::DEVICE_GPU>* p_hamilt, // hamiltonian
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
		const Input_para& inp); // input parameters

template void pw::setup_pot<std::complex<double>, base_device::DEVICE_GPU>(
        const int istep,  // ionic step
		UnitCell& ucell, // unitcell 
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Parallel_Grid &para_grid, // parallel of FFT grids
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* kspw_psi, // electronic wave functions
        hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>* p_hamilt, // hamiltonian
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
		const Input_para& inp); // input parameters

#endif
