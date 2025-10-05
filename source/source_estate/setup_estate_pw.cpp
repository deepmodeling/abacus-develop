#include "source_estate/setup_estate_pw.h"

template <typename T, typename Device>
void elecstate::setup_estate_pw(UnitCell& ucell, // unitcell
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rho, // pw for rho
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
        const ModulePW::PW_Basis_Big* pw_big, // pw for big grid
        surchem &solvent, //  solvent
		const Input_para& inp) // input parameters
{
    ModuleBase::TITLE("elecstate", "setup_estate_pw");

    //! Initialize ElecState, set pelec pointer
    if (pelec == nullptr)
    {
        if (inp.esolver_type == "sdft")
        {
            //! SDFT only supports double precision currently
            pelec = new elecstate::ElecStatePW_SDFT<std::complex<double>, Device>(pw_wfc,
                &chr, &kv, &ucell, &ppcell,
                pw_rhod, pw_rho, pw_big);
        }
        else
        {
            pelec = new elecstate::ElecStatePW<T, Device>(pw_wfc,
                &chr, &kv, &ucell, &ppcell,
                pw_rhod, pw_rho, pw_big);
        }
    }

    //! Set the cell volume variable in pelec
    pelec->omega = ucell.omega;

    //! Inititlize the charge density.
    chr.allocate(inp.nspin);

    //! Initialize DFT-1/2
    if (PARAM.inp.dfthalf_type > 0)
    {
        vsep_cell = new VSep;
        vsep_cell->init_vsep(*pw_rhod, ucell.sep_cell);
    }

    //! Initialize the potential.
    if (pelec->pot == nullptr)
    {
        pelec->pot = new elecstate::Potential(pw_rhod,
              pw_rho, &ucell, &locpp.vloc, &sf,
              &solvent, &(pelec->f_en.etxc), &(pelec->f_en.vtxc), vsep_cell);
    }

    //! Initalize local pseudopotential
    locpp.init_vloc(ucell, pw_rhod);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    //! Initalize non-local pseudopotential
    ppcell.init(ucell, &sf, pw_wfc);
    ppcell.init_vnl(ucell, pw_rhod);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

    //! Setup occupations
    if (inp.ocp)
    {
        elecstate::fixed_weights(inp.ocp_kb,
                                 inp.nbands,
                                 inp.nelec,
                                 pelec->klist,
                                 pelec->wg,
                                 pelec->skip_weights);
    }

    return;
}


void elecstate::teardown_estate_pw(elecstate::ElecState *pelec, VSep* vsep_cell) 
{
    ModuleBase::TITLE("elecstate", "teardown_estate_pw");

    if (vsep_cell != nullptr)
    {
        delete vsep_cell;
    }

    if (pelec != nullptr)
    {
        delete reinterpret_cast<elecstate::ElecStatePW<T, Device>*>(pelec);
        pelec = nullptr;
    }
}


template void elecstate::setup_estate_pw<std::complex<float>, base_device::DEVICE_CPU>(
        UnitCell& ucell, // unitcell
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rho, // pw for rho
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
        const ModulePW::PW_Basis_Big* pw_big, // pw for big grid
        surchem &solvent, //  solvent
		const Input_para& inp); // input parameters

template void elecstate::setup_estate_pw<std::complex<double>, base_device::DEVICE_CPU>(
        UnitCell& ucell, // unitcell
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rho, // pw for rho
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
        const ModulePW::PW_Basis_Big* pw_big, // pw for big grid
        surchem &solvent, //  solvent
		const Input_para& inp); // input parameters



#if ((defined __CUDA) || (defined __ROCM))

template void elecstate::setup_estate_pw<std::complex<float>, base_device::DEVICE_GPU>(
        UnitCell& ucell, // unitcell
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rho, // pw for rho
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
        const ModulePW::PW_Basis_Big* pw_big, // pw for big grid
        surchem &solvent, //  solvent
		const Input_para& inp); // input parameters

template void elecstate::setup_estate_pw<std::complex<double>, base_device::DEVICE_GPU>(
        UnitCell& ucell, // unitcell
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rho, // pw for rho
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
        const ModulePW::PW_Basis_Big* pw_big, // pw for big grid
        surchem &solvent, //  solvent
		const Input_para& inp); // input parameters

#endif
