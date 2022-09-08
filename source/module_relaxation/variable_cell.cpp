#include "variable_cell.h"
#include "../src_pw/global.h"
#include "../input.h"

Variable_Cell::Variable_Cell(){}
Variable_Cell::~Variable_Cell(){}

void Variable_Cell::init_after_vc()
{
	ModuleBase::TITLE("Variable_Cell","init_after_vc");

    GlobalC::ucell.setup_cell_after_vc(GlobalV::ofs_running);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    if(ModuleSymmetry::Symmetry::symm_flag)
    {
        GlobalC::symm.analy_sys(GlobalC::ucell, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    GlobalC::kv.set_after_vc(GlobalC::symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, GlobalC::ucell.G, GlobalC::ucell.latvec);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

    
    //only G-vector and K-vector are changed due to the change of lattice vector
    //FFT grids do not change!!
    GlobalC::rhopw->initgrids(GlobalC::ucell.lat0, GlobalC::ucell.latvec, GlobalC::rhopw->nx, GlobalC::rhopw->ny, GlobalC::rhopw->nz);
    GlobalC::rhopw->collect_local_pw(); 
    GlobalC::rhopw->collect_uniqgg();
    GlobalC::sf.setup_structure_factor(&GlobalC::ucell,GlobalC::rhopw);

    GlobalV::ofs_running << " Setup the Vl+Vh+Vxc according to new structure factor and new charge." << std::endl;
    //=================================
    // initalize local pseudopotential
    //=================================
    GlobalC::ppcell.init_vloc(GlobalC::ppcell.vloc, GlobalC::rhopw);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"LOCAL POTENTIAL");

    //======================================
    // Initalize non local pseudopotential
    //======================================
    if(GlobalV::BASIS_TYPE=="pw")
    {
        GlobalC::ppcell.init_vnl(GlobalC::ucell);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"NON-LOCAL POTENTIAL");
    }

    return;
}
