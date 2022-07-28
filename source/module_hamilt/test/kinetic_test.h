#include "../../module_cell/unitcell_pseudo.h"

UnitCell::UnitCell(){}
UnitCell::~UnitCell(){}
UnitCell_pseudo::UnitCell_pseudo(){}
UnitCell_pseudo::~UnitCell_pseudo(){}
Magnetism::Magnetism(){}
Magnetism::~Magnetism(){}

void UnitCell_pseudo::setup_cell(const std::string &s_pseudopot_dir, 
		const std::string &fn, 
		std::ofstream &log)
{
    std::ifstream ifa(fn.c_str());
	if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_CONSTANT") )
	{
		ModuleBase::GlobalFunc::READ_VALUE(ifa, lat0);
		this->tpiba  = ModuleBase::TWO_PI / lat0;
		this->tpiba2 = tpiba * tpiba;
    }

    if( ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "LATTICE_VECTORS") )
    {
        // Reading lattice vectors. notice
        // here that only one cpu read these
        // parameters.
        ifa >> latvec.e11 >> latvec.e12;
        ModuleBase::GlobalFunc::READ_VALUE(ifa, latvec.e13);
        ifa >> latvec.e21 >> latvec.e22;
        ModuleBase::GlobalFunc::READ_VALUE(ifa, latvec.e23);
        ifa >> latvec.e31 >> latvec.e32;
        ModuleBase::GlobalFunc::READ_VALUE(ifa, latvec.e33);
    }

	//========================================================
	// Calculate unit cell volume
	// the reason to calculate volume here is 
	// Firstly, latvec must be read in.
	//========================================================
	assert(lat0 > 0.0);
	this->omega = abs( latvec.Det() ) * this->lat0 * lat0 * lat0 ;
		
	//==========================================================
	// Calculate recip. lattice vectors and dot products
	// latvec have the unit of lat0, but G has the unit 2Pi/lat0
	//==========================================================
	this->GT = latvec.Inverse();
	this->G  = GT.Transpose();
	this->GGT = G * GT;
	this->invGGT = GGT.Inverse();

	//this->cal_meshx();

	return;
}