//==========================================================
// Used to genarate ucell
//==========================================================
#ifndef SETCELL_H
#define SETCELL_H

#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"
#include "module_io/input.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"

namespace GlobalC
{
	  UnitCell ucell;
    ModulePW::PW_Basis* rhopw;
    Parallel_Grid Pgrid;
}

UnitCell::UnitCell(){};
UnitCell::~UnitCell(){};

Magnetism::Magnetism(){};
Magnetism::~Magnetism(){};

Atom::Atom(){};
Atom::~Atom(){};
Atom_pseudo::Atom_pseudo(){};
Atom_pseudo::~Atom_pseudo(){};
pseudo::pseudo(){};
pseudo::~pseudo(){};
/*
Structure_Factor::Structure_Factor(){};
Structure_Factor::~Structure_Factor(){};
void Structure_Factor::setup_structure_factor(UnitCell* Ucell, const ModulePW::PW_Basis* rho_basis){};
*/
class Setcell
{
  public:
    static void setupcell(UnitCell &ucell)
    {
        ucell.ntype = 2;

        ucell.atoms = new Atom[ucell.ntype];
        ucell.set_atom_flag = true;

        ucell.atoms[0].ncpp.psd = "H";
        ucell.atoms[1].ncpp.psd = "O";

        ucell.lat0 = 1;
        ucell.lat0_angstrom = ucell.lat0 * 0.529177;
        ucell.tpiba = ModuleBase::TWO_PI / ucell.lat0;
        ucell.tpiba2 = ucell.tpiba * ucell.tpiba;

        ucell.latvec.e11 = ucell.latvec.e22 = ucell.latvec.e33 = 10;
        ucell.latvec.e12 = ucell.latvec.e13 = ucell.latvec.e23 = 0;
        ucell.latvec.e21 = ucell.latvec.e31 = ucell.latvec.e32 = 0;

        ucell.a1.x = ucell.latvec.e11;
        ucell.a1.y = ucell.latvec.e12;
        ucell.a1.z = ucell.latvec.e13;

        ucell.a2.x = ucell.latvec.e21;
        ucell.a2.y = ucell.latvec.e22;
        ucell.a2.z = ucell.latvec.e23;

        ucell.a3.x = ucell.latvec.e31;
        ucell.a3.y = ucell.latvec.e32;
        ucell.a3.z = ucell.latvec.e33;

        ucell.nat = 3;
        ucell.atoms[0].na = 2;
        ucell.atoms[1].na = 1;

        ucell.atoms[0].tau = new ModuleBase::Vector3<double>[2];
        ucell.atoms[1].tau = new ModuleBase::Vector3<double>[1];

        ucell.atoms[0].tau[0].set(7.5456, 0, 9.54275);
        ucell.atoms[0].tau[1].set(7.542, 1.8495, 7.34175);
        ucell.atoms[1].tau[0].set(7.54965, 0, 7.48585);

        ucell.atoms[0].ncpp.zv = 1;
        ucell.atoms[1].ncpp.zv = 6;

        ucell.omega = std::abs(ucell.latvec.Det()) * ucell.lat0 * ucell.lat0 * ucell.lat0;

    };
};

#endif