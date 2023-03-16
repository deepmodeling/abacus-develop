#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "../../variable_cell.h"

namespace GlobalC
{
	UnitCell ucell;
    Structure_Factor sf;
    ModulePW::PW_Basis* rhopw;
}

UnitCell::UnitCell(){};
UnitCell::~UnitCell(){};

void UnitCell::remake_cell(){};

void UnitCell::update_pos_taud(double* posd_in)
{
    int iat = 0;
    for (int it = 0; it < this->ntype; it++)
    {
        Atom* atom = &this->atoms[it];
        for (int ia = 0; ia < atom->na; ia++)
        {
            this->atoms[it].taud[ia].x += posd_in[iat*3];
            this->atoms[it].taud[ia].y += posd_in[iat*3 + 1];
            this->atoms[it].taud[ia].z += posd_in[iat*3 + 2];
            iat++;
        }
    }
    assert(iat == this->nat);
    this->periodic_boundary_adjustment();
}

void UnitCell::periodic_boundary_adjustment()
{
    //----------------------------------------------
    // because of the periodic boundary condition
    // we need to adjust the atom positions,
    // first adjust direct coordinates,
    // then update them into cartesian coordinates,
    //----------------------------------------------
    for (int it = 0; it < this->ntype; it++)
    {
        Atom* atom = &this->atoms[it];
        for (int ia = 0; ia < atom->na; ia++)
        {
            // mohan update 2011-03-21
            if (atom->taud[ia].x < 0)
                atom->taud[ia].x += 1.0;
            if (atom->taud[ia].y < 0)
                atom->taud[ia].y += 1.0;
            if (atom->taud[ia].z < 0)
                atom->taud[ia].z += 1.0;
            if (atom->taud[ia].x >= 1.0)
                atom->taud[ia].x -= 1.0;
            if (atom->taud[ia].y >= 1.0)
                atom->taud[ia].y -= 1.0;
            if (atom->taud[ia].z >= 1.0)
                atom->taud[ia].z -= 1.0;

            atom->tau[ia] = atom->taud[ia] * this->latvec;
        }
    }
    return;
}

void UnitCell::print_stru_file(const std::string &fn, const int &type, const int &level)const {};
void UnitCell::print_tau(void)const{};
void UnitCell::setup_cell_after_vc(std::ofstream &log){};
void UnitCell::print_cell_cif(const std::string& fn) const{};

Magnetism::Magnetism(){};
Magnetism::~Magnetism(){};

Structure_Factor::Structure_Factor(){};
Structure_Factor::~Structure_Factor(){};
void Structure_Factor::setup_structure_factor(UnitCell* Ucell, ModulePW::PW_Basis* rho_basis){};

void Variable_Cell::init_after_vc(){};

extern Input INPUT;