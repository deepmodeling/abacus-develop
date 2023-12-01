#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"

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
}

void UnitCell::print_stru_file(const std::string &fn, const int &type, const int &level)const {};
void UnitCell::print_tau(void)const{};
void UnitCell::setup_cell_after_vc(std::ofstream &log){};
void UnitCell::print_cell_cif(const std::string& fn) const{};

Magnetism::Magnetism(){};
Magnetism::~Magnetism(){};

Atom::Atom(){};
Atom::~Atom(){};
Atom_pseudo::Atom_pseudo(){};
Atom_pseudo::~Atom_pseudo(){};
pseudo::pseudo(){};
pseudo::~pseudo(){};
int ModuleSymmetry::Symmetry::symm_flag = 0;
void ModuleSymmetry::Symmetry::symmetrize_mat3(ModuleBase::matrix& sigma, const Lattice& lat)const {};
void ModuleSymmetry::Symmetry::symmetrize_vec3_nat(double* v)const {};
Structure_Factor::Structure_Factor() {};
Structure_Factor::~Structure_Factor(){};
void Structure_Factor::setup_structure_factor(UnitCell* Ucell, const ModulePW::PW_Basis* rho_basis){};

extern Input INPUT;