#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"

namespace GlobalC
{
    Structure_Factor sf;
    ModulePW::PW_Basis* rhopw;
}

UnitCell::UnitCell(){};
UnitCell::~UnitCell(){};


void UnitCell::print_stru_file(const std::string& fn, 
                               const int& nspin,
                               const bool& direct,
                               const bool& vel,
                               const bool& magmom,
                               const bool& orb,
                               const bool& dpks_desc,
                               const int& iproc) const {};
void UnitCell::print_tau()const{};

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
void Structure_Factor::setup_structure_factor(const UnitCell* Ucell, const Parallel_Grid&, const ModulePW::PW_Basis* rho_basis){};