#include "source_pw/module_pwdft/structure_factor.h"

namespace GlobalC
{
    ModulePW::PW_Basis* rhopw;
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
SepPot::SepPot(){}
SepPot::~SepPot(){}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}
int ModuleSymmetry::Symmetry::symm_flag = 0;
void ModuleSymmetry::Symmetry::symmetrize_mat3(ModuleBase::matrix& sigma, const Lattice& lat)const {};
void ModuleSymmetry::Symmetry::symmetrize_vec3_nat(double* v)const {};
Structure_Factor::Structure_Factor() {};
Structure_Factor::~Structure_Factor(){};
void Structure_Factor::setup_structure_factor(const UnitCell* Ucell, const Parallel_Grid&, const ModulePW::PW_Basis* rho_basis){};
