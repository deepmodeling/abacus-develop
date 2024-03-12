#ifndef FOR_TESTING_KLIST_H
#define FOR_TESTING_KLIST_H

#include "module_cell/parallel_kpoints.h"
#include "module_base/parallel_global.h"
#include "module_cell/klist.h"
#include "module_elecstate/magnetism.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "module_cell/atom_pseudo.h"
#include "module_cell/atom_spec.h"
#include "module_cell/unitcell.h"
#include "module_cell/pseudo.h"
#include "module_cell/setup_nonlocal.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#include "module_cell/parallel_kpoints.h"
#include "module_io/berryphase.h"
#include "module_basis/module_ao/ORB_gaunt_table.h"

bool berryphase::berry_phase_flag=0;

pseudo::pseudo(){}
pseudo::~pseudo(){}
Atom::Atom(){}
Atom::~Atom(){}
Atom_pseudo::Atom_pseudo(){}
Atom_pseudo::~Atom_pseudo(){}
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
UnitCell::UnitCell(){}
UnitCell::~UnitCell(){}
Magnetism::Magnetism(){}
Magnetism::~Magnetism(){}
ORB_gaunt_table::ORB_gaunt_table(){}
ORB_gaunt_table::~ORB_gaunt_table(){}
pseudopot_cell_vl::pseudopot_cell_vl(){}
pseudopot_cell_vl::~pseudopot_cell_vl(){}
pseudopot_cell_vnl::pseudopot_cell_vnl(){}
pseudopot_cell_vnl::~pseudopot_cell_vnl(){}
Soc::~Soc()
{
}
Fcoef::~Fcoef()
{
}

namespace GlobalC
{
	Parallel_Kpoints Pkpoints;
	UnitCell ucell;
}

#endif
