#ifndef FOR_TESTING_KLIST_H
#define FOR_TESTING_KLIST_H

#include "source_base/parallel_global.h"
#include "source_basis/module_ao/ORB_gaunt_table.h"
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/klist.h"
#include "source_cell/parallel_kpoints.h"
#include "source_cell/pseudo.h"
#include "source_cell/setup_nonlocal.h"
#include "source_cell/unitcell.h"
#include "source_estate/magnetism.h"
#include "source_pw/hamilt_pwdft/VL_in_pw.h"
#include "source_pw/hamilt_pwdft/VNL_in_pw.h"
#include "source_pw/hamilt_pwdft/parallel_grid.h"
#include "module_io/berryphase.h"

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



#endif
