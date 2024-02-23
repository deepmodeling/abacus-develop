#include "module_cell/module_symmetry/symmetry.h"
#include "module_cell/atom_spec.h"
#include "module_cell/unitcell.h"


//Add json objects to init
namespace Json
{
#ifdef __RAPIDJSON
// void gen_init(ModuleSymmetry::Symmetry *symm,Atom *atoms);
void gen_init(UnitCell *ucell);

void add_nkstot(int nkstot,int nkstot_ibz);


#endif
}