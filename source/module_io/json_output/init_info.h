#include "module_cell/module_symmetry/symmetry.h"
#include "module_cell/atom_spec.h"
#include "module_cell/unitcell.h"


/**
* @brief In this part of the code to complete the init part of the json tree.
*/
namespace Json
{
#ifdef __RAPIDJSON
// void gen_init(ModuleSymmetry::Symmetry *symm,Atom *atoms);

/**
* @param ucell: ucell for reading json parameters.
*/
void gen_init(UnitCell *ucell);

/**
* @param nkstot,nkstot_ibz: two param in json tree
*/
void add_nkstot(int nkstot,int nkstot_ibz);

/**
* @param ucell: ucell for reading structure init in abacus.
*/
void gen_stru(UnitCell *ucell);
#endif
}