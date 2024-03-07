#include <ctime>
#include <string>

#include "module_cell/module_symmetry/symmetry.h"
#include "module_cell/atom_spec.h"
#include "module_io/input.h"
#include "module_cell/unitcell.h"
namespace Json
{

// void create_Json(ModuleSymmetry::Symmetry *symm,Atom *atoms,Input *input);

void create_Json(UnitCell *ucell,Input *input);

// Output the json to abacus.json file
void json_output();

// Convert time_t to string
void convert_time(std::time_t time_now, std::string& time_str);

// generate struture wrapper function
void gen_stru_wrapper(UnitCell *ucell);
} // namespace Json
