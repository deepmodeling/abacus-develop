#include "module_parameter/parameter.h"
#include "unitcell.h"

void UnitCell::print_unitcell_pseudo(const std::string& fn)
{
    if (PARAM.inp.test_pseudo_cell) {
        ModuleBase::TITLE("UnitCell", "print_unitcell_pseudo");
}
    std::ofstream ofs(fn.c_str());

    this->print_cell(ofs);
    for (int i = 0; i < ntype; i++)
    {
        atoms[i].print_Atom(ofs);
    }

    ofs.close();
    return;
}
