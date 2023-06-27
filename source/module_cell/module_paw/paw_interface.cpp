#include "paw_element.h"
#include "module_base/tool_title.h"
#include "module_base/tool_quit.h"

void Paw_Element::init_paw(const double ecutwfc_in, const double cell_factor_in)
{
    this -> ecutwfc = ecutwfc_in;
    this -> cell_factor = cell_factor_in;
}