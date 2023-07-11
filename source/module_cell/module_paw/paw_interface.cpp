#include "paw_element.h"
#include "module_base/tool_title.h"
#include "module_base/tool_quit.h"

void Paw_Element::init_paw_element(const double ecutwfc_in, const double cell_factor_in)
{
    this -> ecutwfc = ecutwfc_in;
    this -> cell_factor = cell_factor_in;
}

double Paw_Element::get_ptilde(const int istate_in, const double q_in)
{
    return this->splint(qgrid, ptilde_q[istate_in], d2ptilde_q[istate_in], q_in);
}