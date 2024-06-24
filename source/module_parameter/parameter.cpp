#include "parameter.h"

Parameter PARAM;

const Input_para& Parameter::get() const
{
    return this->input;
}
const MD_para& Parameter::get_mdp() const
{
    return this->input.mdp;
}
const Other_para& Parameter::globalV() const
{
    return this->gv;
}