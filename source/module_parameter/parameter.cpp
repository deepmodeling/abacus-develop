#include "parameter.h"

Parameter PARAM;

Parameter::Parameter()
{
}
Parameter::~Parameter()
{
}
const Parameter& Parameter::get() const
{
    return *this;
}
const MD_para& Parameter::get_mdp() const
{
    return this->mdp;
}