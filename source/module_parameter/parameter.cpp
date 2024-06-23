#include "parameter.h"

Parameter PARAM;

Parameter::Parameter()
{
}
Parameter::~Parameter()
{
}

std::string Parameter::get_suffix() const
{
    return suffix;
}
std::string Parameter::get_latname() const
{
    return latname;
}
std::string Parameter::get_stru_file() const
{
    return stru_file;
}
std::string Parameter::get_kpoint_file() const
{
    return kpoint_file;
}
std::string Parameter::get_pseudo_dir() const
{
    return pseudo_dir;
}
std::string Parameter::get_orbital_dir() const
{
    return orbital_dir;
}
double Parameter::get_pseudo_rcut() const
{
    return pseudo_rcut;
}
bool Parameter::get_pseudo_mesh() const
{
    return pseudo_mesh;
}
int Parameter::get_lmaxmax() const
{
    return lmaxmax;
}
std::string Parameter::get_dft_functional() const
{
    return dft_functional;
}
double Parameter::get_xc_temperature() const
{
    return xc_temperature;
}
std::string Parameter::get_calculation() const
{
    return calculation;
}
std::string Parameter::get_esolver_type() const
{
    return esolver_type;
}

