#include "spin_constrain.h"

SpinConstrain& SpinConstrain::getInstance() {
    static SpinConstrain instance; // Guaranteed to be created and destroyed only once
    return instance;
}

// set atomCounts
void SpinConstrain::set_atomCounts(const std::map<int, int>& atomCounts) {
    this->atomCounts = atomCounts;
}

// get atomCounts
const std::map<int, int>& SpinConstrain::get_atomCounts() const
{
    return this->atomCounts;
}

int SpinConstrain::get_nat()
{
    int nat = 0;
    for (std::map<int, int>::iterator it = this->atomCounts.begin(); it != this->atomCounts.end(); ++it) {
        nat += it->second;
    }
    return nat;
}

int SpinConstrain::get_ntype()
{
    return this->atomCounts.size();
}

void SpinConstrain::check_atomCounts()
{
    if (!this->atomCounts.size())
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::check_atomCounts","atomCounts is not set");
    }
    if (this->get_nat() <= 0)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::check_atomCounts","nat <= 0");
    }
    for (std::map<int, int>::iterator it = this->atomCounts.begin(); it != this->atomCounts.end(); ++it) {
        int itype = it->first;
        if (itype < 0 || itype >= this->get_ntype())
        {
            ModuleBase::WARNING_QUIT("SpinConstrain::check_atomCounts","itype out of range [0, ntype)");
        }
        int inat = it->second;
        if (inat <= 0)
        {
            ModuleBase::WARNING_QUIT("SpinConstrain::check_atomCounts","number of atoms <= 0 for some element");
        }
    }
}

// get iat
int SpinConstrain::get_iat(int itype, int index)
{
    if (itype < 0 || itype >= this->get_ntype())
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::get_iat","itype out of range [0, ntype)");
    }
    if (index < 0 || index >= this->atomCounts[itype])
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::get_iat","index out of range [0, nat)");
    }
    int iat = 0;
    for (std::map<int, int>::iterator it = this->atomCounts.begin(); it != this->atomCounts.end(); ++it) {
        if (it->first == itype)
        {
            break;
        }
        iat += it->second;
    }
    iat += index;
    return iat;
}

// clear atomCounts
void SpinConstrain::clear_atomCounts()
{
    this->atomCounts.clear();
}

// get sc_lambda from ScData
std::vector<ModuleBase::Vector3<double>> SpinConstrain::get_sc_lambda()
{
    this->check_atomCounts();
    int nat = this->get_nat();
    std::vector<ModuleBase::Vector3<double>> sc_lambda(nat);
    for (auto& itype_data : this->ScData) {
        int itype = itype_data.first;
        for (auto& element_data : itype_data.second) {
            int index = element_data.index;
            int iat = this->get_iat(itype, index);
            ModuleBase::Vector3<double> lambda;
            lambda.x = element_data.lambda[0];
            lambda.y = element_data.lambda[1];
            lambda.z = element_data.lambda[2];
            sc_lambda[iat] = lambda;
        }
    }
    return sc_lambda;
}

// get sc_mag from ScData
std::vector<ModuleBase::Vector3<double>> SpinConstrain::get_sc_mag()
{
    this->check_atomCounts();
    int nat = this->get_nat();
    std::vector<ModuleBase::Vector3<double>> sc_mag(nat);
    for (auto& itype_data : this->ScData) {
        int itype = itype_data.first;
        for (auto& element_data : itype_data.second) {
            int index = element_data.index;
            int iat = this->get_iat(itype, index);
            ModuleBase::Vector3<double> mag;
            mag.x = element_data.sc_mag[0];
            mag.y = element_data.sc_mag[1];
            mag.z = element_data.sc_mag[2];
            sc_mag[iat] = mag;
        }
    }
    return sc_mag;
}