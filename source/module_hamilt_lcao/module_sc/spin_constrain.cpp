#include "spin_constrain.h"

SpinConstrain& SpinConstrain::getInstance() {
    static SpinConstrain instance; // Guaranteed to be created and destroyed only once
    return instance;
}

// set atomCounts
void SpinConstrain::set_atomCounts(const std::map<int, int>& atomCounts_in) {
    this->atomCounts = atomCounts_in;
}

// get atomCounts
const std::map<int, int>& SpinConstrain::get_atomCounts() const
{
    return this->atomCounts;
}

/// set npol
void SpinConstrain::set_npol(int npol)
{
    this->npol_ = npol;
}

/// get npol
int SpinConstrain::get_npol()
{
    return this->npol_;
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
int SpinConstrain::get_iat(int itype, int atom_index)
{
    if (itype < 0 || itype >= this->get_ntype())
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::get_iat","itype out of range [0, ntype)");
    }
    if (atom_index < 0 || atom_index >= this->atomCounts[itype])
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::get_iat","atom index out of range [0, nat)");
    }
    int iat = 0;
    for (std::map<int, int>::iterator it = this->atomCounts.begin(); it != this->atomCounts.end(); ++it) {
        if (it->first == itype)
        {
            break;
        }
        iat += it->second;
    }
    iat += atom_index;
    return iat;
}

// clear atomCounts
void SpinConstrain::clear_atomCounts()
{
    this->atomCounts.clear();
}

// clear orbitalCounts
void SpinConstrain::clear_orbitalCounts()
{
    this->orbitalCounts.clear();
}

// set orbitalCounts
void SpinConstrain::set_orbitalCounts(const std::map<int, int>& orbitalCounts_in) {
    this->orbitalCounts = orbitalCounts_in;
}

// get orbitalCounts
const std::map<int, int>& SpinConstrain::get_orbitalCounts() const
{
    return this->orbitalCounts;
}

int SpinConstrain::get_nw()
{
    this->check_atomCounts();
    int nw = 0;
    for (std::map<int, int>::iterator it = this->orbitalCounts.begin(); it != this->orbitalCounts.end(); ++it) {
        nw += (it->second)*this->atomCounts[it->first]*this->npol_;
    }
    return nw;
}

int SpinConstrain::get_iwt(int itype, int iat, int orbital_index)
{
    this->check_atomCounts();
    if (itype < 0 || itype >= this->get_ntype())
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::get_iwt","itype out of range [0, ntype)");
    }
    if (iat < 0 || iat >= this->get_nat())
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::get_iwt","iat out of range [0, nat)");
    }
    if (orbital_index < 0 || orbital_index >= this->orbitalCounts[itype])
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::get_iwt","orbital index out of range [0, atom_nw)");
    }
    int iwt = 0;
    for (std::map<int, int>::iterator it = this->orbitalCounts.begin(); it != this->orbitalCounts.end(); ++it) {
        if (it->first == itype)
        {
            break;
        }
        iwt += (it->second)*this->atomCounts[it->first]*this->npol_;
    }
    for (int i = 0; i < iat; ++i) {
        iwt += this->orbitalCounts[itype]*this->npol_;
    }
    iwt += orbital_index*this->npol_;
    return iwt;
}

// set sc_lambda from ScData
void SpinConstrain::set_sc_lambda()
{
    this->check_atomCounts();
    int nat = this->get_nat();
    this->lambda_.resize(nat);
    for (auto& itype_data : this->ScData) {
        int itype = itype_data.first;
        for (auto& element_data : itype_data.second) {
            int index = element_data.index;
            int iat = this->get_iat(itype, index);
            ModuleBase::Vector3<double> lambda;
            lambda.x = element_data.lambda[0];
            lambda.y = element_data.lambda[1];
            lambda.z = element_data.lambda[2];
            this->lambda_[iat] = lambda;
        }
    }
}

// set sc_mag from ScData
void SpinConstrain::set_sc_mag()
{
    this->check_atomCounts();
    int nat = this->get_nat();
    this->sc_mag_.resize(nat);
    for (auto& itype_data : this->ScData) {
        int itype = itype_data.first;
        for (auto& element_data : itype_data.second) {
            int index = element_data.index;
            int iat = this->get_iat(itype, index);
            ModuleBase::Vector3<double> mag;
            mag.x = element_data.sc_mag[0];
            mag.y = element_data.sc_mag[1];
            mag.z = element_data.sc_mag[2];
            this->sc_mag_[iat] = mag;
        }
    }
}

const std::vector<ModuleBase::Vector3<double>>& SpinConstrain::get_sc_lambda() const
{
    return this->lambda_;
}

const std::vector<ModuleBase::Vector3<double>>& SpinConstrain::get_sc_mag() const
{
    return this->sc_mag_;
}