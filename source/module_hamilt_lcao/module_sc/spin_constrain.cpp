#include "spin_constrain.h"

template<typename FPTYPE, typename Device>
SpinConstrain<FPTYPE, Device>& SpinConstrain<FPTYPE, Device>::getInstance() {
    static SpinConstrain<FPTYPE, Device> instance; // Guaranteed to be created and destroyed only once
    return instance;
}

template<typename FPTYPE, typename Device>
double SpinConstrain<FPTYPE, Device>::cal_escon()
{
    this->escon_ = 0.0;
    int nat = this->get_nat();
    for (int iat = 0; iat < nat; iat++)
    {
        this->escon_ += this->lambda_[iat].x * this->Mi_[iat].x;
        this->escon_ += this->lambda_[iat].y * this->Mi_[iat].y;
        this->escon_ += this->lambda_[iat].z * this->Mi_[iat].z;
    }
    //std::cout << "this->escon_ " << this->escon_ << std::endl;
    return this->escon_;
}

// set atomCounts
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_atomCounts(const std::map<int, int>& atomCounts_in) {
    this->atomCounts = atomCounts_in;
}

// get atomCounts
template<typename FPTYPE, typename Device>
const std::map<int, int>& SpinConstrain<FPTYPE, Device>::get_atomCounts() const
{
    return this->atomCounts;
}

/// set npol
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_npol(int npol)
{
    this->npol_ = npol;
}

/// get npol
template<typename FPTYPE, typename Device>
int SpinConstrain<FPTYPE, Device>::get_npol()
{
    return this->npol_;
}

/// set nspin
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_nspin(int nspin_in)
{
    if (nspin_in != 4)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_nspin","nspin must be 4 now");
    }
    this->nspin_ = nspin_in;
}

/// get nspin
template<typename FPTYPE, typename Device>
int SpinConstrain<FPTYPE, Device>::get_nspin()
{
    return this->nspin_;
}

template<typename FPTYPE, typename Device>
int SpinConstrain<FPTYPE, Device>::get_nat()
{
    int nat = 0;
    for (std::map<int, int>::iterator it = this->atomCounts.begin(); it != this->atomCounts.end(); ++it) {
        nat += it->second;
    }
    return nat;
}

template<typename FPTYPE, typename Device>
int SpinConstrain<FPTYPE, Device>::get_ntype()
{
    return this->atomCounts.size();
}

template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::check_atomCounts()
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
template<typename FPTYPE, typename Device>
int SpinConstrain<FPTYPE, Device>::get_iat(int itype, int atom_index)
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
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::clear_atomCounts()
{
    this->atomCounts.clear();
}

// clear orbitalCounts
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::clear_orbitalCounts()
{
    this->orbitalCounts.clear();
}

// set orbitalCounts
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_orbitalCounts(const std::map<int, int>& orbitalCounts_in) {
    this->orbitalCounts = orbitalCounts_in;
}

// get orbitalCounts
template<typename FPTYPE, typename Device>
const std::map<int, int>& SpinConstrain<FPTYPE, Device>::get_orbitalCounts() const
{
    return this->orbitalCounts;
}

template<typename FPTYPE, typename Device>
int SpinConstrain<FPTYPE, Device>::get_nw()
{
    this->check_atomCounts();
    int nw = 0;
    for (std::map<int, int>::iterator it = this->orbitalCounts.begin(); it != this->orbitalCounts.end(); ++it) {
        nw += (it->second)*this->atomCounts[it->first]*this->npol_;
    }
    return nw;
}

template<typename FPTYPE, typename Device>
int SpinConstrain<FPTYPE, Device>::get_iwt(int itype, int iat, int orbital_index)
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
    if (orbital_index < 0 || orbital_index >= this->orbitalCounts[itype]*this->npol_)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::get_iwt","orbital index out of range [0, atom_nw*npol)");
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
    iwt += orbital_index;
    return iwt;
}

// set sc_lambda from ScData
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_sc_lambda()
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
            this->lambda_[iat] = lambda*this->meV_to_Ry;
        }
    }
}

// set init_mag from ScData
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_init_mag()
{
    this->check_atomCounts();
    int nat = this->get_nat();
    this->init_mag_.resize(nat);
    for (auto& itype_data : this->ScData) {
        int itype = itype_data.first;
        for (auto& element_data : itype_data.second) {
            int index = element_data.index;
            int iat = this->get_iat(itype, index);
            ModuleBase::Vector3<double> init_mag;
            init_mag.x = element_data.init_mag[0];
            init_mag.y = element_data.init_mag[1];
            init_mag.z = element_data.init_mag[2];
            this->init_mag_[iat] = init_mag;
        }
    }
}

// set sc_mag from ScData
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_sc_mag()
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

// set sc_mag from ScData
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_constrain()
{
    this->check_atomCounts();
    int nat = this->get_nat();
    this->constrain_.resize(nat);
    for (auto& itype_data : this->ScData) {
        int itype = itype_data.first;
        for (auto& element_data : itype_data.second) {
            int index = element_data.index;
            int iat = this->get_iat(itype, index);
            ModuleBase::Vector3<int> constr;
            constr.x = element_data.constrain[0];
            constr.y = element_data.constrain[1];
            constr.z = element_data.constrain[2];
            this->constrain_[iat] = constr;
        }
    }
}

// set sc_lambda from variable
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_sc_lambda(const ModuleBase::Vector3<double>* lambda_in, int nat_in)
{
    this->check_atomCounts();
    int nat = this->get_nat();
    if (nat_in != nat)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_sc_lambda","lambda_in size mismatch with nat");
    }
    this->lambda_.resize(nat);
    for (int iat=0; iat < nat; ++iat)
    {
        this->lambda_[iat] = lambda_in[iat];
    }
}

// set init_mag from variable
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_init_mag(const ModuleBase::Vector3<double>* init_mag_in, int nat_in)
{
    this->check_atomCounts();
    int nat = this->get_nat();
    if (nat_in != nat)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_init_mag","init_mag_in size mismatch with nat");
    }
    this->init_mag_.resize(nat);
    for (int iat=0; iat < nat; ++iat)
    {
        this->init_mag_[iat] = init_mag_in[iat];
    }
}

// set sc_mag from variable
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_sc_mag(const ModuleBase::Vector3<double>* sc_mag_in, int nat_in)
{
    this->check_atomCounts();
    int nat = this->get_nat();
    if (nat_in != nat)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_sc_mag","sc_mag_in size mismatch with nat");
    }
    this->sc_mag_.resize(nat);
    for (int iat=0; iat < nat; ++iat)
    {
        this->sc_mag_[iat] = sc_mag_in[iat];
    }
}

/// set constrain from variable
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_constrain(const ModuleBase::Vector3<int>* constrain_in, int nat_in)
{
    this->check_atomCounts();
    int nat = this->get_nat();
    if (nat_in != nat)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_constrain","constrain_in size mismatch with nat");
    }
    this->constrain_.resize(nat);
    for (int iat=0; iat < nat; ++iat)
    {
        this->constrain_[iat] = constrain_in[iat];
    }
}

template<typename FPTYPE, typename Device>
const std::vector<ModuleBase::Vector3<double>>& SpinConstrain<FPTYPE, Device>::get_sc_lambda() const
{
    return this->lambda_;
}

/// get init_mag
template<typename FPTYPE, typename Device>
const std::vector<ModuleBase::Vector3<double>>& SpinConstrain<FPTYPE, Device>::get_init_mag() const
{
    return this->init_mag_;
}

template<typename FPTYPE, typename Device>
const std::vector<ModuleBase::Vector3<double>>& SpinConstrain<FPTYPE, Device>::get_sc_mag() const
{
    return this->sc_mag_;
}

/// get_constrain
template<typename FPTYPE, typename Device>
const std::vector<ModuleBase::Vector3<int>>& SpinConstrain<FPTYPE, Device>::get_constrain() const
{
    return this->constrain_;
}

template class SpinConstrain<double, psi::DEVICE_CPU>;