#include "spin_constrain.h"

#include <cmath>

template<typename FPTYPE, typename Device>
SpinConstrain<FPTYPE, Device>& SpinConstrain<FPTYPE, Device>::getScInstance() {
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
    return this->escon_;
}

template<typename FPTYPE, typename Device>
double SpinConstrain<FPTYPE, Device>::get_escon()
{
    return this->escon_;
}

// set atomCounts
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_atomCounts(const std::map<int, int>& atomCounts_in) {
    this->atomCounts.clear();
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
    if (nspin_in != 4 && nspin_in != 2)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_nspin","nspin must be 2 or 4");
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

// set orbitalCounts
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_orbitalCounts(const std::map<int, int>& orbitalCounts_in) {
    this->orbitalCounts.clear();
    this->orbitalCounts = orbitalCounts_in;
}

// get orbitalCounts
template<typename FPTYPE, typename Device>
const std::map<int, int>& SpinConstrain<FPTYPE, Device>::get_orbitalCounts() const
{
    return this->orbitalCounts;
}

// set lnchiCounts
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_lnchiCounts(const std::map<int, std::map<int, int>>& lnchiCounts_in) {
    this->lnchiCounts.clear();
    this->lnchiCounts = lnchiCounts_in;
}

// get lnchiCounts
template<typename FPTYPE, typename Device>
const std::map<int, std::map<int, int>>& SpinConstrain<FPTYPE, Device>::get_lnchiCounts() const
{
    return this->lnchiCounts;
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

// set target_mag from ScData
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_target_mag()
{
    this->check_atomCounts();
    int nat = this->get_nat();
    this->target_mag_.resize(nat, 0.0);
    for (auto& itype_data : this->ScData) {
        int itype = itype_data.first;
        for (auto& element_data : itype_data.second) {
            int index = element_data.index;
            int iat = this->get_iat(itype, index);
            ModuleBase::Vector3<double> mag(0.0, 0.0, 0.0);
            if (element_data.mag_type == 0)
            {
                mag.x = element_data.target_mag[0];
                mag.y = element_data.target_mag[1];
                mag.z = element_data.target_mag[2];
            }
            else if (element_data.mag_type == 1)
            {
                double radian_angle1 = element_data.target_mag_angle1 * M_PI / 180.0;
                double radian_angle2 = element_data.target_mag_angle2 * M_PI / 180.0;
                mag.x = element_data.target_mag_val * std::sin(radian_angle1) * std::cos(radian_angle2);
                mag.y = element_data.target_mag_val * std::sin(radian_angle1) * std::sin(radian_angle2);
                mag.z = element_data.target_mag_val * std::cos(radian_angle1);
                if (std::abs(mag.x) < 1e-14)
                    mag.x = 0.0;
                if (std::abs(mag.y) < 1e-14)
                    mag.y = 0.0;
                if (std::abs(mag.z) < 1e-14)
                    mag.z = 0.0;
            }
            this->target_mag_[iat] = mag;
        }
    }
}

// set constrain from ScData
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_constrain()
{
    this->check_atomCounts();
    int nat = this->get_nat();
    this->constrain_.resize(nat);
    // constrain is 0 by default, which means no constrain
    // and the corresponding mag moments should be determined
    // by the physical nature of the system
    for (int iat = 0; iat < nat; iat++)
    {
            this->constrain_[iat].x = 0;
            this->constrain_[iat].y = 0;
            this->constrain_[iat].z = 0;
    }
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

// set target_mag from variable
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_target_mag(const ModuleBase::Vector3<double>* target_mag_in, int nat_in)
{
    this->check_atomCounts();
    int nat = this->get_nat();
    if (nat_in != nat)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_target_mag","target_mag_in size mismatch with nat");
    }
    this->target_mag_.resize(nat);
    for (int iat=0; iat < nat; ++iat)
    {
        this->target_mag_[iat] = target_mag_in[iat];
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

template<typename FPTYPE, typename Device>
const std::vector<ModuleBase::Vector3<double>>& SpinConstrain<FPTYPE, Device>::get_target_mag() const
{
    return this->target_mag_;
}

/// get_constrain
template<typename FPTYPE, typename Device>
const std::vector<ModuleBase::Vector3<int>>& SpinConstrain<FPTYPE, Device>::get_constrain() const
{
    return this->constrain_;
}

/// zero atomic magnetic moment
template<typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::zero_Mi()
{
    this->check_atomCounts();
    int nat = this->get_nat();
    this->Mi_.resize(nat);
    for (int iat=0; iat < nat; ++iat)
    {
        this->Mi_[iat].x = 0.0;
        this->Mi_[iat].y = 0.0;
        this->Mi_[iat].z = 0.0;
    }
}

/// get grad_decay
/// this function can only be called by the root process because only
/// root process reads the ScDecayGrad from json file
template <typename FPTYPE, typename Device>
double SpinConstrain<FPTYPE, Device>::get_decay_grad(int itype)
{
    return this->ScDecayGrad[itype];
}

/// set grad_decy
template <typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_decay_grad()
{
    this->check_atomCounts();
    int ntype = this->get_ntype();
    this->decay_grad_.resize(ntype);
    for (int itype = 0; itype < ntype; ++itype)
    {
        this->decay_grad_[itype] = 0.0;
    }
    if (this->decay_grad_switch_)
    {
        for (auto& itype_data: this->ScDecayGrad)
        {
            int itype = itype_data.first;
            this->decay_grad_[itype] = itype_data.second * ModuleBase::Ry_to_eV;
        }
    }
}

/// get decay_grad
template <typename FPTYPE, typename Device>
const std::vector<double>& SpinConstrain<FPTYPE, Device>::get_decay_grad()
{
    return this->decay_grad_;
}

/// set grad_decy from variable
template <typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_decay_grad(const double* decay_grad_in, int ntype_in)
{
    this->check_atomCounts();
    int ntype = this->get_ntype();
    if (ntype_in != ntype)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_decay_grad", "decay_grad_in size mismatch with ntype");
    }
    this->decay_grad_.resize(ntype);
    for (int itype = 0; itype < ntype; ++itype)
    {
        this->decay_grad_[itype] = decay_grad_in[itype];
    }
}

/// @brief  set input parameters
template <typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_input_parameters(double sc_thr_in,
                                                         int nsc_in,
                                                         int nsc_min_in,
                                                         double alpha_trial_in,
                                                         double sccut_in,
                                                         bool decay_grad_switch_in)
{
    this->sc_thr_ = sc_thr_in;
    this->nsc_ = nsc_in;
    this->nsc_min_ = nsc_min_in;
    this->alpha_trial_ = alpha_trial_in / ModuleBase::Ry_to_eV;
    this->restrict_current_ = sccut_in / ModuleBase::Ry_to_eV;
    this->decay_grad_switch_ = decay_grad_switch_in;
}

/// get sc_thr
template <typename FPTYPE, typename Device>
double SpinConstrain<FPTYPE, Device>::get_sc_thr()
{
    return this->sc_thr_;
}

/// get nsc
template <typename FPTYPE, typename Device>
int SpinConstrain<FPTYPE, Device>::get_nsc()
{
    return this->nsc_;
}

/// get nsc_min
template <typename FPTYPE, typename Device>
int SpinConstrain<FPTYPE, Device>::get_nsc_min()
{
    return this->nsc_min_;
}

/// get alpha_trial
template <typename FPTYPE, typename Device>
double SpinConstrain<FPTYPE, Device>::get_alpha_trial()
{
    return this->alpha_trial_;
}

/// get sccut
template <typename FPTYPE, typename Device>
double SpinConstrain<FPTYPE, Device>::get_sccut()
{
    return this->restrict_current_;
}

/// set decay_grad_switch
template <typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_decay_grad_switch(bool decay_grad_switch_in)
{
    this->decay_grad_switch_ = decay_grad_switch_in;
}

/// get decay_grad_switch
template <typename FPTYPE, typename Device>
bool SpinConstrain<FPTYPE, Device>::get_decay_grad_switch()
{
    return this->decay_grad_switch_;
}

template <typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_solver_parameters(K_Vectors kv_in,
                                                          hsolver::HSolver<FPTYPE, Device>* phsol_in,
                                                          hamilt::Hamilt<FPTYPE, Device>* p_hamilt_in,
                                                          psi::Psi<FPTYPE>* psi_in,
                                                          elecstate::ElecState* pelec_in,
                                                          std::string KS_SOLVER_in,
                                                          LCAO_Matrix* LM_in)
{
    this->kv_ = kv_in;
    this->phsol = phsol_in;
    this->p_hamilt = p_hamilt_in;
    this->psi = psi_in;
    this->pelec = pelec_in;
    this->KS_SOLVER = KS_SOLVER_in;
    this->LM = LM_in;
}

/// @brief  set ParaV
template <typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::set_ParaV(Parallel_Orbitals* ParaV_in)
{
    this->ParaV = ParaV_in;
    int nloc = this->ParaV->nloc;
    if (nloc <= 0)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_ParaV", "nloc <= 0");
    }
}

/// print Mi
template <typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::print_Mi(bool print)
{
    this->check_atomCounts();
    int nat = this->get_nat();
    if (print)
    {
        for (int iat = 0; iat < nat; ++iat)
        {
            if (this->nspin_ == 2)
            {
                std::cout << "Total Magnetism on atom: " << iat << " " << std::setprecision(10) << " (" << Mi_[iat].z << ")" << std::endl;
            }
            else if (this->nspin_ ==4)
            {
                std::cout << "Total Magnetism on atom: " << iat << " " << std::setprecision(10) << " (" << Mi_[iat].x
                        << ", " << Mi_[iat].y << ", " << Mi_[iat].z << ")" << std::endl;
            }
        }
    }
}

/// save mi in log file
template <typename FPTYPE, typename Device>
void SpinConstrain<FPTYPE, Device>::log_Mi_info(std::ofstream& ofs)
{
    this->check_atomCounts();
    int nat = this->get_nat();
    ofs << std::endl;
    ofs << "-----------------------------------------------------------------------" << std::endl;
    for (int iat = 0; iat < nat; ++iat)
    {
        if (this->nspin_ == 2)
        {
            ofs << "Total Magnetism on atom (with weight function): " << iat << " " << std::setprecision(10) << " ("
                << Mi_[iat].z << ")" << std::endl;
        }
        else if (this->nspin_ == 4)
        {
            ofs << "Total Magnetism on atom (with weight function): " << iat << " " << std::setprecision(10) << " ("
                << Mi_[iat].x << ", " << Mi_[iat].y << ", " << Mi_[iat].z << ")" << std::endl;
        }
    }
    ofs << "-----------------------------------------------------------------------" << std::endl;
    ofs << std::endl;
}

template class SpinConstrain<std::complex<double>, psi::DEVICE_CPU>;
template class SpinConstrain<double, psi::DEVICE_CPU>;