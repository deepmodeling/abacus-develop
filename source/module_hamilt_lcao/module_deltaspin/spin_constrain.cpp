#include "spin_constrain.h"

#include "module_base/formatter.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/dspin_lcao.h"

#include <cmath>

namespace spinconstrain
{

template <typename FPTYPE>
SpinConstrain<FPTYPE>& SpinConstrain<FPTYPE>::getScInstance()
{
    static SpinConstrain<FPTYPE> instance; // Guaranteed to be created and destroyed only once
    return instance;
}

template <typename FPTYPE>
double SpinConstrain<FPTYPE>::cal_escon()
{
    this->escon_ = 0.0;
    if (!this->is_Mi_converged)
    {
        return this->escon_;
    }
    int nat = this->get_nat();
    for (int iat = 0; iat < nat; iat++)
    {
        this->escon_ -= this->lambda_[iat].x * this->Mi_[iat].x;
        this->escon_ -= this->lambda_[iat].y * this->Mi_[iat].y;
        this->escon_ -= this->lambda_[iat].z * this->Mi_[iat].z;
    }
    return this->escon_;
}

template <typename FPTYPE>
double SpinConstrain<FPTYPE>::get_escon()
{
    return this->escon_;
}

// set atomCounts
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_atomCounts(const std::map<int, int>& atomCounts_in)
{
    this->atomCounts.clear();
    this->atomCounts = atomCounts_in;
}

// get atomCounts
template <typename FPTYPE>
const std::map<int, int>& SpinConstrain<FPTYPE>::get_atomCounts() const
{
    return this->atomCounts;
}

/// set nspin
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_nspin(int nspin_in)
{
    if (nspin_in != 4 && nspin_in != 2)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_nspin", "nspin must be 2 or 4");
    }
    this->nspin_ = nspin_in;
}

/// get nspin
template <typename FPTYPE>
int SpinConstrain<FPTYPE>::get_nspin()
{
    return this->nspin_;
}

template <typename FPTYPE>
int SpinConstrain<FPTYPE>::get_nat()
{
    int nat = 0;
    for (std::map<int, int>::iterator it = this->atomCounts.begin(); it != this->atomCounts.end(); ++it)
    {
        nat += it->second;
    }
    return nat;
}

template <typename FPTYPE>
int SpinConstrain<FPTYPE>::get_ntype()
{
    return this->atomCounts.size();
}

template <typename FPTYPE>
void SpinConstrain<FPTYPE>::check_atomCounts()
{
    if (!this->atomCounts.size())
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::check_atomCounts", "atomCounts is not set");
    }
    if (this->get_nat() <= 0)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::check_atomCounts", "nat <= 0");
    }
    for (std::map<int, int>::iterator it = this->atomCounts.begin(); it != this->atomCounts.end(); ++it)
    {
        int itype = it->first;
        if (itype < 0 || itype >= this->get_ntype())
        {
            ModuleBase::WARNING_QUIT("SpinConstrain::check_atomCounts", "itype out of range [0, ntype)");
        }
        int inat = it->second;
        if (inat <= 0)
        {
            ModuleBase::WARNING_QUIT("SpinConstrain::check_atomCounts", "number of atoms <= 0 for some element");
        }
    }
}

// get iat
template <typename FPTYPE>
int SpinConstrain<FPTYPE>::get_iat(int itype, int atom_index)
{
    if (itype < 0 || itype >= this->get_ntype())
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::get_iat", "itype out of range [0, ntype)");
    }
    if (atom_index < 0 || atom_index >= this->atomCounts[itype])
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::get_iat", "atom index out of range [0, nat)");
    }
    int iat = 0;
    for (std::map<int, int>::iterator it = this->atomCounts.begin(); it != this->atomCounts.end(); ++it)
    {
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
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_orbitalCounts(const std::map<int, int>& orbitalCounts_in)
{
    this->orbitalCounts.clear();
    this->orbitalCounts = orbitalCounts_in;
}

// get orbitalCounts
template <typename FPTYPE>
const std::map<int, int>& SpinConstrain<FPTYPE>::get_orbitalCounts() const
{
    return this->orbitalCounts;
}

// set lnchiCounts
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_lnchiCounts(const std::map<int, std::map<int, int>>& lnchiCounts_in)
{
    this->lnchiCounts.clear();
    this->lnchiCounts = lnchiCounts_in;
}

// get lnchiCounts
template <typename FPTYPE>
const std::map<int, std::map<int, int>>& SpinConstrain<FPTYPE>::get_lnchiCounts() const
{
    return this->lnchiCounts;
}

// set sc_lambda from ScData
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_sc_lambda()
{
    this->check_atomCounts();
    int nat = this->get_nat();
    this->lambda_.resize(nat);
    for (auto& itype_data: this->ScData)
    {
        int itype = itype_data.first;
        for (auto& element_data: itype_data.second)
        {
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

// set target_mag from ScData
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_target_mag()
{
    this->check_atomCounts();
    int nat = this->get_nat();
    this->target_mag_.resize(nat, 0.0);
    for (auto& itype_data: this->ScData)
    {
        int itype = itype_data.first;
        for (auto& element_data: itype_data.second)
        {
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
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_constrain()
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
    for (auto& itype_data: this->ScData)
    {
        int itype = itype_data.first;
        for (auto& element_data: itype_data.second)
        {
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
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_sc_lambda(const ModuleBase::Vector3<double>* lambda_in, int nat_in)
{
    this->check_atomCounts();
    int nat = this->get_nat();
    if (nat_in != nat)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_sc_lambda", "lambda_in size mismatch with nat");
    }
    this->lambda_.resize(nat);
    for (int iat = 0; iat < nat; ++iat)
    {
        this->lambda_[iat] = lambda_in[iat];
    }
}

// set target_mag from variable
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_target_mag(const ModuleBase::Vector3<double>* target_mag_in, int nat_in)
{
    this->check_atomCounts();
    int nat = this->get_nat();
    if (nat_in != nat)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_target_mag", "target_mag_in size mismatch with nat");
    }
    this->target_mag_.resize(nat);
    for (int iat = 0; iat < nat; ++iat)
    {
        this->target_mag_[iat] = target_mag_in[iat];
    }
}

template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_target_mag(const std::vector<ModuleBase::Vector3<double>>& target_mag_in)
{
    int nat = this->get_nat();
    assert(target_mag_in.size() == nat);
    if (this->nspin_ == 2)
    {
        this->target_mag_.resize(nat, 0.0);
        for (int iat = 0; iat < nat; iat++)
        {
            this->target_mag_[iat].z
                = target_mag_in[iat].x; /// this is wired because the UnitCell class set in x direction
        }
    }
    else if (this->nspin_ == 4)
    {
        this->target_mag_ = target_mag_in;
    }
    else
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_target_mag", "nspin must be 2 or 4");
    }
}

/// set constrain from variable
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_constrain(const ModuleBase::Vector3<int>* constrain_in, int nat_in)
{
    this->check_atomCounts();
    int nat = this->get_nat();
    if (nat_in != nat)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_constrain", "constrain_in size mismatch with nat");
    }
    this->constrain_.resize(nat);
    for (int iat = 0; iat < nat; ++iat)
    {
        this->constrain_[iat] = constrain_in[iat];
    }
}

template <typename FPTYPE>
const std::vector<ModuleBase::Vector3<double>>& SpinConstrain<FPTYPE>::get_sc_lambda() const
{
    return this->lambda_;
}

template <typename FPTYPE>
const std::vector<ModuleBase::Vector3<double>>& SpinConstrain<FPTYPE>::get_target_mag() const
{
    return this->target_mag_;
}

/// get_constrain
template <typename FPTYPE>
const std::vector<ModuleBase::Vector3<int>>& SpinConstrain<FPTYPE>::get_constrain() const
{
    return this->constrain_;
}

/// zero atomic magnetic moment
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::zero_Mi()
{
    this->check_atomCounts();
    int nat = this->get_nat();
    this->Mi_.resize(nat);
    for (int iat = 0; iat < nat; ++iat)
    {
        this->Mi_[iat].x = 0.0;
        this->Mi_[iat].y = 0.0;
        this->Mi_[iat].z = 0.0;
    }
}

/// get grad_decay
/// this function can only be called by the root process because only
/// root process reads the ScDecayGrad from json file
template <typename FPTYPE>
double SpinConstrain<FPTYPE>::get_decay_grad(int itype)
{
    return this->ScDecayGrad[itype];
}

/// set grad_decy
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_decay_grad()
{
    this->check_atomCounts();
    int ntype = this->get_ntype();
    this->decay_grad_.resize(ntype);
    for (int itype = 0; itype < ntype; ++itype)
    {
        this->decay_grad_[itype] = 0.0;
    }
}

/// get decay_grad
template <typename FPTYPE>
const std::vector<double>& SpinConstrain<FPTYPE>::get_decay_grad()
{
    return this->decay_grad_;
}

/// set grad_decy from variable
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_decay_grad(const double* decay_grad_in, int ntype_in)
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
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_input_parameters(double sc_thr_in,
                                                 int nsc_in,
                                                 int nsc_min_in,
                                                 double alpha_trial_in,
                                                 double sccut_in,
                                                 double sc_drop_thr_in)
{
    this->sc_thr_ = sc_thr_in;
    this->nsc_ = nsc_in;
    this->nsc_min_ = nsc_min_in;
    this->alpha_trial_ = alpha_trial_in / ModuleBase::Ry_to_eV;
    this->restrict_current_ = sccut_in / ModuleBase::Ry_to_eV;
    this->sc_drop_thr_ = sc_drop_thr_in;
}

/// get sc_thr
template <typename FPTYPE>
double SpinConstrain<FPTYPE>::get_sc_thr()
{
    return this->sc_thr_;
}

/// get nsc
template <typename FPTYPE>
int SpinConstrain<FPTYPE>::get_nsc()
{
    return this->nsc_;
}

/// get nsc_min
template <typename FPTYPE>
int SpinConstrain<FPTYPE>::get_nsc_min()
{
    return this->nsc_min_;
}

/// get alpha_trial
template <typename FPTYPE>
double SpinConstrain<FPTYPE>::get_alpha_trial()
{
    return this->alpha_trial_;
}

/// get sccut
template <typename FPTYPE>
double SpinConstrain<FPTYPE>::get_sccut()
{
    return this->restrict_current_;
}

/// set sc_drop_thr
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_sc_drop_thr(double sc_drop_thr_in)
{
    this->sc_drop_thr_ = sc_drop_thr_in;
}

/// get sc_drop_thr
template <typename FPTYPE>
double SpinConstrain<FPTYPE>::get_sc_drop_thr()
{
    return this->sc_drop_thr_;
}

template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_solver_parameters(K_Vectors& kv_in,
                                                  void* p_hamilt_in,
                                                  void* psi_in,
                                                  elecstate::ElecState* pelec_in,
                                                  std::string KS_SOLVER_in)
{
    this->kv_ = kv_in;
    this->p_hamilt = p_hamilt_in;
    this->psi = psi_in;
    this->pelec = pelec_in;
    this->KS_SOLVER = KS_SOLVER_in;
}

/// @brief  set ParaV
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::set_ParaV(Parallel_Orbitals* ParaV_in)
{
    this->ParaV = ParaV_in;
    int nloc = this->ParaV->nloc;
    if (nloc <= 0)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_ParaV", "nloc <= 0");
    }
}

/// print Mi
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::print_Mi(std::ofstream& ofs_running)
{
    this->check_atomCounts();
    int nat = this->get_nat();
    std::vector<double> mag_x(nat, 0.0);
    std::vector<double> mag_y(nat, 0.0);
    std::vector<double> mag_z(nat, 0.0);
    if (this->nspin_ == 2)
    {
        const std::vector<std::string> title = {"Total Magnetism (uB)", ""};
        const std::vector<std::string> fmts = {"%-26s", "%20.10f"};
        FmtTable table(title, nat, fmts, {FmtTable::Align::RIGHT, FmtTable::Align::LEFT});
        for (int iat = 0; iat < nat; ++iat)
        {
            mag_z[iat] = Mi_[iat].z;
        }
        table << this->atomLabels_ << mag_z;
        ofs_running << table.str() << std::endl;
    }
    else if (this->nspin_ == 4)
    {
        const std::vector<std::string> title = {"Total Magnetism (uB)", "", "", ""};
        const std::vector<std::string> fmts = {"%-26s", "%20.10f", "%20.10f", "%20.10f"};
        FmtTable table(title, nat, fmts, {FmtTable::Align::RIGHT, FmtTable::Align::LEFT});
        for (int iat = 0; iat < nat; ++iat)
        {
            mag_x[iat] = Mi_[iat].x;
            mag_y[iat] = Mi_[iat].y;
            mag_z[iat] = Mi_[iat].z;
        }
        table << this->atomLabels_ << mag_x << mag_y << mag_z;
        ofs_running << table.str() << std::endl;
    }
}

/// print magnetic force (defined as \frac{\delta{L}}/{\delta{Mi}} = -lambda[iat])
template <typename FPTYPE>
void SpinConstrain<FPTYPE>::print_Mag_Force(std::ofstream& ofs_running)
{
    this->check_atomCounts();
    int nat = this->get_nat();
    std::vector<double> mag_force_x(nat, 0.0);
    std::vector<double> mag_force_y(nat, 0.0);
    std::vector<double> mag_force_z(nat, 0.0);
    if (this->nspin_ == 2)
    {
        const std::vector<std::string> title = {"Magnetic force (eV/uB)", ""};
        const std::vector<std::string> fmts = {"%-26s", "%20.10f"};
        FmtTable table(title, nat, fmts, {FmtTable::Align::RIGHT, FmtTable::Align::LEFT});
        for (int iat = 0; iat < nat; ++iat)
        {
            mag_force_z[iat] = lambda_[iat].z * ModuleBase::Ry_to_eV;
        }
        table << this->atomLabels_ << mag_force_z;
        ofs_running << table.str() << std::endl;
    }
    else if (this->nspin_ == 4)
    {
        const std::vector<std::string> title = {"Magnetic force (eV/uB)", "", "", ""};
        const std::vector<std::string> fmts = {"%-26s", "%20.10f", "%20.10f", "%20.10f"};
        FmtTable table(title, nat, fmts, {FmtTable::Align::RIGHT, FmtTable::Align::LEFT});
        for (int iat = 0; iat < nat; ++iat)
        {
            mag_force_x[iat] = lambda_[iat].x * ModuleBase::Ry_to_eV;
            mag_force_y[iat] = lambda_[iat].y * ModuleBase::Ry_to_eV;
            mag_force_z[iat] = lambda_[iat].z * ModuleBase::Ry_to_eV;
        }
        table << this->atomLabels_ << mag_force_x << mag_force_y << mag_force_z;
        ofs_running << table.str() << std::endl;
    }
}

template class SpinConstrain<std::complex<double>>;
template class SpinConstrain<double>;

} // namespace spinconstrain