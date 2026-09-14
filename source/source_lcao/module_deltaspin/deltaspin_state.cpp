#include "deltaspin_state.h"

#include "source_base/constants.h"
#include "source_base/tool_quit.h"

#include <cmath>

namespace spinconstrain
{

/**
 * @brief Calculate the spin constraint energy: E_scon = -sum_i (lambda_i . Mi_i).
 *
 * @details The constraint energy is the Lagrange multiplier term in the
 * constrained DFT functional:
 *   E'[rho] = E_DFT[rho] - sum_i lambda_i . (Mi_i - M_target_i)
 *
 * @return Constraint energy in Ry
 */
double ScState::cal_escon()
{
    this->escon_ = 0.0;
    if (this->lambda_.empty() || this->Mi_.empty())
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

double ScState::get_escon() const
{
    return this->escon_;
}

// set atomCounts
void ScState::set_atomCounts(const std::map<int, int>& atomCounts_in)
{
    this->atomCounts.clear();
    this->atomCounts = atomCounts_in;
}

// get atomCounts
const std::map<int, int>& ScState::get_atomCounts() const
{
    return this->atomCounts;
}

/// set nspin
void ScState::set_nspin(int nspin_in)
{
    if (nspin_in != 4 && nspin_in != 2)
    {
        ModuleBase::WARNING_QUIT("ScState::set_nspin", "nspin must be 2 or 4");
    }
    this->nspin_ = nspin_in;
}

/// get nspin
int ScState::get_nspin() const
{
    return this->nspin_;
}

void ScState::set_npol(int npol)
{
    this->npol_ = npol;
}

int ScState::get_npol() const
{
    return this->npol_;
}

int ScState::get_nw() const
{
    int nw = 0;
    for (const auto& pair : this->orbitalCounts)
    {
        nw += pair.second;
    }
    return nw;
}

/**
 * @brief Convert (itype, local_atom_index, orbital_index) to global orbital index.
 *
 * @details The global orbital index is used to access elements in distributed
 * matrices (ScaLAPACK format). The mapping is:
 *   iwt = sum_{t < itype} orbitalCounts[t]  +  iat * orbitalCounts[itype]  +  orbital_index
 * where iat = get_iat(itype, local_atom_index).
 *
 * @return Global orbital index, or 0 if itype not found
 */
int ScState::get_iwt(int itype, int iat, int orbital_index) const
{
    auto it1 = this->orbitalCounts.find(itype);
    if (it1 == this->orbitalCounts.end())
    {
        return 0;
    }
    int offset = 0;
    for (auto it = this->orbitalCounts.begin(); it != it1; ++it)
    {
        offset += it->second;
    }
    auto it2 = this->atomCounts.find(itype);
    if (it2 == this->atomCounts.end())
    {
        return offset;
    }
    return offset + iat * it1->second + orbital_index;
}

/// @brief Get total number of atoms across all element types
int ScState::get_nat() const
{
    int nat = 0;
    for (std::map<int, int>::const_iterator it = this->atomCounts.begin(); it != this->atomCounts.end(); ++it)
    {
        nat += it->second;
    }
    return nat;
}

/// @brief Get number of element types
int ScState::get_ntype() const
{
    return this->atomCounts.size();
}

/**
 * @brief Validate atom count data integrity.
 *
 * @details Checks that atomCounts has been properly initialized and contains
 * valid data. Called before any operation that depends on atom indexing.
 *
 * @par Error conditions
 * - "atomCounts is not set": init_sc() was not called
 * - "nat <= 0": no atoms in the system
 * - "itype out of range": element type index exceeds ntype
 * - "number of atoms <= 0": some element type has no atoms
 */
void ScState::check_atomCounts() const
{
    if (!this->atomCounts.size())
    {
        ModuleBase::WARNING_QUIT("ScState::check_atomCounts", "atomCounts is not set");
    }
    if (this->get_nat() <= 0)
    {
        ModuleBase::WARNING_QUIT("ScState::check_atomCounts", "nat <= 0");
    }
    for (std::map<int, int>::const_iterator it = this->atomCounts.begin(); it != this->atomCounts.end(); ++it)
    {
        int itype = it->first;
        if (itype < 0 || itype >= this->get_ntype())
        {
            ModuleBase::WARNING_QUIT("ScState::check_atomCounts", "itype out of range [0, ntype)");
        }
        int inat = it->second;
        if (inat <= 0)
        {
            ModuleBase::WARNING_QUIT("ScState::check_atomCounts", "number of atoms <= 0 for some element");
        }
    }
}

/**
 * @brief Convert (element_type, local_atom_index) to global atom index.
 *
 * @details Atoms in ABACUS are organized by element type. Within each type,
 * atoms are indexed locally (0, 1, ..., nat_itype-1). This function maps
 * to the global index that runs across all atoms (0, 1, ..., nat-1).
 *
 * Example: If type 0 has 2 Fe atoms and type 1 has 3 O atoms:
 *   get_iat(0, 0) -> 0 (Fe_0)
 *   get_iat(0, 1) -> 1 (Fe_1)
 *   get_iat(1, 0) -> 2 (O_0)
 *   get_iat(1, 1) -> 3 (O_1)
 *   get_iat(1, 2) -> 4 (O_2)
 *
 * @param itype Element type index (0 to ntype-1)
 * @param atom_index Local index within the element type
 * @return Global atom index
 */
int ScState::get_iat(int itype, int atom_index) const
{
    if (itype < 0 || itype >= this->get_ntype())
    {
        ModuleBase::WARNING_QUIT("ScState::get_iat", "itype out of range [0, ntype)");
    }
    if (atom_index < 0 || atom_index >= this->atomCounts.at(itype))
    {
        ModuleBase::WARNING_QUIT("ScState::get_iat", "atom index out of range [0, nat)");
    }
    int iat = 0;
    for (std::map<int, int>::const_iterator it = this->atomCounts.begin(); it != this->atomCounts.end(); ++it)
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
void ScState::set_orbitalCounts(const std::map<int, int>& orbitalCounts_in)
{
    this->orbitalCounts.clear();
    this->orbitalCounts = orbitalCounts_in;
}

// get orbitalCounts
const std::map<int, int>& ScState::get_orbitalCounts() const
{
    return this->orbitalCounts;
}

// set lnchiCounts
void ScState::set_lnchiCounts(const std::map<int, std::map<int, int>>& lnchiCounts_in)
{
    this->lnchiCounts.clear();
    this->lnchiCounts = lnchiCounts_in;
}

// get lnchiCounts
const std::map<int, std::map<int, int>>& ScState::get_lnchiCounts() const
{
    return this->lnchiCounts;
}

// set sc_lambda from ScData (parsed from STRU file)
// ScData is organized by element type; this function flattens it to per-atom arrays
void ScState::set_sc_lambda()
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

/**
 * @brief Set target magnetic moments from ScData (parsed from STRU file).
 *
 * @details Supports two specification modes:
 * - mag_type=0: Direct Cartesian (mx, my, mz) in uB
 * - mag_type=1: Spherical (|M|, theta, phi) converted to Cartesian:
 *   Mx = |M| * sin(theta) * cos(phi)
 *   My = |M| * sin(theta) * sin(phi)
 *   Mz = |M| * cos(theta)
 *   Angles are in degrees and converted to radians.
 *
 * Near-zero components (< 1e-14) are explicitly set to 0.0 to avoid
 * floating-point noise in constraint checks.
 */
void ScState::set_target_mag()
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

/**
 * @brief Set constraint flags from ScData.
 *
 * @details The constrain array determines which components of each atom's
 * magnetic moment are actively constrained:
 * - constrain[ia].x = 1: Mx is constrained to target_mag[ia].x
 * - constrain[ia].y = 1: My is constrained to target_mag[ia].y
 * - constrain[ia].z = 1: Mz is constrained to target_mag[ia].z
 * - constrain[ia].c = 0: component is free (determined by the system)
 *
 * Default is all zeros (no constraints). Components with constrain=0
 * are excluded from the lambda optimization and RMS error calculation.
 */
void ScState::set_constrain()
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
void ScState::set_sc_lambda(const ModuleBase::Vector3<double>* lambda_in, int nat_in)
{
    this->check_atomCounts();
    int nat = this->get_nat();
    if (nat_in != nat)
    {
        ModuleBase::WARNING_QUIT("ScState::set_sc_lambda", "lambda_in size mismatch with nat");
    }
    this->lambda_.resize(nat);
    for (int iat = 0; iat < nat; ++iat)
    {
        this->lambda_[iat] = lambda_in[iat];
    }
}

// set target_mag from variable
void ScState::set_target_mag(const ModuleBase::Vector3<double>* target_mag_in, int nat_in)
{
    this->check_atomCounts();
    int nat = this->get_nat();
    if (nat_in != nat)
    {
        ModuleBase::WARNING_QUIT("ScState::set_target_mag", "target_mag_in size mismatch with nat");
    }
    this->target_mag_.resize(nat);
    for (int iat = 0; iat < nat; ++iat)
    {
        this->target_mag_[iat] = target_mag_in[iat];
    }
}

void ScState::set_target_mag(const std::vector<ModuleBase::Vector3<double>>& target_mag_in)
{
    int nat = this->get_nat();
    assert(target_mag_in.size() == nat);
    if (this->nspin_ == 2)
    {
        this->target_mag_.resize(nat, 0.0);
        for (int iat = 0; iat < nat; iat++)
        {
            this->target_mag_[iat].z
                = target_mag_in[iat].z;
        }
    }
    else if (this->nspin_ == 4)
    {
        this->target_mag_ = target_mag_in;
    }
    else
    {
        ModuleBase::WARNING_QUIT("ScState::set_target_mag", "nspin must be 2 or 4");
    }
}

/// set constrain from variable
void ScState::set_constrain(const ModuleBase::Vector3<int>* constrain_in, int nat_in)
{
    this->check_atomCounts();
    int nat = this->get_nat();
    if (nat_in != nat)
    {
        ModuleBase::WARNING_QUIT("ScState::set_constrain", "constrain_in size mismatch with nat");
    }
    this->constrain_.resize(nat);
    for (int iat = 0; iat < nat; ++iat)
    {
        this->constrain_[iat] = constrain_in[iat];
    }
}

const std::vector<ModuleBase::Vector3<double>>& ScState::get_sc_lambda() const
{
    return this->lambda_;
}

const std::vector<ModuleBase::Vector3<double>>& ScState::get_target_mag() const
{
    return this->target_mag_;
}

/// get_constrain
const std::vector<ModuleBase::Vector3<int>>& ScState::get_constrain() const
{
    return this->constrain_;
}

/// @brief Reset all atomic magnetic moments to zero. Called before each Mi calculation.
void ScState::zero_Mi()
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
double ScState::get_decay_grad(int itype) const
{
    std::map<int, double>::const_iterator it = this->ScDecayGrad.find(itype);
    return it != this->ScDecayGrad.end() ? it->second : 0.0;
}

/// set grad_decy
void ScState::set_decay_grad()
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
const std::vector<double>& ScState::get_decay_grad() const
{
    return this->decay_grad_;
}

/// set grad_decy from variable
void ScState::set_decay_grad(const double* decay_grad_in, int ntype_in)
{
    this->check_atomCounts();
    int ntype = this->get_ntype();
    if (ntype_in != ntype)
    {
        ModuleBase::WARNING_QUIT("ScState::set_decay_grad", "decay_grad_in size mismatch with ntype");
    }
    this->decay_grad_.resize(ntype);
    for (int itype = 0; itype < ntype; ++itype)
    {
        this->decay_grad_[itype] = decay_grad_in[itype];
    }
}

/// @brief  set input parameters
void ScState::set_input_parameters(double sc_thr_in,
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
double ScState::get_sc_thr() const
{
    return this->sc_thr_;
}

/// get current adaptive sc threshold
double ScState::get_current_sc_thr() const
{
    return this->current_sc_thr_;
}

/// get computed magnetic moments Mi per atom
const std::vector<ModuleBase::Vector3<double>>& ScState::get_Mi() const
{
    return this->Mi_;
}

/// get human-readable atom labels for table printing
const std::vector<std::string>& ScState::get_atomLabels() const
{
    return this->atomLabels_;
}

/// get nsc
int ScState::get_nsc() const
{
    return this->nsc_;
}

/// get nsc_min
int ScState::get_nsc_min() const
{
    return this->nsc_min_;
}

/// get alpha_trial
double ScState::get_alpha_trial() const
{
    return this->alpha_trial_;
}

/// get sccut
double ScState::get_sccut() const
{
    return this->restrict_current_;
}

/// set sc_drop_thr
void ScState::set_sc_drop_thr(double sc_drop_thr_in)
{
    this->sc_drop_thr_ = sc_drop_thr_in;
}

/// get sc_drop_thr
double ScState::get_sc_drop_thr() const
{
    return this->sc_drop_thr_;
}

} // namespace spinconstrain
