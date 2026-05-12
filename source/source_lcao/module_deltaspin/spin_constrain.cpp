#include "spin_constrain.h"

#include "source_base/formatter.h"
#include "source_lcao/module_operator_lcao/dspin_lcao.h"

#include <cmath>

namespace spinconstrain
{

/**
 * @brief Singleton instance accessor.
 *
 * @details Uses Meyers' Singleton pattern (local static variable).
 * Guaranteed thread-safe initialization in C++11 and later.
 * Each template instantiation (complex<double>, double) gets its own instance.
 */
template <typename TK>
SpinConstrain<TK>& SpinConstrain<TK>::getScInstance()
{
    static SpinConstrain<TK> instance; // Guaranteed to be created and destroyed only once
    return instance;
}

/**
 * @brief Calculate the spin constraint energy: E_scon = -sum_i (lambda_i . Mi_i).
 *
 * @details The constraint energy is the Lagrange multiplier term in the
 * constrained DFT functional:
 *   E'[rho] = E_DFT[rho] - sum_i lambda_i . (Mi_i - M_target_i)
 *
 * IMPORTANT: Returns 0.0 if magnetic moments are NOT yet converged.
 * This is because the constraint energy is only physically meaningful
 * when Mi ≈ M_target. Before convergence, the lambda values are still
 * adjusting and the energy would be misleading.
 *
 * @par Output meaning
 * - E_scon < 0: lambda and Mi are aligned (system resists the constraint)
 * - E_scon > 0: lambda and Mi are anti-aligned (constraint assists the system)
 * - E_scon = 0: not converged OR all lambda = 0 (no constraint needed)
 *
 * @return Constraint energy in Ry (0.0 if not converged)
 */
template <typename TK>
double SpinConstrain<TK>::cal_escon()
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

template <typename TK>
double SpinConstrain<TK>::get_escon() const
{
    return this->escon_;
}

// set atomCounts
template <typename TK>
void SpinConstrain<TK>::set_atomCounts(const std::map<int, int>& atomCounts_in)
{
    this->atomCounts.clear();
    this->atomCounts = atomCounts_in;
}

// get atomCounts
template <typename TK>
const std::map<int, int>& SpinConstrain<TK>::get_atomCounts() const
{
    return this->atomCounts;
}

/// set nspin
template <typename TK>
void SpinConstrain<TK>::set_nspin(int nspin_in)
{
    if (nspin_in != 4 && nspin_in != 2)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_nspin", "nspin must be 2 or 4");
    }
    this->nspin_ = nspin_in;
}

/// get nspin
template <typename TK>
int SpinConstrain<TK>::get_nspin() const
{
    return this->nspin_;
}

template <typename TK>
void SpinConstrain<TK>::set_npol(int npol)
{
    this->npol_ = npol;
}

template <typename TK>
int SpinConstrain<TK>::get_npol() const
{
    return this->npol_;
}

/**
 * @brief Get spin sign for k-point: determines whether this k-point is
 * spin-up (+1) or spin-down (-1) in collinear (nspin=2) calculations.
 *
 * @details In collinear spin, the wavefunction is split into two k-point pools:
 * - isk[ik] == 0: spin-up channel (majority spin) -> sign = +1
 * - isk[ik] == 1: spin-down channel (minority spin) -> sign = -1
 * For non-collinear (npol=2), always returns +1 since both components
 * are handled together.
 *
 * @return +1 for spin-up, -1 for spin-down, +1 for non-collinear
 */
template <typename TK>
int SpinConstrain<TK>::get_spin_sign(int ik) const
{
    if (this->npol_ == 2) return 1;
    // npol == 1 (nspin == 2): isk[ik]==0 => spin-up (+1), isk[ik]==1 => spin-down (-1)
    return (this->pelec->klist->isk[ik] == 0) ? 1 : -1;
}

/**
 * @brief Accumulate magnetic moments from projector coefficients (becp) for one k-point.
 *
 * @par Algorithm (npol=2, non-collinear):
 * For each atom, compute the 2x2 occupation matrix from becp coefficients:
 *   occ[0] = sum_ih becp_up^*(ih) * becp_up(ih)   = <psi_up|P_at|psi_up>
 *   occ[1] = sum_ih becp_up^*(ih) * becp_dn(ih)   = <psi_up|P_at|psi_dn>
 *   occ[2] = sum_ih becp_dn^*(ih) * becp_up(ih)   = <psi_dn|P_at|psi_up>
 *   occ[3] = sum_ih becp_dn^*(ih) * becp_dn(ih)   = <psi_dn|P_at|psi_dn>
 * where P_at = sum_{l,m} |alpha_{l,m}><alpha_{l,m}| is the atomic projector.
 *
 * The magnetic moment is extracted via Pauli matrix traces:
 *   Mx = Re(occ[1] + occ[2]), My = Im(occ[1] - occ[2]), Mz = Re(occ[0] - occ[3])
 *
 * @par Algorithm (npol=1, collinear):
 * Only the z-component (spin projection) is computed:
 *   occ = sum_ih |becp(ih)|^2 = <psi|P_at|psi>
 *   Mz += weight * occ * spin_sign
 * where spin_sign = +1 for spin-up, -1 for spin-down.
 *
 * @param becp Projector coefficients, layout: [ib * npol * nkb + spin * nkb + ih]
 * @param nkb Total number of projectors across all atoms
 * @param nbands Number of bands (occupied + unoccupied in the subspace)
 * @param npol Number of spinor components (1 for collinear, 2 for non-collinear)
 * @param ik K-point index (used for spin_sign lookup in collinear mode)
 * @param wg_ik Band occupation weights for this k-point (from Fermi-Dirac)
 * @param nh_iat Array of projector counts per atom: nh_iat[iat] = nproj for atom iat
 */
template <typename TK>
void SpinConstrain<TK>::accumulate_Mi_from_becp(const std::complex<double>* becp,
                                                  int nkb,
                                                  int nbands,
                                                  int npol,
                                                  int ik,
                                                  const double* wg_ik,
                                                  const int* nh_iat)
{
    if (npol == 2)
    {
        for (int ib = 0; ib < nbands; ib++)
        {
            const double weight = wg_ik[ib];
            int begin_ih = 0;
            for (int iat = 0; iat < static_cast<int>(this->Mi_.size()); iat++)
            {
                std::complex<double> occ[4] = {ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO};
                const int nh = nh_iat[iat];
                for (int ih = 0; ih < nh; ih++)
                {
                    const int index = ib * 2 * nkb + begin_ih + ih;
                    occ[0] += conj(becp[index]) * becp[index];
                    occ[1] += conj(becp[index]) * becp[index + nkb];
                    occ[2] += conj(becp[index + nkb]) * becp[index];
                    occ[3] += conj(becp[index + nkb]) * becp[index + nkb];
                }
                this->Mi_[iat] += pauli_to_moment(occ, weight);
                begin_ih += nh;
            }
        }
    }
    else // npol == 1
    {
        const int sign = this->get_spin_sign(ik);
        for (int ib = 0; ib < nbands; ib++)
        {
            const double weight = wg_ik[ib];
            int begin_ih = 0;
            for (int iat = 0; iat < static_cast<int>(this->Mi_.size()); iat++)
            {
                double occ = 0.0;
                const int nh = nh_iat[iat];
                for (int ih = 0; ih < nh; ih++)
                {
                    const int index = ib * nkb + begin_ih + ih;
                    occ += (conj(becp[index]) * becp[index]).real();
                }
                this->Mi_[iat].z += weight * occ * sign;
                begin_ih += nh;
            }
        }
    }
}

template <typename TK>
int SpinConstrain<TK>::get_nw() const
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
template <typename TK>
int SpinConstrain<TK>::get_iwt(int itype, int iat, int orbital_index) const
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
template <typename TK>
int SpinConstrain<TK>::get_nat()
{
    int nat = 0;
    for (std::map<int, int>::iterator it = this->atomCounts.begin(); it != this->atomCounts.end(); ++it)
    {
        nat += it->second;
    }
    return nat;
}

/// @brief Get number of element types
template <typename TK>
int SpinConstrain<TK>::get_ntype()
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
template <typename TK>
void SpinConstrain<TK>::check_atomCounts()
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
template <typename TK>
int SpinConstrain<TK>::get_iat(int itype, int atom_index)
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
template <typename TK>
void SpinConstrain<TK>::set_orbitalCounts(const std::map<int, int>& orbitalCounts_in)
{
    this->orbitalCounts.clear();
    this->orbitalCounts = orbitalCounts_in;
}

// get orbitalCounts
template <typename TK>
const std::map<int, int>& SpinConstrain<TK>::get_orbitalCounts() const
{
    return this->orbitalCounts;
}

// set lnchiCounts
template <typename TK>
void SpinConstrain<TK>::set_lnchiCounts(const std::map<int, std::map<int, int>>& lnchiCounts_in)
{
    this->lnchiCounts.clear();
    this->lnchiCounts = lnchiCounts_in;
}

// get lnchiCounts
template <typename TK>
const std::map<int, std::map<int, int>>& SpinConstrain<TK>::get_lnchiCounts() const
{
    return this->lnchiCounts;
}

// set sc_lambda from ScData (parsed from STRU file)
// ScData is organized by element type; this function flattens it to per-atom arrays
template <typename TK>
void SpinConstrain<TK>::set_sc_lambda()
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
template <typename TK>
void SpinConstrain<TK>::set_target_mag()
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
template <typename TK>
void SpinConstrain<TK>::set_constrain()
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
template <typename TK>
void SpinConstrain<TK>::set_sc_lambda(const ModuleBase::Vector3<double>* lambda_in, int nat_in)
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
template <typename TK>
void SpinConstrain<TK>::set_target_mag(const ModuleBase::Vector3<double>* target_mag_in, int nat_in)
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

template <typename TK>
void SpinConstrain<TK>::set_target_mag(const std::vector<ModuleBase::Vector3<double>>& target_mag_in)
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
        ModuleBase::WARNING_QUIT("SpinConstrain::set_target_mag", "nspin must be 2 or 4");
    }
}

/// set constrain from variable
template <typename TK>
void SpinConstrain<TK>::set_constrain(const ModuleBase::Vector3<int>* constrain_in, int nat_in)
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

template <typename TK>
const std::vector<ModuleBase::Vector3<double>>& SpinConstrain<TK>::get_sc_lambda() const
{
    return this->lambda_;
}

template <typename TK>
const std::vector<ModuleBase::Vector3<double>>& SpinConstrain<TK>::get_target_mag() const
{
    return this->target_mag_;
}

/// get_constrain
template <typename TK>
const std::vector<ModuleBase::Vector3<int>>& SpinConstrain<TK>::get_constrain() const
{
    return this->constrain_;
}

/// @brief Reset all atomic magnetic moments to zero. Called before each Mi calculation.
template <typename TK>
void SpinConstrain<TK>::zero_Mi()
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
template <typename TK>
double SpinConstrain<TK>::get_decay_grad(int itype)
{
    return this->ScDecayGrad[itype];
}

/// set grad_decy
template <typename TK>
void SpinConstrain<TK>::set_decay_grad()
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
template <typename TK>
const std::vector<double>& SpinConstrain<TK>::get_decay_grad()
{
    return this->decay_grad_;
}

/// set grad_decy from variable
template <typename TK>
void SpinConstrain<TK>::set_decay_grad(const double* decay_grad_in, int ntype_in)
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
template <typename TK>
void SpinConstrain<TK>::set_input_parameters(double sc_thr_in,
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
template <typename TK>
double SpinConstrain<TK>::get_sc_thr() const
{
    return this->sc_thr_;
}

/// get nsc
template <typename TK>
int SpinConstrain<TK>::get_nsc() const
{
    return this->nsc_;
}

/// get nsc_min
template <typename TK>
int SpinConstrain<TK>::get_nsc_min() const
{
    return this->nsc_min_;
}

/// get alpha_trial
template <typename TK>
double SpinConstrain<TK>::get_alpha_trial() const
{
    return this->alpha_trial_;
}

/// get sccut
template <typename TK>
double SpinConstrain<TK>::get_sccut() const
{
    return this->restrict_current_;
}

/// set sc_drop_thr
template <typename TK>
void SpinConstrain<TK>::set_sc_drop_thr(double sc_drop_thr_in)
{
    this->sc_drop_thr_ = sc_drop_thr_in;
}

/// get sc_drop_thr
template <typename TK>
double SpinConstrain<TK>::get_sc_drop_thr() const
{
    return this->sc_drop_thr_;
}

template <typename TK>
void SpinConstrain<TK>::set_solver_parameters(const K_Vectors& kv_in,
                                                  void* p_hamilt_in,
                                                  void* psi_in,
                                                  elecstate::ElecState* pelec_in)
{
    this->kv_ = kv_in;
    this->p_hamilt = p_hamilt_in;
    this->psi = psi_in;
    this->pelec = pelec_in;
}

/// @brief  set ParaV
template <typename TK>
void SpinConstrain<TK>::set_ParaV(Parallel_Orbitals* ParaV_in)
{
    this->ParaV = ParaV_in;
    int nloc = this->ParaV->nloc;
    if (nloc <= 0)
    {
        ModuleBase::WARNING_QUIT("SpinConstrain::set_ParaV", "nloc <= 0");
    }
}

/**
 * @brief Print magnetic moments per atom in formatted table.
 *
 * @par Output format
 * - nspin=2: "ATOM   1    2.0000000000" (z-component only)
 * - nspin=4: "ATOM   1    0.0010000000    0.0020000000    1.9990000000" (x, y, z)
 *
 * @par Interpretation
 * - Positive Mi.z: spin aligned with z-axis (spin-up character)
 * - Negative Mi.z: spin anti-aligned with z-axis (spin-down character)
 * - Non-zero Mi.x/Mi.y: non-collinear spin components
 * - Mi close to target_mag: constraint is well-satisfied
 * - Mi far from target_mag: constraint is not yet converged
 */
template <typename TK>
void SpinConstrain<TK>::print_Mi(std::ofstream& ofs_running)
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
        FmtTable table(/*titles=*/title, 
                       /*nrows=*/nat, 
                       /*formats=*/fmts, 
                       /*indent=*/0,
                       /*align=*/{/*value*/FmtTable::Align::RIGHT, /*title*/FmtTable::Align::LEFT});
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
        FmtTable table(/*titles=*/title, 
                       /*nrows=*/nat, 
                       /*formats=*/fmts, 
                       /*indent=*/0,
                       /*align=*/{/*value*/FmtTable::Align::RIGHT, /*title*/FmtTable::Align::LEFT});
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

/**
 * @brief Print the magnetic force (-lambda) per atom in eV/uB.
 *
 * @par Physical meaning
 * The "magnetic force" is the derivative of the constrained Lagrangian
 * with respect to the magnetic moment: dL/dMi = -lambda_i.
 * It represents how much energy would change if the constraint were relaxed.
 *
 * @par Interpretation
 * - Large |lambda|: The system strongly resists the target moment constraint
 * - lambda ≈ 0: The system naturally has the target moment (no constraint needed)
 * - Positive lambda.z: The constraint pushes the moment in the +z direction
 * - Negative lambda.z: The constraint pushes the moment in the -z direction
 *
 * @par Typical values
 * - Well-converged SCF: lambda ~ 0.01-1 eV/uB
 * - Strongly constrained: lambda ~ 1-10 eV/uB
 * - Diverging SCF: lambda growing without bound (check target_mag合理性)
 */
template <typename TK>
void SpinConstrain<TK>::print_Mag_Force(std::ofstream& ofs_running)
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
        FmtTable table(/*titles=*/title, 
                       /*nrows=*/nat, 
                       /*formats=*/fmts, 
                       /*indent=*/0,
                       /*align=*/{/*value*/FmtTable::Align::RIGHT, /*title*/FmtTable::Align::LEFT});
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
        FmtTable table(/*titles=*/title, 
                       /*nrows=*/nat, 
                       /*formats=*/fmts, 
                       /*indent=*/0,
                       /*align=*/{/*value*/FmtTable::Align::RIGHT, /*title*/FmtTable::Align::LEFT});
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

/**
 * @brief Reset DeltaSpin operator initialization state.
 *
 * @details The DeltaSpin operator caches internal state (projector matrices, etc.)
 * from a previous SCF iteration. When the constraint parameters change (e.g., new
 * target moments or lambda values), the cached state may be invalid. This function
 * forces the operator to reinitialize on the next call.
 *
 * @par When to call
 * - After changing target_mag_ or constrain_ arrays
 * - When restarting from a previous SCF calculation with different constraints
 * - When switching between LCAO and PW basis sets
 */
template <typename TK>
void SpinConstrain<TK>::reset_dspin_operator()
{
#ifdef __LCAO
    if (this->p_operator == nullptr)
    {
        return;
    }
    if (this->nspin_ == 4)
    {
        auto* dspin = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>*>(this->p_operator);
        if (dspin)
        {
            dspin->reset_initialized();
        }
    }
    else if (this->nspin_ == 2)
    {
        auto* dspin = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator);
        if (dspin)
        {
            dspin->reset_initialized();
        }
    }
#endif
}

template class SpinConstrain<std::complex<double>>;
template class SpinConstrain<double>;

} // namespace spinconstrain
