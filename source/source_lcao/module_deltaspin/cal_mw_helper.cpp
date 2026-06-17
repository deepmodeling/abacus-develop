#ifdef __LCAO
#include "spin_constrain.h"

/**
 * @file cal_mw_helper.cpp
 * @brief LCAO-specific helper functions for magnetic moment calculation from orbital matrices.
 *
 * @par Purpose
 * Provides alternative paths for computing magnetic moments from the orbital
 * multiplication matrix (orbMulP) and the mu*density matrix (mud). These are
 * used when the DeltaSpin operator path is not available or for debugging.
 *
 * @par Data flow
 * 1. convert(): Flatten orbMulP into nested vector [nspin][iat][iw]
 * 2. calculate_MW(): Sum orbital contributions per atom, compute Mi
 * 3. collect_MW(): Accumulate mu*dm contributions into MecMulP matrix
 */

/**
 * @brief Convert flat orbital matrix to nested vector format.
 *
 * @details The orbMulP matrix stores orbital contributions in a flat layout:
 *   orbMulP(is, num) where num runs through all orbitals of all atoms.
 * This function reorganizes it into a nested structure:
 *   AorbMulP[is][iat][iw] = orbMulP(is, num)
 *
 * Values below 1e-10 are set to 0.0 to avoid floating-point noise.
 *
 * @param orbMulP Flat matrix of orbital contributions [nspin x ntotal_orbitals]
 * @return Nested vector [nspin][iat][iw]
 */
template <>
std::vector<std::vector<std::vector<double>>> spinconstrain::SpinConstrain<std::complex<double>>::convert(
    const ModuleBase::matrix& orbMulP)
{
    std::vector<std::vector<std::vector<double>>> AorbMulP;
    AorbMulP.resize(this->nspin_);
    int nat = this->get_nat();
    for (int is = 0; is < this->nspin_; ++is)
    {
        int num = 0;
        AorbMulP[is].resize(nat);
        for (const auto& sc_elem: this->get_atomCounts())
        {
            int it = sc_elem.first;
            int nat_it = sc_elem.second;
            int nw_it = this->get_orbitalCounts().at(it);
            for (int ia = 0; ia < nat_it; ia++)
            {
                int iat = this->get_iat(it, ia);
                AorbMulP[is][iat].resize(nw_it, 0.0);
                for (int iw = 0; iw < nw_it; iw++)
                {
                    AorbMulP[is][iat][iw] = std::abs(orbMulP(is, num))< 1e-10 ? 0.0 : orbMulP(is, num);
                    num++;
                }
            }
        }
    }
    return AorbMulP;
}

/**
 * @brief Calculate magnetic moments from converted orbital matrix.
 *
 * @par Algorithm (nspin=2):
 *   atom_mag = sum(orbMulP[0][iat]) - sum(orbMulP[1][iat])
 *   Mi[iat].z = atom_mag (z-component only)
 *
 * @par Algorithm (nspin=4):
 * The 4 spinor components are mapped to magnetic moments:
 *   total_charge_soc[0] = Tr(rho * I) / 2      (charge)
 *   total_charge_soc[1] = Tr(rho * sigma_x)    (Mx)
 *   total_charge_soc[2] = Tr(rho * sigma_y)    (My)
 *   total_charge_soc[3] = Tr(rho * sigma_z)    (Mz)
 * Components below sc_thr_ are set to 0.0 to avoid noise.
 *
 * @param AorbMulP Nested vector [nspin][iat][iw] from convert()
 */
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::calculate_MW(
    const std::vector<std::vector<std::vector<double>>>& AorbMulP)
{
    size_t nw = this->get_nw();
    int nat = this->get_nat();

    this->zero_Mi();

    const int nlocal = (this->nspin_ == 4) ? nw / 2 : nw;
    for (const auto& sc_elem: this->get_atomCounts())
    {
        int it = sc_elem.first;
        int nat_it = sc_elem.second;
        for (int ia = 0; ia < nat_it; ia++)
        {
            int num = 0;
            int iat = this->get_iat(it, ia);
            double atom_mag = 0.0;
            std::vector<double> total_charge_soc(this->nspin_, 0.0);
            for (const auto& lnchi: this->get_lnchiCounts().at(it))
            {
                std::vector<double> sum_l(this->nspin_, 0.0);
                int L = lnchi.first;
                int nchi = lnchi.second;
                for (int Z = 0; Z < nchi; ++Z)
                {
                    std::vector<double> sum_m(this->nspin_, 0.0);
                    for (int M = 0; M < (2 * L + 1); ++M)
                    {
                        for (int j = 0; j < this->nspin_; j++)
                        {
                            sum_m[j] += AorbMulP[j][iat][num];
                        }
                        num++;
                    }
                    for (int j = 0; j < this->nspin_; j++)
                    {
                        sum_l[j] += sum_m[j];
                    }
                }
                if (this->nspin_ == 2)
                {
                    atom_mag += sum_l[0] - sum_l[1];
                }
                else if (this->nspin_ == 4)
                {
                    for (int j = 0; j < this->nspin_; j++)
                    {
                        total_charge_soc[j] += sum_l[j];
                    }
                }
            }
            if (this->nspin_ == 2)
            {
                this->Mi_[iat].x = 0.0;
                this->Mi_[iat].y = 0.0;
                this->Mi_[iat].z = atom_mag;
            }
            else if (this->nspin_ == 4)
            {
                this->Mi_[iat].x = (std::abs(total_charge_soc[1]) < this->sc_thr_)? 0.0 : total_charge_soc[1];
                this->Mi_[iat].y = (std::abs(total_charge_soc[2]) < this->sc_thr_)? 0.0 : total_charge_soc[2];
                this->Mi_[iat].z = (std::abs(total_charge_soc[3]) < this->sc_thr_)? 0.0 : total_charge_soc[3];
            }
        }
    }
}

/**
 * @brief Accumulate magnetic moment contributions from mu*density matrix.
 *
 * @details For distributed matrices (ScaLAPACK), only the local processor's
 * elements are accumulated. The ParaV mapping converts global indices to
 * local row/column indices.
 *
 * @par nspin=4 spinor decomposition
 * The mud matrix stores the 2x2 spinor blocks interleaved:
 *   Global index 2j -> spin-up component
 *   Global index 2j+1 -> spin-down component
 * The Pauli matrix traces are:
 *   M0 (charge): mud(k1,k1).real + mud(k2,k2).real
 *   M3 (Mz):     mud(k1,k1).real - mud(k2,k2).real
 *   M1 (Mx):     mud(k1,k2).real + mud(k2,k1).real
 *   M2 (My):    -mud(k1,k2).imag + mud(k2,k1).imag
 *
 * @param MecMulP Output matrix [4 x nw/2]: MecMulP[0]=charge, [1]=Mx, [2]=My, [3]=Mz
 * @param mud Input mu*density matrix (column-major)
 * @param nw Total number of orbitals
 * @param isk Spin index (0 or 1 for nspin=2)
 */
template <>
void spinconstrain::SpinConstrain<std::complex<double>>::collect_MW(ModuleBase::matrix& MecMulP,
                                                      const ModuleBase::ComplexMatrix& mud,
                                                      int nw,
                                                      int isk)
{
    if (this->nspin_ == 2)
    {
        for (size_t i=0; i < nw; ++i)
        {
            if (this->ParaV->in_this_processor(i, i))
            {
                const int ir = this->ParaV->global2local_row(i);
                const int ic = this->ParaV->global2local_col(i);
                MecMulP(isk, i) += mud(ic, ir).real();
            }
        }
    }
    else if (this->nspin_ == 4)
    {
        for (size_t i = 0; i < nw; ++i)
        {
            const int index = i % 2;
            if (!index)
            {
                const int j = i / 2;
                const int k1 = 2 * j;
                const int k2 = 2 * j + 1;
                if (this->ParaV->in_this_processor(k1, k1))
                {
                    const int ir = this->ParaV->global2local_row(k1);
                    const int ic = this->ParaV->global2local_col(k1);
                    MecMulP(0, j) += mud(ic, ir).real();
                    MecMulP(3, j) += mud(ic, ir).real();
                }
                if (this->ParaV->in_this_processor(k1, k2))
                {
                    const int ir = this->ParaV->global2local_row(k1);
                    const int ic = this->ParaV->global2local_col(k2);
                    // note that mud is column major
                    MecMulP(1, j) += mud(ic, ir).real();
                    // M_y = i(M_{up,down} - M_{down,up}) = -(M_{up,down} - M_{down,up}).imag()
                    MecMulP(2, j) -= mud(ic, ir).imag();
                }
                if (this->ParaV->in_this_processor(k2, k1))
                {
                    const int ir = this->ParaV->global2local_row(k2);
                    const int ic = this->ParaV->global2local_col(k1);
                    MecMulP(1, j) += mud(ic, ir).real();
                    // M_y = i(M_{up,down} - M_{down,up}) = -(M_{up,down} - M_{down,up}).imag()
                    MecMulP(2, j) += mud(ic, ir).imag();
                }
                if (this->ParaV->in_this_processor(k2, k2))
                {
                    const int ir = this->ParaV->global2local_row(k2);
                    const int ic = this->ParaV->global2local_col(k2);
                    MecMulP(0, j) += mud(ic, ir).real();
                    MecMulP(3, j) -= mud(ic, ir).real();
                }
            }
        }
    }
}

#endif
