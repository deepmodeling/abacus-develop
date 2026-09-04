/**
 * @file deltaspin_lcao_mi.cpp
 * @brief LCAO-specific magnetic-moment calculation for DeltaSpin.
 *
 * @par Calculation methods
 * - cal_mi_lcao(): Uses the DeltaSpin operator to compute Tr(rho * mu) directly
 *   from the density matrix. This is the primary method.
 * - convert_orbital_matrix()/calculate_mw_from_orbitals(): Alternative path via the
 *   orbital multiplication matrix, mainly for debugging.
 * - collect_mw(): Accumulates mu*density-matrix contributions for distributed
 *   (ScaLAPACK) matrices.
 *
 * @par nspin=2 (collinear)
 * Density matrix dmr[0] is spin-up, dmr[1] is spin-down. Mi is the
 * difference (up - down), giving only the z-component.
 *
 * @par nspin=4 (non-collinear)
 * Density matrix has 4 interleaved spinor components. The DeltaSpin
 * operator decomposes them into charge and 3 magnetic components via
 * Pauli matrix traces.
 */
#ifdef __LCAO

#include "deltaspin_lcao_mi.h"

#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_base/timer.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/module_operator_lcao/dspin_lcao.h"

#include <cmath>

namespace spinconstrain
{
namespace lcao
{

void cal_mi_lcao(ScState& state,
                 hamilt::Operator<std::complex<double>>* p_operator,
                 elecstate::DensityMatrix<std::complex<double>, double>* dm,
                 const int& step,
                 bool print)
{
    ModuleBase::TITLE("module_deltaspin", "cal_mi_lcao");
    ModuleBase::timer::start("spinconstrain::SpinConstrain", "cal_mi_lcao");
    // Reset Mi before calculation
    state.zero_Mi();
    const hamilt::HContainer<double>* dmr = dm->get_DMR_pointer(1);
    std::vector<double> moments;
    if (state.nspin_ == 2)
    {
        // Switch to spin-difference density matrix (rho_up - rho_dn)
        dm->switch_dmr(2);

        // Compute moments via DeltaSpin operator
        moments = static_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(p_operator)->cal_moment(dmr, state.get_constrain());

        // Switch back to total density matrix
        dm->switch_dmr(0);

        // For nspin=2, only z-component is meaningful
        for (int iat = 0; iat < state.Mi_.size(); iat++)
        {
            state.Mi_[iat].x = 0.0;
            state.Mi_[iat].y = 0.0;
            state.Mi_[iat].z = moments[iat];
        }
    }
    else if (state.nspin_ == 4)
    {
        // For nspin=4, moments array contains interleaved [Mx, My, Mz] per atom
        moments = static_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>*>(p_operator)->cal_moment(dmr, state.get_constrain());
        for (int iat = 0; iat < state.Mi_.size(); iat++)
        {
            state.Mi_[iat].x = moments[iat * 3];
            state.Mi_[iat].y = moments[iat * 3 + 1];
            state.Mi_[iat].z = moments[iat * 3 + 2];
        }
    }

    ModuleBase::timer::end("spinconstrain::SpinConstrain", "cal_mi_lcao");
}

std::vector<std::vector<std::vector<double>>> convert_orbital_matrix(
    const ModuleBase::matrix& orbMulP,
    const ScState& state)
{
    std::vector<std::vector<std::vector<double>>> AorbMulP;
    AorbMulP.resize(state.nspin_);
    int nat = state.get_nat();
    for (int is = 0; is < state.nspin_; ++is)
    {
        int num = 0;
        AorbMulP[is].resize(nat);
        for (const auto& sc_elem: state.get_atomCounts())
        {
            int it = sc_elem.first;
            int nat_it = sc_elem.second;
            int nw_it = state.get_orbitalCounts().at(it);
            for (int ia = 0; ia < nat_it; ia++)
            {
                int iat = state.get_iat(it, ia);
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

void calculate_mw_from_orbitals(const std::vector<std::vector<std::vector<double>>>& AorbMulP,
                                ScState& state)
{
    size_t nw = state.get_nw();
    int nat = state.get_nat();

    state.zero_Mi();

    for (const auto& sc_elem: state.get_atomCounts())
    {
        int it = sc_elem.first;
        int nat_it = sc_elem.second;
        for (int ia = 0; ia < nat_it; ia++)
        {
            int num = 0;
            int iat = state.get_iat(it, ia);
            double atom_mag = 0.0;
            std::vector<double> total_charge_soc(state.nspin_, 0.0);
            for (const auto& lnchi: state.get_lnchiCounts().at(it))
            {
                std::vector<double> sum_l(state.nspin_, 0.0);
                int L = lnchi.first;
                int nchi = lnchi.second;
                for (int Z = 0; Z < nchi; ++Z)
                {
                    std::vector<double> sum_m(state.nspin_, 0.0);
                    for (int M = 0; M < (2 * L + 1); ++M)
                    {
                        for (int j = 0; j < state.nspin_; j++)
                        {
                            sum_m[j] += AorbMulP[j][iat][num];
                        }
                        num++;
                    }
                    for (int j = 0; j < state.nspin_; j++)
                    {
                        sum_l[j] += sum_m[j];
                    }
                }
                if (state.nspin_ == 2)
                {
                    atom_mag += sum_l[0] - sum_l[1];
                }
                else if (state.nspin_ == 4)
                {
                    for (int j = 0; j < state.nspin_; j++)
                    {
                        total_charge_soc[j] += sum_l[j];
                    }
                }
            }
            if (state.nspin_ == 2)
            {
                state.Mi_[iat].x = 0.0;
                state.Mi_[iat].y = 0.0;
                state.Mi_[iat].z = atom_mag;
            }
            else if (state.nspin_ == 4)
            {
                state.Mi_[iat].x = (std::abs(total_charge_soc[1]) < state.sc_thr_)? 0.0 : total_charge_soc[1];
                state.Mi_[iat].y = (std::abs(total_charge_soc[2]) < state.sc_thr_)? 0.0 : total_charge_soc[2];
                state.Mi_[iat].z = (std::abs(total_charge_soc[3]) < state.sc_thr_)? 0.0 : total_charge_soc[3];
            }
        }
    }
}

void collect_mw(ModuleBase::matrix& MecMulP,
                const ModuleBase::ComplexMatrix& mud,
                int nw,
                int isk,
                const ScState& state,
                const Parallel_Orbitals* pv)
{
    if (state.nspin_ == 2)
    {
        for (size_t i=0; i < nw; ++i)
        {
            if (pv->in_this_processor(i, i))
            {
                const int ir = pv->global2local_row(i);
                const int ic = pv->global2local_col(i);
                MecMulP(isk, i) += mud(ic, ir).real();
            }
        }
    }
    else if (state.nspin_ == 4)
    {
        for (size_t i = 0; i < nw; ++i)
        {
            const int index = i % 2;
            if (!index)
            {
                const int j = i / 2;
                const int k1 = 2 * j;
                const int k2 = 2 * j + 1;
                if (pv->in_this_processor(k1, k1))
                {
                    const int ir = pv->global2local_row(k1);
                    const int ic = pv->global2local_col(k1);
                    MecMulP(0, j) += mud(ic, ir).real();
                    MecMulP(3, j) += mud(ic, ir).real();
                }
                if (pv->in_this_processor(k1, k2))
                {
                    const int ir = pv->global2local_row(k1);
                    const int ic = pv->global2local_col(k2);
                    // note that mud is column major
                    MecMulP(1, j) += mud(ic, ir).real();
                    // M_y = i(M_{up,down} - M_{down,up}) = -(M_{up,down} - M_{down,up}).imag()
                    MecMulP(2, j) -= mud(ic, ir).imag();
                }
                if (pv->in_this_processor(k2, k1))
                {
                    const int ir = pv->global2local_row(k2);
                    const int ic = pv->global2local_col(k1);
                    MecMulP(1, j) += mud(ic, ir).real();
                    // M_y = i(M_{up,down} - M_{down,up}) = -(M_{up,down} - M_{down,up}).imag()
                    MecMulP(2, j) += mud(ic, ir).imag();
                }
                if (pv->in_this_processor(k2, k2))
                {
                    const int ir = pv->global2local_row(k2);
                    const int ic = pv->global2local_col(k2);
                    MecMulP(0, j) += mud(ic, ir).real();
                    MecMulP(3, j) -= mud(ic, ir).real();
                }
            }
        }
    }
}

} // namespace lcao
} // namespace spinconstrain

#endif // __LCAO
