#include "source_lcao/module_extrap/wf_history_lcao.h"
#include "source_lcao/module_extrap/wf_orthonormalize_lcao.h"
#include "source_lcao/module_extrap/wf_snapshot_lcao.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace
{

constexpr int nstate = 1;
constexpr int nbands = 3;
constexpr int nbasis = 3;

void initialize_serial_orbitals(Parallel_Orbitals& pv)
{
    pv.set_serial(nbasis, nbasis);
    pv.ncol_bands = nbands;
    pv.nrow_bands = nbasis;
    pv.nbands = nbands;
#ifdef __MPI
    // get_wfc_global_nbands() reads desc_wfc in MPI-enabled builds even when
    // this test deliberately uses the serial distribution.
    pv.desc_wfc[2] = nbasis;
    pv.desc_wfc[3] = nbands;
#endif
}

psi::Psi<double> make_psi(const std::vector<double>& coeff)
{
    psi::Psi<double> wfc(nstate, nbands, nbasis, nbasis, true);
    wfc.set_all_psi(coeff.data(), coeff.size());
    return wfc;
}

ModuleBase::matrix make_occupations()
{
    ModuleBase::matrix occupations(nstate, 2);
    occupations(0, 0) = 1.0;
    occupations(0, 1) = 1.0;
    return occupations;
}

double metric_element(const psi::Psi<double>& wfc,
                      const std::vector<double>& overlap,
                      const int ib,
                      const int jb)
{
    double value = 0.0;
    for (int mu = 0; mu < nbasis; ++mu)
    {
        for (int nu = 0; nu < nbasis; ++nu)
        {
            value += wfc(0, ib, mu) * overlap[mu * nbasis + nu] * wfc(0, jb, nu);
        }
    }
    return value;
}

const std::vector<double> overlap = {
    2.0, 0.2, 0.0,
    0.2, 1.5, 0.1,
    0.0, 0.1, 1.0,
};

const std::vector<double> coeff = {
    1.0, 0.2, 0.1,
    0.3, 1.0, 0.2,
    7.0, 8.0, 9.0,
};

} // namespace

TEST(WfSnapshotLCAO, OwnsAndRestoresData)
{
    psi::Psi<double> wfc = make_psi(coeff);
    ModuleBase::matrix occupations = make_occupations();
    ModuleExtrap::WfSnapshotLCAO<double> snapshot(7, wfc, occupations);

    wfc(0, 0, 0) = -100.0;
    occupations(0, 0) = 0.0;

    EXPECT_EQ(snapshot.istep, 7);
    EXPECT_DOUBLE_EQ(snapshot.coeff.front(), coeff.front());
    EXPECT_DOUBLE_EQ(snapshot.wg(0, 0), 1.0);

    psi::Psi<double> restored(nstate, nbands, nbasis, nbasis, true);
    ASSERT_TRUE(snapshot.load_to(restored));
    for (std::size_t i = 0; i < coeff.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(restored.get_pointer()[i], coeff[i]);
    }
}

TEST(WfOrthonormalizeLCAO, ReorthonormalizesOccupiedSubspace)
{
    Parallel_Orbitals pv;
    initialize_serial_orbitals(pv);
    psi::Psi<double> wfc = make_psi(coeff);
    const ModuleBase::matrix occupations = make_occupations();

    const ModuleExtrap::WfOrthonormalizeResult result
        = ModuleExtrap::reorthonormalize_gamma_lcao(overlap.data(), pv, wfc, occupations);

    ASSERT_TRUE(result.ok());
    EXPECT_EQ(result.nstate, nstate);
    EXPECT_EQ(result.nbasis, nbasis);
    EXPECT_EQ(result.nactive_bands, 2);
    EXPECT_LT(result.max_deviation, 1.0e-10);

    EXPECT_NEAR(metric_element(wfc, overlap, 0, 0), 1.0, 1.0e-10);
    EXPECT_NEAR(metric_element(wfc, overlap, 1, 1), 1.0, 1.0e-10);
    EXPECT_NEAR(metric_element(wfc, overlap, 0, 1), 0.0, 1.0e-10);

    // The third band has no occupation entry and must not be modified.
    EXPECT_DOUBLE_EQ(wfc(0, 2, 0), 7.0);
    EXPECT_DOUBLE_EQ(wfc(0, 2, 1), 8.0);
    EXPECT_DOUBLE_EQ(wfc(0, 2, 2), 9.0);
}

TEST(WfOrthonormalizeLCAO, DoesNotModifyPsiOnFailure)
{
    Parallel_Orbitals pv;
    initialize_serial_orbitals(pv);
    const std::vector<double> rank_deficient_coeff = {
        1.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        7.0, 8.0, 9.0,
    };
    psi::Psi<double> wfc = make_psi(rank_deficient_coeff);
    const ModuleBase::matrix occupations = make_occupations();

    const ModuleExtrap::WfOrthonormalizeResult result
        = ModuleExtrap::reorthonormalize_gamma_lcao(overlap.data(), pv, wfc, occupations);

    EXPECT_EQ(result.status, ModuleExtrap::WfcExtrapStatus::OrthogonalizationFailed);
    for (std::size_t i = 0; i < rank_deficient_coeff.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(wfc.get_pointer()[i], rank_deficient_coeff[i]);
    }
}

TEST(WfHistoryLCAO, UsesLatestOwnedSnapshot)
{
    Parallel_Orbitals pv;
    initialize_serial_orbitals(pv);
    const ModuleBase::matrix occupations = make_occupations();
    ModuleExtrap::WfHistoryLCAO<double> history(ModuleExtrap::WfcExtrapMethod::UsePrevWf, 2);

    psi::Psi<double> first = make_psi(coeff);
    history.update_after_scf(3, first, occupations);

    std::vector<double> second_coeff = coeff;
    second_coeff[0] = 0.8;
    second_coeff[1] = 0.4;
    psi::Psi<double> second = make_psi(second_coeff);
    history.update_after_scf(4, second, occupations);

    // Mutating the source after insertion must not alter the stored snapshot.
    second.zero_out();

    psi::Psi<double> predicted(nstate, nbands, nbasis, nbasis, true);
    const ModuleExtrap::WfExtrapApplyResult result
        = history.try_use_prev_wf_gamma(overlap.data(), pv, predicted, occupations);

    ASSERT_TRUE(result.ok());
    EXPECT_EQ(result.snapshot_istep, 4);
    EXPECT_EQ(history.size(), 2U);
    EXPECT_NEAR(metric_element(predicted, overlap, 0, 0), 1.0, 1.0e-10);
    EXPECT_NEAR(metric_element(predicted, overlap, 1, 1), 1.0, 1.0e-10);
    EXPECT_NEAR(metric_element(predicted, overlap, 0, 1), 0.0, 1.0e-10);

    history.set_max_depth(1);
    EXPECT_EQ(history.size(), 1U);
}
