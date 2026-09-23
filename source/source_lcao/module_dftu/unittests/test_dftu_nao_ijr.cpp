#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <vector>

/***********************************************************************
 * Unit tests for the IJR (atom-pair) helper functions in dftu_nao_ijr.h:
 * cal_coeff_lambda (spin encoding) and accumulate_hr_for_iat0 (HR assembly).
 ***********************************************************************/

// =====================================================================
// 1. cal_coeff_lambda: Lambda coefficient encoding
// Collinear (nspin=2): coeff[0]=lambda_z, coeff[1]=-lambda_z
// Non-collinear (nspin=4): coeff[0]=lambda_z, coeff[1]=lambda_x+i*lambda_y,
//   coeff[2]=lambda_x-i*lambda_y, coeff[3]=-lambda_z
// =====================================================================

static void cal_coeff_lambda_collinear(const std::vector<double>& lambda,
                                        std::vector<double>& coeff)
{ coeff[0] = lambda[0]; coeff[1] = -lambda[0]; }

static void cal_coeff_lambda_noncollinear(const std::vector<double>& lambda,
                                           std::vector<std::complex<double>>& coeff)
{
    coeff[0] = std::complex<double>(lambda[2], 0.0);
    coeff[1] = std::complex<double>(lambda[0], lambda[1]);
    coeff[2] = std::complex<double>(lambda[0], -lambda[1]);
    coeff[3] = std::complex<double>(-lambda[2], 0.0);
}

class CalCoeffLambdaTest : public ::testing::Test { protected: void SetUp() override {} };

TEST_F(CalCoeffLambdaTest, Collinear_PositiveLambdaZ)
{
    std::vector<double> lambda = {2.5}, coeff(2);
    cal_coeff_lambda_collinear(lambda, coeff);
    EXPECT_DOUBLE_EQ(coeff[0], 2.5); EXPECT_DOUBLE_EQ(coeff[1], -2.5);
}

TEST_F(CalCoeffLambdaTest, NonCollinear_General)
{
    std::vector<double> lambda = {1.0, 2.0, 3.0};
    std::vector<std::complex<double>> coeff(4);
    cal_coeff_lambda_noncollinear(lambda, coeff);
    EXPECT_NEAR(coeff[0].real(), 3.0, 1e-15);
    EXPECT_NEAR(coeff[1].real(), 1.0, 1e-15); EXPECT_NEAR(coeff[1].imag(), 2.0, 1e-15);
    EXPECT_NEAR(coeff[2].real(), 1.0, 1e-15); EXPECT_NEAR(coeff[2].imag(), -2.0, 1e-15);
    EXPECT_NEAR(coeff[3].real(), -3.0, 1e-15);
}

// =====================================================================
// 2. accumulate_hr_for_iat0: HR assembly for one Hubbard atom
// =====================================================================

#include "../dftu_nao_ijr.h"

#include <algorithm>
#include <memory>

class AccumulateHrIat0Test : public ::testing::Test
{
  protected:
    void SetUp() override
    {
#ifdef __MPI
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
        ucell.ntype = 1;
        ucell.nat = nat;
        ucell.atoms = atoms_buf;
        iat2it_buf.assign(nat, 0);
        iat2ia_buf.resize(nat);
        ucell.iat2it = iat2it_buf.data();
        ucell.iat2ia = iat2ia_buf.data();
        ucell.atoms[0].tau.resize(nat);
        ucell.lat0 = 1.0;
        ucell.itia2iat.create(ucell.ntype, nat);
        for (int iat = 0; iat < nat; iat++)
        {
            iat2ia_buf[iat] = iat;
            ucell.atoms[0].tau[iat] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
            ucell.itia2iat(0, iat) = iat;
        }
        ucell.atoms[0].na = nat;
        ucell.atoms[0].nw = nw;
        ucell.atoms[0].iw2l = {2, 2};
        ucell.atoms[0].iw2m = {0, 0};
        ucell.atoms[0].iw2n = {0, 0};
        ucell.set_iat2iwt(1);
        paraV.reset(new Parallel_Orbitals());
#ifdef __MPI
        paraV->init(nat * nw, nat * nw, nat * nw, MPI_COMM_WORLD);
        paraV->set_atomic_trace(ucell.get_iat2iwt(), nat, nat * nw);
#endif
        HR.reset(new hamilt::HContainer<double>(ucell, paraV.get()));
        std::fill(HR->get_wrapper(), HR->get_wrapper() + HR->get_nnr(), 0.0);
    }

    void TearDown() override
    {
        HR.reset();
        paraV.reset();
        ucell.atoms = nullptr;
        ucell.iat2it = nullptr;
        ucell.iat2ia = nullptr;
    }

    AdjacentAtomInfo make_adjs() const
    {
        AdjacentAtomInfo adjs;
        adjs.adj_num = 1;
        adjs.ntype = {0, 0};
        adjs.natom = {0, 1};
        adjs.box = {ModuleBase::Vector3<int>(0, 0, 0), ModuleBase::Vector3<int>(0, 0, 0)};
        return adjs;
    }

    DFTU_LCAO::NlmTot make_nlm_tot() const
    {
        DFTU_LCAO::NlmTot nlm_tot(nat);
        for (int iat = 0; iat < nat; ++iat)
        {
            nlm_tot[iat].resize(2);
            for (int ad = 0; ad < 2; ++ad)
            {
                for (int iw = 0; iw < nw; ++iw)
                {
                    for (int m = 0; m < 5; ++m)
                    {
                        nlm_tot[iat][ad][iw * 5 + m] = std::vector<double>(5, 1.0);
                    }
                }
            }
        }
        return nlm_tot;
    }

    const int nat = 2;
    const int nw = 2;
    int dsize = 1;
    int my_rank = 0;
    UnitCell ucell;
    Atom atoms_buf[1];
    std::vector<int> iat2it_buf;
    std::vector<int> iat2ia_buf;
    std::unique_ptr<Parallel_Orbitals> paraV;
    std::unique_ptr<hamilt::HContainer<double>> HR;
};

TEST_F(AccumulateHrIat0Test, AccumulatesAllPairs)
{
    AdjacentAtomInfo adjs = make_adjs();
    DFTU_LCAO::NlmTot nlm_tot = make_nlm_tot();
    std::vector<double> pot_onsite(25, 0.0);
    for (int m = 0; m < 5; ++m)
    {
        pot_onsite[m * 5 + m] = 1.0;
    }

    DFTU_LCAO::accumulate_hr_for_iat0<double>(ucell, HR.get(), nlm_tot, 0, adjs, *paraV, pot_onsite);

    for (int iap = 0; iap < HR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = HR->get_atom_pair(iap);
        const int iat1 = tmp.get_atom_i();
        const int iat2 = tmp.get_atom_j();
        std::vector<int> indexes1 = paraV->get_indexes_row(iat1);
        std::vector<int> indexes2 = paraV->get_indexes_col(iat2);
        const int nwt = indexes1.size() * indexes2.size();
        for (int i = 0; i < nwt; ++i)
        {
            EXPECT_NEAR(tmp.get_pointer(0)[i], 5.0, 1e-12);
        }
    }
}

TEST_F(AccumulateHrIat0Test, MissingPairIsSkipped)
{
    AdjacentAtomInfo adjs = make_adjs();
    adjs.box[1] = ModuleBase::Vector3<int>(5, 5, 5);
    ASSERT_EQ(HR->find_matrix(0, 1, 5, 5, 5), nullptr);
    DFTU_LCAO::NlmTot nlm_tot = make_nlm_tot();
    std::vector<double> pot_onsite(25, 1.0);

    DFTU_LCAO::accumulate_hr_for_iat0<double>(ucell, HR.get(), nlm_tot, 0, adjs, *paraV, pot_onsite);

    for (int iap = 0; iap < HR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = HR->get_atom_pair(iap);
        const int iat1 = tmp.get_atom_i();
        const int iat2 = tmp.get_atom_j();
        std::vector<int> indexes1 = paraV->get_indexes_row(iat1);
        std::vector<int> indexes2 = paraV->get_indexes_col(iat2);
        const int nwt = indexes1.size() * indexes2.size();
        const double expected = (iat1 == iat2) ? 25.0 : 0.0;
        for (int i = 0; i < nwt; ++i)
        {
            EXPECT_DOUBLE_EQ(tmp.get_pointer(0)[i], expected);
        }
    }
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
