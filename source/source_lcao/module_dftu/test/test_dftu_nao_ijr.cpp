#include "../dftu_nao_ijr.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <memory>
#include <vector>

/// @file test_dftu_nao_ijr.cpp
/// @brief Focused unit test for the free function
///        DFTU_LCAO::accumulate_hr_for_iat0 (dftu_nao_ijr.h), which
///        accumulates the real-space HR blocks of one Hubbard atom.
///        Two atoms with 2 d-type orbitals each (nw=2, iw2l=2); the
///        Hubbard center (iat0=0) sees both atoms as neighbors at R=0.
///        With nlm=1, pot_onsite(m,m')=delta_{m,m'} and the full (I,J,R)
///        pair grid present in HR, every HR entry accumulates 5.
///        iat2it/iat2ia are borrowed pointers into iat2it_buf/iat2ia_buf;
///        their lifetime is managed by the fixture vectors (and ultimately
///        released by UnitCell's Statistics destructor).
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
        // HR holds the full 3x3x3 R grid of atom pairs
        HR.reset(new hamilt::HContainer<double>(ucell, paraV.get()));
        std::fill(HR->get_wrapper(), HR->get_wrapper() + HR->get_nnr(), 0.0);
    }

    void TearDown() override
    {
        // reset HR/paraV before unitcell buffers go out of scope;
        // iat2it/iat2ia are released by UnitCell's Statistics member
        HR.reset();
        paraV.reset();
        ucell.atoms = nullptr;
        ucell.iat2it = nullptr;
        ucell.iat2ia = nullptr;
    }

    // neighbor list of Hubbard center iat0=0: both atoms at box (0,0,0)
    AdjacentAtomInfo make_adjs() const
    {
        AdjacentAtomInfo adjs;
        adjs.adj_num = 1; // one neighbor besides the center itself
        adjs.ntype = {0, 0};
        adjs.natom = {0, 1};
        adjs.box = {ModuleBase::Vector3<int>(0, 0, 0), ModuleBase::Vector3<int>(0, 0, 0)};
        return adjs;
    }

    // nlm_tot[iat0][ad][orbital-index(iw*5+m)] = 1 for both neighbors
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
    Atom atoms_buf[1]; // borrowed by ucell.atoms (raw Atom* member)
    std::vector<int> iat2it_buf; // borrowed by ucell.iat2it
    std::vector<int> iat2ia_buf; // borrowed by ucell.iat2ia
    std::unique_ptr<Parallel_Orbitals> paraV;
    std::unique_ptr<hamilt::HContainer<double>> HR;
};

TEST_F(AccumulateHrIat0Test, AccumulatesAllPairs)
{
    AdjacentAtomInfo adjs = make_adjs();
    DFTU_LCAO::NlmTot nlm_tot = make_nlm_tot();
    // 5x5 identity: pot_onsite(m,m') = delta_{m,m'}
    std::vector<double> pot_onsite(25, 0.0);
    for (int m = 0; m < 5; ++m)
    {
        pot_onsite[m * 5 + m] = 1.0;
    }

    DFTU_LCAO::accumulate_hr_for_iat0<double>(ucell, HR.get(), nlm_tot, 0, adjs, *paraV, pot_onsite);

    // every local HR entry gets sum_m 1*1*1 = 5 per orbital pair,
    // i.e. each element of an orbital-pair block equals 5;
    // each atom pair holds 2x2 orbital pairs => each HR value is 5
    for (int iap = 0; iap < HR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = HR->get_atom_pair(iap);
        const int iat1 = tmp.get_atom_i();
        const int iat2 = tmp.get_atom_j();
        std::vector<int> indexes1 = paraV->get_indexes_row(iat1);
        std::vector<int> indexes2 = paraV->get_indexes_col(iat2);
        const int nwt = indexes1.size() * indexes2.size();
        // only pairs with R = (0,0,0) are touched; the R grid of
        // HContainer(ucell) only contains the R=0 block per pair
        for (int i = 0; i < nwt; ++i)
        {
            EXPECT_NEAR(tmp.get_pointer(0)[i], 5.0, 1e-12);
        }
    }
}

TEST_F(AccumulateHrIat0Test, MissingPairIsSkipped)
{
    AdjacentAtomInfo adjs = make_adjs();
    // move neighbor ad2 to R=(5,5,5): the cross pairs (ad1,ad2)=(0,1)
    // and (1,0) then map to R_vector = +/-(5,5,5), for which HR holds
    // no block, so accumulate_hr_for_iat0 must skip them. The diagonal
    // pairs (0,0) and (1,1) still have R_vector = 0 and are accumulated.
    adjs.box[1] = ModuleBase::Vector3<int>(5, 5, 5);
    // precondition: HR really has no (0,1) block at R=(5,5,5)
    ASSERT_EQ(HR->find_matrix(0, 1, 5, 5, 5), nullptr);
    DFTU_LCAO::NlmTot nlm_tot = make_nlm_tot();
    // pot_onsite(m,m') = 1 everywhere: each touched entry gets
    // sum_{m,m'} 1*1*1 = 25
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
        // diagonal pairs are accumulated (25); cross pairs whose only
        // contribution sits at the missing R block stay zero
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
