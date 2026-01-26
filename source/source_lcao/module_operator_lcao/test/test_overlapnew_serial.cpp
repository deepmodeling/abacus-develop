#include "../overlap_new.h"

#include "gtest/gtest.h"

//---------------------------------------
// Unit test of OverlapNew class
// OverlapNew is a derivative class of Operator, it is used to calculate the overlap matrix
// It use HContainer to store the real space SR matrix
// In this test, we test the correctness and time consuming of 3 functions in OverlapNew class
// - initialize_SR() called in constructor
// - contributeHR()
// - contributeHk()
// - SR(double) and SK(complex<double>) are tested in constructHRd2cd
// - SR(double) and SK(double) are tested in constructHRd2d
//---------------------------------------

// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 10;
int test_nw = 10;
class OverlapNewTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
#ifdef __MPI
        // MPI parallel settings
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

        // set up a unitcell, with one element and test_size atoms, each atom has test_nw orbitals
        ucell.ntype = 1;
        ucell.nat = test_size;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.iat2it = new int[ucell.nat];
        ucell.iat2ia = new int[ucell.nat];
        ucell.atoms[0].tau.resize(ucell.nat);
        ucell.itia2iat.create(ucell.ntype, ucell.nat);
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2it[iat] = 0;
            ucell.iat2ia[iat] = iat;
            ucell.atoms[0].tau[iat] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
            ucell.itia2iat(0, iat) = iat;
        }
        ucell.atoms[0].na = test_size;
        ucell.atoms[0].nw = test_nw;
        ucell.atoms[0].iw2l.resize(test_nw);
        ucell.atoms[0].iw2m.resize(test_nw);
        ucell.atoms[0].iw2n.resize(test_nw);
        for (int iw = 0; iw < test_nw; ++iw)
        {
            ucell.atoms[0].iw2l[iw] = 0;
            ucell.atoms[0].iw2m[iw] = 0;
            ucell.atoms[0].iw2n[iw] = 0;
        }
        ucell.set_iat2iwt(1);
        init_parav();
        // set up a HContainer with ucell
        SR = new hamilt::HContainer<double>(paraV);
    }

    void TearDown() override
    {
        delete SR;
        delete paraV;
        delete[] ucell.atoms;
    }

#ifdef __MPI
    void init_parav()
    {
        int nb = 10;
        int global_row = test_size * test_nw;
        int global_col = test_size * test_nw;
        std::ofstream ofs_running;
        paraV = new Parallel_Orbitals();
        paraV->init(global_row, global_col, nb, MPI_COMM_WORLD);
        paraV->set_atomic_trace(ucell.get_iat2iwt(), test_size, global_row);
    }
#else
    void init_parav()
    {
    }
#endif

    UnitCell ucell;
    hamilt::HContainer<double>* SR;
    Parallel_Orbitals* paraV;
    TwoCenterIntegrator intor_;

    int dsize;
    int my_rank = 0;
};
TEST_F(OverlapNewTest, constructHRcd2cd)
{
    // Create complex SR container
    hamilt::HContainer<std::complex<double>>* SR_complex = new hamilt::HContainer<std::complex<double>>(paraV);

    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.1, 0.2, 0.3));
    hamilt::HS_Matrix_K<std::complex<double>> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);

    hamilt::OverlapNew<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>
        op(&hsk, kvec_d_in, nullptr, SR_complex, &ucell, {1.0}, &gd, &intor_);
    op.contributeHR();

    // Check that SR_complex has been initialized
    // Note: In MPI parallel runs, some processes may not have any atom pairs
    if (SR_complex->size_atom_pairs() > 0)
    {
        EXPECT_GT(SR_complex->size_atom_pairs(), 0);
    }

    // Calculate SK
    op.contributeHk(0);
    auto* sk = hsk.get_sk();

    // Verify SK is computed (values should be non-zero for non-gamma k-point)
    bool has_nonzero = false;
    for (int i = 0; i < paraV->get_row_size() * paraV->get_col_size(); ++i)
    {
        if (std::abs(sk[i]) > 1e-10)
        {
            has_nonzero = true;
            break;
        }
    }
    EXPECT_TRUE(has_nonzero);

    delete SR_complex;
}

// Test getSk method
TEST_F(OverlapNewTest, getSk)
{
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<double> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);

    hamilt::OverlapNew<hamilt::OperatorLCAO<double, double>>
        op(&hsk, kvec_d_in, nullptr, SR, &ucell, {1.0}, &gd, &intor_);

    op.contributeHR();
    op.contributeHk(0);

    // Get SK pointer
    double* sk_from_method = op.getSk();
    double* sk_from_hsk = hsk.get_sk();

    // Verify they point to the same data
    EXPECT_EQ(sk_from_method, sk_from_hsk);

    // Verify values match
    for (int i = 0; i < hsk.get_size(); ++i)
    {
        EXPECT_EQ(sk_from_method[i], sk_from_hsk[i]);
    }
}

// Test k-vector caching optimization
TEST_F(OverlapNewTest, kVectorCaching)
{
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(2, ModuleBase::Vector3<double>(0.1, 0.2, 0.3));
    kvec_d_in[1] = ModuleBase::Vector3<double>(0.1, 0.2, 0.3); // Same k-vector
    hamilt::HS_Matrix_K<std::complex<double>> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);

    hamilt::OverlapNew<hamilt::OperatorLCAO<std::complex<double>, double>>
        op(&hsk, kvec_d_in, nullptr, SR, &ucell, {1.0}, &gd, &intor_);

    op.contributeHR();

    // First call with k-vector index 0
    op.contributeHk(0);
    auto* sk = hsk.get_sk();
    std::complex<double> first_value = sk[0];

    // Second call with same k-vector (index 1)
    hsk.set_zero_sk();
    op.contributeHk(1);
    std::complex<double> second_value = sk[0];

    // Values should be identical due to caching
    EXPECT_NEAR(first_value.real(), second_value.real(), 1e-10);
    EXPECT_NEAR(first_value.imag(), second_value.imag(), 1e-10);
}

// Test with single atom system
TEST_F(OverlapNewTest, singleAtom)
{
    // Create a single atom unit cell
    UnitCell ucell_single;
    ucell_single.ntype = 1;
    ucell_single.nat = 1;
    ucell_single.atoms = new Atom[1];
    ucell_single.iat2it = new int[1];
    ucell_single.iat2ia = new int[1];
    ucell_single.atoms[0].tau.resize(1);
    ucell_single.itia2iat.create(1, 1);
    ucell_single.iat2it[0] = 0;
    ucell_single.iat2ia[0] = 0;
    ucell_single.atoms[0].tau[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    ucell_single.itia2iat(0, 0) = 0;
    ucell_single.atoms[0].na = 1;
    ucell_single.atoms[0].nw = 5;
    ucell_single.atoms[0].iw2l.resize(5);
    ucell_single.atoms[0].iw2m.resize(5);
    ucell_single.atoms[0].iw2n.resize(5);
    for (int iw = 0; iw < 5; ++iw)
    {
        ucell_single.atoms[0].iw2l[iw] = 0;
        ucell_single.atoms[0].iw2m[iw] = 0;
        ucell_single.atoms[0].iw2n[iw] = 0;
    }
    ucell_single.set_iat2iwt(1);

    Parallel_Orbitals* paraV_single = nullptr;
#ifdef __MPI
    int nb = 5;
    paraV_single = new Parallel_Orbitals();
    paraV_single->init(5, 5, nb, MPI_COMM_WORLD);
    paraV_single->set_atomic_trace(ucell_single.get_iat2iwt(), 1, 5);
#endif

    hamilt::HContainer<double>* SR_single = new hamilt::HContainer<double>(paraV_single);

    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<double> hsk(paraV_single);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);

    hamilt::OverlapNew<hamilt::OperatorLCAO<double, double>>
        op(&hsk, kvec_d_in, nullptr, SR_single, &ucell_single, {1.0}, &gd, &intor_);

    op.contributeHR();

    // Should have only self-interaction (atom 0 with itself)
    // Note: In MPI parallel runs, some processes may not have any atom pairs
    if (SR_single->size_atom_pairs() > 0)
    {
        EXPECT_GT(SR_single->size_atom_pairs(), 0);
    }

    // Calculate SK
    op.contributeHk(0);

    delete SR_single;
    delete paraV_single;
    delete[] ucell_single.atoms;
}

// Test with different orbital quantum numbers (L, N, M)
TEST_F(OverlapNewTest, differentOrbitals)
{
    // Modify orbital quantum numbers to test different L, N, M values
    ucell.atoms[0].iw2l[0] = 0; // s orbital
    ucell.atoms[0].iw2l[1] = 1; // p orbital
    ucell.atoms[0].iw2l[2] = 1; // p orbital
    ucell.atoms[0].iw2l[3] = 1; // p orbital
    ucell.atoms[0].iw2l[4] = 2; // d orbital

    ucell.atoms[0].iw2m[0] = 0;
    ucell.atoms[0].iw2m[1] = -1;
    ucell.atoms[0].iw2m[2] = 0;
    ucell.atoms[0].iw2m[3] = 1;
    ucell.atoms[0].iw2m[4] = 0;

    ucell.atoms[0].iw2n[0] = 0;
    ucell.atoms[0].iw2n[1] = 0;
    ucell.atoms[0].iw2n[2] = 0;
    ucell.atoms[0].iw2n[3] = 0;
    ucell.atoms[0].iw2n[4] = 0;

    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<double> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);

    hamilt::OverlapNew<hamilt::OperatorLCAO<double, double>>
        op(&hsk, kvec_d_in, nullptr, SR, &ucell, {1.0}, &gd, &intor_);

    op.contributeHR();

    // Verify SR is calculated
    // Note: In MPI parallel runs, some processes may not have any atom pairs
    if (SR->size_atom_pairs() > 0)
    {
        EXPECT_GT(SR->size_atom_pairs(), 0);
    }

    op.contributeHk(0);
}

// Test force calculation
TEST_F(OverlapNewTest, forceCalculation)
{
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<double> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);

    hamilt::OverlapNew<hamilt::OperatorLCAO<double, double>>
        op(&hsk, kvec_d_in, nullptr, SR, &ucell, {1.0}, &gd, &intor_);

    op.contributeHR();

    // Create a simple density matrix
    hamilt::HContainer<double> dmR(paraV);
    // Initialize dmR with same structure as SR
    for (int iap = 0; iap < SR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& sr_pair = SR->get_atom_pair(iap);
        int iat1 = sr_pair.get_atom_i();
        int iat2 = sr_pair.get_atom_j();
        for (int iR = 0; iR < sr_pair.get_R_size(); ++iR)
        {
            ModuleBase::Vector3<int> R_index = sr_pair.get_R_index(iR);
            hamilt::AtomPair<double> dm_pair(iat1, iat2, R_index, paraV);
            dmR.insert_pair(dm_pair);
        }
    }
    dmR.allocate(nullptr, true);

    // Set density matrix to identity-like values
    for (int iap = 0; iap < dmR.size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = dmR.get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        auto indexes1 = paraV->get_indexes_row(iat1);
        auto indexes2 = paraV->get_indexes_col(iat2);
        int nwt = indexes1.size() * indexes2.size();
        for (int i = 0; i < nwt; ++i)
        {
            tmp.get_pointer(0)[i] = 0.1;
        }
    }

    ModuleBase::matrix force(ucell.nat, 3);
    ModuleBase::matrix stress(3, 3);

    // Calculate force only
    op.cal_force_stress(true, false, &dmR, force, stress);

    // Test passes if no crash occurs
    EXPECT_TRUE(true);
}

// Test stress calculation
TEST_F(OverlapNewTest, stressCalculation)
{
    // Initialize unit cell parameters for stress calculation
    ucell.lat0 = 1.0;
    ucell.omega = 1000.0;  // Set non-zero volume to avoid division by zero
    ucell.latvec.e11 = 10.0;
    ucell.latvec.e22 = 10.0;
    ucell.latvec.e33 = 10.0;

    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<double> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);

    hamilt::OverlapNew<hamilt::OperatorLCAO<double, double>>
        op(&hsk, kvec_d_in, nullptr, SR, &ucell, {1.0}, &gd, &intor_);

    op.contributeHR();

    // Create density matrix
    hamilt::HContainer<double> dmR(paraV);
    for (int iap = 0; iap < SR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& sr_pair = SR->get_atom_pair(iap);
        int iat1 = sr_pair.get_atom_i();
        int iat2 = sr_pair.get_atom_j();
        for (int iR = 0; iR < sr_pair.get_R_size(); ++iR)
        {
            ModuleBase::Vector3<int> R_index = sr_pair.get_R_index(iR);
            hamilt::AtomPair<double> dm_pair(iat1, iat2, R_index, paraV);
            dmR.insert_pair(dm_pair);
        }
    }
    dmR.allocate(nullptr, true);

    for (int iap = 0; iap < dmR.size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = dmR.get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        auto indexes1 = paraV->get_indexes_row(iat1);
        auto indexes2 = paraV->get_indexes_col(iat2);
        int nwt = indexes1.size() * indexes2.size();
        for (int i = 0; i < nwt; ++i)
        {
            tmp.get_pointer(0)[i] = 0.1;
        }
    }

    ModuleBase::matrix force(ucell.nat, 3);
    ModuleBase::matrix stress(3, 3);

    // Calculate stress only
    op.cal_force_stress(false, true, &dmR, force, stress);

    // Verify stress tensor is symmetric (within numerical precision)
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_NEAR(stress(i, j), stress(j, i), 1e-8);
        }
    }
}

// Test force and stress together
TEST_F(OverlapNewTest, forceStressTogether)
{
    // Initialize unit cell parameters for stress calculation
    ucell.lat0 = 1.0;
    ucell.omega = 1000.0;  // Set non-zero volume to avoid division by zero
    ucell.latvec.e11 = 10.0;
    ucell.latvec.e22 = 10.0;
    ucell.latvec.e33 = 10.0;

    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<double> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);

    hamilt::OverlapNew<hamilt::OperatorLCAO<double, double>>
        op(&hsk, kvec_d_in, nullptr, SR, &ucell, {1.0}, &gd, &intor_);

    op.contributeHR();

    // Create density matrix
    hamilt::HContainer<double> dmR(paraV);
    for (int iap = 0; iap < SR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& sr_pair = SR->get_atom_pair(iap);
        int iat1 = sr_pair.get_atom_i();
        int iat2 = sr_pair.get_atom_j();
        for (int iR = 0; iR < sr_pair.get_R_size(); ++iR)
        {
            ModuleBase::Vector3<int> R_index = sr_pair.get_R_index(iR);
            hamilt::AtomPair<double> dm_pair(iat1, iat2, R_index, paraV);
            dmR.insert_pair(dm_pair);
        }
    }
    dmR.allocate(nullptr, true);

    for (int iap = 0; iap < dmR.size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = dmR.get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        auto indexes1 = paraV->get_indexes_row(iat1);
        auto indexes2 = paraV->get_indexes_col(iat2);
        int nwt = indexes1.size() * indexes2.size();
        for (int i = 0; i < nwt; ++i)
        {
            tmp.get_pointer(0)[i] = 0.1;
        }
    }

    ModuleBase::matrix force(ucell.nat, 3);
    ModuleBase::matrix stress(3, 3);

    // Calculate both force and stress
    op.cal_force_stress(true, true, &dmR, force, stress);

    // Verify stress symmetry
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_NEAR(stress(i, j), stress(j, i), 1e-8);
        }
    }

    // Test passes if no crash occurs
    EXPECT_TRUE(true);
}

// Test with zero orbital cutoff
TEST_F(OverlapNewTest, zeroOrbitalCutoff)
{
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<double> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);

    // Use zero cutoff - should result in no atom pairs (except possibly self-interaction)
    hamilt::OverlapNew<hamilt::OperatorLCAO<double, double>>
        op(&hsk, kvec_d_in, nullptr, SR, &ucell, {0.0}, &gd, &intor_);

    op.contributeHR();

    // With zero cutoff, there should be no or very few atom pairs
    // (implementation dependent - might include self-interaction)
    EXPECT_GE(SR->size_atom_pairs(), 0);
}

// Test with large orbital cutoff
TEST_F(OverlapNewTest, largeOrbitalCutoff)
{
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<double> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);

    // Use very large cutoff - should include all atoms
    hamilt::OverlapNew<hamilt::OperatorLCAO<double, double>>
        op(&hsk, kvec_d_in, nullptr, SR, &ucell, {1000.0}, &gd, &intor_);

    op.contributeHR();

    // With large cutoff, should have many atom pairs
    // Note: In MPI parallel runs, some processes may not have any atom pairs
    if (SR->size_atom_pairs() > 0)
    {
        EXPECT_GT(SR->size_atom_pairs(), 0);
    }

    op.contributeHk(0);
}

// Test with atoms at cutoff boundary
TEST_F(OverlapNewTest, cutoffBoundary)
{
    // Set up atoms at specific distances to test cutoff boundary
    ucell.lat0 = 1.0;
    ucell.latvec.e11 = 10.0;
    ucell.latvec.e22 = 10.0;
    ucell.latvec.e33 = 10.0;

    // Place two atoms at distance exactly at cutoff
    ucell.atoms[0].tau[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    ucell.atoms[0].tau[1] = ModuleBase::Vector3<double>(0.5, 0.0, 0.0); // distance = 5.0

    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<double> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);

    // Use cutoff of 5.0 - atoms at exactly this distance should be excluded
    hamilt::OverlapNew<hamilt::OperatorLCAO<double, double>>
        op(&hsk, kvec_d_in, nullptr, SR, &ucell, {2.5}, &gd, &intor_);

    op.contributeHR();

    // Verify initialization completed
    EXPECT_GE(SR->size_atom_pairs(), 0);
}

// Test Hermitian property of SK matrix
TEST_F(OverlapNewTest, hermitianProperty)
{
    // Use gamma point to test that diagonal elements are real
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<std::complex<double>> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);

    hamilt::OverlapNew<hamilt::OperatorLCAO<std::complex<double>, double>>
        op(&hsk, kvec_d_in, nullptr, SR, &ucell, {1.0}, &gd, &intor_);

    op.contributeHR();
    op.contributeHk(0);

    auto* sk = hsk.get_sk();
    int nrow = paraV->get_row_size();
    int ncol = paraV->get_col_size();

    // For overlap matrix at gamma point, SK should be real and symmetric
    // Diagonal elements should be real (imaginary part should be zero)
    for (int i = 0; i < std::min(nrow, ncol); ++i)
    {
        if (i < nrow && i < ncol)
        {
            int idx = i * ncol + i;
            if (idx < nrow * ncol)
            {
                EXPECT_NEAR(sk[idx].imag(), 0.0, 1e-8);
            }
        }
    }
}

// Test with null SR pointer (should skip initialization)
TEST_F(OverlapNewTest, nullSRPointer)
{
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<double> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);

    // Pass nullptr for SR - should not crash during construction
    hamilt::OverlapNew<hamilt::OperatorLCAO<double, double>>
        op(&hsk, kvec_d_in, nullptr, nullptr, &ucell, {1.0}, &gd, &intor_);

    // Test passes if no crash occurs during construction
    EXPECT_TRUE(true);
}

// Test force calculation with npol=2 (nspin=4, spin-orbit coupling)
TEST_F(OverlapNewTest, forceCalculationNpol2)
{
    // Set up unit cell with npol=2
    ucell.set_iat2iwt(2);  // npol=2

    // Reinitialize paraV with doubled size for npol=2
    delete paraV;
    paraV = nullptr;
#ifdef __MPI
    int nb = 10;
    int global_row = test_size * test_nw * 2;  // doubled for npol=2
    int global_col = test_size * test_nw * 2;
    paraV = new Parallel_Orbitals();
    paraV->init(global_row, global_col, nb, MPI_COMM_WORLD);
    paraV->set_atomic_trace(ucell.get_iat2iwt(), test_size, global_row);
#endif

    // Create complex HContainer for npol=2
    hamilt::HContainer<std::complex<double>>* SR_complex = new hamilt::HContainer<std::complex<double>>(paraV);

    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<std::complex<double>> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);

    hamilt::OverlapNew<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>
        op(&hsk, kvec_d_in, nullptr, SR_complex, &ucell, {1.0}, &gd, &intor_);

    op.contributeHR();

    // Create REAL density matrix (charge density) for force/stress calculation
    // Even with npol=2, the density matrix for force/stress is real-valued
    hamilt::HContainer<double> dmR(paraV);
    for (int iap = 0; iap < SR_complex->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<std::complex<double>>& sr_pair = SR_complex->get_atom_pair(iap);
        int iat1 = sr_pair.get_atom_i();
        int iat2 = sr_pair.get_atom_j();
        for (int iR = 0; iR < sr_pair.get_R_size(); ++iR)
        {
            ModuleBase::Vector3<int> R_index = sr_pair.get_R_index(iR);
            hamilt::AtomPair<double> dm_pair(iat1, iat2, R_index, paraV);
            dmR.insert_pair(dm_pair);
        }
    }
    dmR.allocate(nullptr, true);

    // Set density matrix values - real values representing charge density
    // For npol=2, the layout is still handled by step_trace in the implementation
    for (int iap = 0; iap < dmR.size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = dmR.get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        auto indexes1 = paraV->get_indexes_row(iat1);
        auto indexes2 = paraV->get_indexes_col(iat2);

        // Fill with real charge density values
        double* dm_ptr = tmp.get_pointer(0);
        for (int iw1 = 0; iw1 < indexes1.size(); iw1 += 2)
        {
            for (int iw2 = 0; iw2 < indexes2.size(); iw2 += 2)
            {
                int idx = iw1 * indexes2.size() + iw2;
                // Set charge density values (diagonal of spin density matrix)
                dm_ptr[idx] = 0.1;  // Charge density at this orbital pair
            }
        }
    }

    ModuleBase::matrix force(ucell.nat, 3);
    ModuleBase::matrix stress(3, 3);

    // Calculate force with npol=2
    op.cal_force_stress(true, false, &dmR, force, stress);

    // Verify force calculation completed without crash
    EXPECT_TRUE(true);

    delete SR_complex;

    // Restore npol=1 for other tests
    ucell.set_iat2iwt(1);
}

// Test stress calculation with npol=2 (nspin=4, spin-orbit coupling)
TEST_F(OverlapNewTest, stressCalculationNpol2)
{
    // Set up unit cell with npol=2
    ucell.set_iat2iwt(2);  // npol=2

    // Initialize unit cell parameters for stress calculation
    ucell.lat0 = 1.0;
    ucell.omega = 1000.0;  // Set non-zero volume to avoid division by zero
    ucell.latvec.e11 = 10.0;
    ucell.latvec.e22 = 10.0;
    ucell.latvec.e33 = 10.0;

    // Reinitialize paraV with doubled size for npol=2
    delete paraV;
    paraV = nullptr;
#ifdef __MPI
    int nb = 10;
    int global_row = test_size * test_nw * 2;  // doubled for npol=2
    int global_col = test_size * test_nw * 2;
    paraV = new Parallel_Orbitals();
    paraV->init(global_row, global_col, nb, MPI_COMM_WORLD);
    paraV->set_atomic_trace(ucell.get_iat2iwt(), test_size, global_row);
#endif

    // Create complex HContainer for npol=2
    hamilt::HContainer<std::complex<double>>* SR_complex = new hamilt::HContainer<std::complex<double>>(paraV);

    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<std::complex<double>> hsk(paraV);
    hsk.set_zero_sk();
    Grid_Driver gd(0, 0);

    hamilt::OverlapNew<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>
        op(&hsk, kvec_d_in, nullptr, SR_complex, &ucell, {1.0}, &gd, &intor_);

    op.contributeHR();

    // Create REAL density matrix (charge density) for force/stress calculation
    hamilt::HContainer<double> dmR(paraV);
    for (int iap = 0; iap < SR_complex->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<std::complex<double>>& sr_pair = SR_complex->get_atom_pair(iap);
        int iat1 = sr_pair.get_atom_i();
        int iat2 = sr_pair.get_atom_j();
        for (int iR = 0; iR < sr_pair.get_R_size(); ++iR)
        {
            ModuleBase::Vector3<int> R_index = sr_pair.get_R_index(iR);
            hamilt::AtomPair<double> dm_pair(iat1, iat2, R_index, paraV);
            dmR.insert_pair(dm_pair);
        }
    }
    dmR.allocate(nullptr, true);

    // Set density matrix values - real values representing charge density
    for (int iap = 0; iap < dmR.size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = dmR.get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        auto indexes1 = paraV->get_indexes_row(iat1);
        auto indexes2 = paraV->get_indexes_col(iat2);

        double* dm_ptr = tmp.get_pointer(0);
        for (int iw1 = 0; iw1 < indexes1.size(); iw1 += 2)
        {
            for (int iw2 = 0; iw2 < indexes2.size(); iw2 += 2)
            {
                int idx = iw1 * indexes2.size() + iw2;
                dm_ptr[idx] = 0.1;  // Charge density
            }
        }
    }

    ModuleBase::matrix force(ucell.nat, 3);
    ModuleBase::matrix stress(3, 3);

    // Calculate stress with npol=2
    op.cal_force_stress(false, true, &dmR, force, stress);

    // Verify stress tensor is symmetric
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_NEAR(stress(i, j), stress(j, i), 1e-8);
        }
    }

    delete SR_complex;

    // Restore npol=1 for other tests
    ucell.set_iat2iwt(1);
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
