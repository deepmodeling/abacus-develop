#include "../spin_constrain.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "prepare_unitcell.h"

#include "mpi.h"

// mock functions
#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
#endif
Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}
// mock functions

/************************************************
 *  unit test of init_sc
 ***********************************************/

/**
 * - Tested functions:
 *  - SpinConstrain::init_sc()
 *    - initialize the SpinConstrain class
 */

K_Vectors::K_Vectors()
{
}
K_Vectors::~K_Vectors()
{
}

template <typename T>
class SpinConstrainTest : public testing::Test
{
  protected:
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["SiO"];
    void SetUp() override
    {
        ucell = utp.SetUcellInfo();
    }
    SpinConstrain<T, psi::DEVICE_CPU>& sc = SpinConstrain<T, psi::DEVICE_CPU>::getScInstance();
};

using MyTypes = ::testing::Types<double, std::complex<double>>;
TYPED_TEST_SUITE(SpinConstrainTest, MyTypes);

TYPED_TEST(SpinConstrainTest, InitSc)
{
    double sc_thr = 1e-6;
    int nsc = 100;
    int nsc_min = 2;
    double alpha_trial = 0.01;
    double sccut = 3.0;
    bool decay_grad_switch = 1;
    K_Vectors kv;
    Parallel_Orbitals paraV;
    paraV.nloc = 2;
    std::string sc_file = "./support/sc_f2.json";
    std::string KS_SOLVER = "genelpa";
    this->sc.init_sc(sc_thr,
                     nsc,
                     nsc_min,
                     alpha_trial,
                     sccut,
                     decay_grad_switch,
                     *(this->ucell),
                     sc_file,
                     2,
                     &paraV,
                     4,
                     kv,
                     KS_SOLVER,
                     nullptr,
                     nullptr,
                     nullptr,
                     nullptr,
                     nullptr);
    EXPECT_EQ(this->sc.get_nat(), 6);
    EXPECT_EQ(this->sc.get_ntype(), 2);
    EXPECT_EQ(this->sc.get_npol(), 2);
    EXPECT_EQ(this->sc.get_nsc(), 100);
    std::map<int, int> atomCounts = this->sc.get_atomCounts();
    EXPECT_EQ(atomCounts[0], 1);
    EXPECT_EQ(atomCounts[1], 5);
    std::map<int, int> orbitalCounts = this->sc.get_orbitalCounts();
    EXPECT_EQ(orbitalCounts[0], 9);
    EXPECT_EQ(orbitalCounts[1], 9);
    std::map<int, std::map<int, int>> LnchiCounts = this->sc.get_lnchiCounts();
    EXPECT_EQ(LnchiCounts[0][0], 1);
    EXPECT_EQ(LnchiCounts[0][1], 1);
    EXPECT_EQ(LnchiCounts[0][2], 1);
    EXPECT_EQ(LnchiCounts[1][0], 1);
    EXPECT_EQ(LnchiCounts[1][1], 1);
    EXPECT_EQ(LnchiCounts[1][2], 1);
}

int main(int argc, char** argv)
{

    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);

    int result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}
