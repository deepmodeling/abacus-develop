#include "source_base/global_variable.h"
#include "source_base/parallel_global.h"
#include "source_cell/klist.h"
#include "source_estate/elecstate_tools.h"
#include "source_estate/occupy.h"
#include "source_io/module_parameter/parameter.h"

#include "gtest/gtest.h"
#include <mpi.h>
#include <vector>

class TestParameters
{
  public:
    static void configure_weights()
    {
        PARAM.input.nbands = 5;
        PARAM.input.nelec = 6.0;
        PARAM.input.nspin = 1;
        PARAM.sys.two_fermi = false;
    }

    static void configure_band_parallelism()
    {
        PARAM.input.bndpar = 2;
    }
};

namespace
{
constexpr int GLOBAL_NBANDS = 5;

void initialize_klist(K_Vectors* klist)
{
    klist->set_nks(1);
    klist->set_nkstot(1);
    klist->ik2iktot.push_back(0);
    klist->wk.push_back(1.0);
    klist->isk.push_back(0);
}

class BandParallelWeightsTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        TestParameters::configure_weights();
        Occupy::use_gaussian_broadening = false;
        Occupy::fixed_occupations = false;
    }
};

TEST_F(BandParallelWeightsTest, FixedWeightsReplicated)
{
    K_Vectors klist;
    initialize_klist(&klist);
    ModuleBase::matrix weights(1, GLOBAL_NBANDS);
    const std::vector<double> occupations{0.5, 0.4, 0.3, 0.2, 0.1};
    bool skip_weights = false;

    elecstate::fixed_weights(occupations, GLOBAL_NBANDS, 1.5, &klist, weights, skip_weights);

    for (int band = 0; band < GLOBAL_NBANDS; ++band)
    {
        EXPECT_DOUBLE_EQ(weights(0, band), occupations[band]);
    }
    EXPECT_TRUE(skip_weights);
}

TEST_F(BandParallelWeightsTest, FixedWeightsDistributed)
{
    K_Vectors klist;
    initialize_klist(&klist);
    const int local_nbands = GlobalV::MY_BNDGROUP == 0 ? 3 : 2;
    const int band_offset = GlobalV::MY_BNDGROUP == 0 ? 0 : 3;
    ModuleBase::matrix weights(1, local_nbands);
    const std::vector<double> occupations{0.5, 0.4, 0.3, 0.2, 0.1};
    bool skip_weights = false;

    elecstate::fixed_weights(occupations, GLOBAL_NBANDS, 1.5, &klist, weights, skip_weights);

    for (int band = 0; band < local_nbands; ++band)
    {
        EXPECT_DOUBLE_EQ(weights(0, band), occupations[band_offset + band]);
    }
    EXPECT_TRUE(skip_weights);
}

TEST_F(BandParallelWeightsTest, IntegerWeightsReplicated)
{
    K_Vectors klist;
    initialize_klist(&klist);
    ModuleBase::matrix eigenvalues(1, GLOBAL_NBANDS);
    ModuleBase::matrix weights(1, GLOBAL_NBANDS);
    for (int band = 0; band < GLOBAL_NBANDS; ++band)
    {
        eigenvalues(0, band) = static_cast<double>(band);
    }
    elecstate::Efermi eferm;
    elecstate::fenergy energy;
    std::vector<double> nelec_spin;

    elecstate::calculate_weights(eigenvalues, weights, &klist, eferm, energy, nelec_spin, GLOBAL_NBANDS, false);

    for (int band = 0; band < GLOBAL_NBANDS; ++band)
    {
        EXPECT_DOUBLE_EQ(weights(0, band), band < 3 ? 1.0 : 0.0);
    }
    EXPECT_DOUBLE_EQ(eferm.ef, 2.0);
}

TEST_F(BandParallelWeightsTest, IntegerWeightsDistributed)
{
    K_Vectors klist;
    initialize_klist(&klist);
    const int local_nbands = GlobalV::MY_BNDGROUP == 0 ? 3 : 2;
    const int band_offset = GlobalV::MY_BNDGROUP == 0 ? 0 : 3;
    ModuleBase::matrix eigenvalues(1, local_nbands);
    ModuleBase::matrix weights(1, local_nbands);
    for (int band = 0; band < local_nbands; ++band)
    {
        eigenvalues(0, band) = static_cast<double>(band_offset + band);
    }
    elecstate::Efermi eferm;
    elecstate::fenergy energy;
    std::vector<double> nelec_spin;

    elecstate::calculate_weights(eigenvalues, weights, &klist, eferm, energy, nelec_spin, GLOBAL_NBANDS, false);

    for (int band = 0; band < local_nbands; ++band)
    {
        EXPECT_DOUBLE_EQ(weights(0, band), band_offset + band < 3 ? 1.0 : 0.0);
    }
    EXPECT_DOUBLE_EQ(eferm.ef, 2.0);
}
} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);

    TestParameters::configure_band_parallelism();
    GlobalV::KPAR = 1;
    Parallel_Global::init_pools(GlobalV::NPROC,
                                GlobalV::MY_RANK,
                                PARAM.inp.bndpar,
                                GlobalV::KPAR,
                                GlobalV::NPROC_IN_BNDGROUP,
                                GlobalV::RANK_IN_BPGROUP,
                                GlobalV::MY_BNDGROUP,
                                GlobalV::NPROC_IN_POOL,
                                GlobalV::RANK_IN_POOL,
                                GlobalV::MY_POOL);

    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
