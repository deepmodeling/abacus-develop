#include "source_base/global_variable.h"

#define private public
#include "module_parameter/parameter.h"
#undef private
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <streambuf>
#ifdef __MPI
#include "source_base/parallel_global.h"
#include "source_cell/parallel_kpoints.h"
#include "mpi.h"
#endif
#include "../write_eig_occ.h"
#include "for_testing_klist.h"

/************************************************
 *  unit test of write_eig_occ
 ***********************************************/

/**
 * - Tested Functions:
 *   - write_eig_occ()
 *     - print out electronic eigen energies and
 *     - occupation
 */

class IstateInfoTest : public ::testing::Test
{
  protected:
    K_Vectors* kv = nullptr;
    ModuleBase::matrix ekb;
    ModuleBase::matrix wg;
    void SetUp()
    {
        kv = new K_Vectors;
    }
    void TearDown()
    {
        delete kv;
    }
};

TEST_F(IstateInfoTest, OutIstateInfoS1)
{
    // Global variables 
    GlobalV::KPAR = 1;
    PARAM.input.nbands = 4;
    PARAM.sys.nbands_l = 4;
    PARAM.input.nspin = 1;
    PARAM.sys.global_out_dir = "./";

    // MPI setting
    Parallel_Global::init_pools(GlobalV::NPROC,
                                GlobalV::MY_RANK,
                                PARAM.input.bndpar,
                                GlobalV::KPAR,
                                GlobalV::NPROC_IN_BNDGROUP,
                                GlobalV::RANK_IN_BPGROUP,
                                GlobalV::MY_BNDGROUP,
                                GlobalV::NPROC_IN_POOL,
                                GlobalV::RANK_IN_POOL,
                                GlobalV::MY_POOL);

    const int nkstot_init = 10;
    kv->set_nkstot(nkstot_init);
    int nkstot = kv->get_nkstot();
    kv->para_k.kinfo(nkstot, GlobalV::KPAR, GlobalV::MY_POOL, GlobalV::RANK_IN_POOL, 
    GlobalV::NPROC_IN_POOL, PARAM.input.nspin);
    kv->set_nks(kv->para_k.nks_pool[GlobalV::MY_POOL]);

    // The number of plane waves for each k point
    kv->ngk.resize(nkstot);
    kv->ik2iktot.resize(nkstot);
    for(int i=0; i<nkstot; ++i)
    {
        kv->ngk[i]=299;
        kv->ik2iktot[i]=i;
    }

    // Initialize the number of bands
    ekb.create(kv->get_nks(), PARAM.input.nbands);
    wg.create(kv->get_nks(), PARAM.input.nbands);

    // fill the eigenvalues
    ekb.fill_out(0.15);
    
    // fill the weights
    wg.fill_out(0.0);

    // setup coordinates of k-points
    kv->kvec_c.resize(kv->get_nkstot());
    int i = 0;
    for (auto& kd: kv->kvec_c)
    {
        kd.set(0.01 * i, 0.01 * i, 0.01 * i);
        ++i;
    }
   
    // write eigenvalues and occupations
    ModuleIO::write_eig_file(ekb, wg, *kv);

    // check the output files
    std::ifstream ifs;
    ifs.open("eig.txt");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("Electronic state energy (eV) and occupations"));
    EXPECT_THAT(str, testing::HasSubstr("spin=1 k-point=1/10 Cartesian=0.0000000 0.0000000 0.0000000 (299 plane wave)"));
    EXPECT_THAT(str, testing::HasSubstr("1 2.040854700000000 0.000000000000000"));
    ifs.close();
    remove("eig.txt");
}

#ifdef __MPI
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
#endif
