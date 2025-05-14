#include "module_base/global_variable.h"

#define private public
#include "module_parameter/parameter.h"
#undef private
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <streambuf>
#ifdef __MPI
#include "module_base/parallel_global.h"
#include "module_cell/parallel_kpoints.h"
#include "mpi.h"
#endif
#include "../write_istate_info.h"
#include "for_testing_klist.h"

/************************************************
 *  unit test of write_istate_info
 ***********************************************/

/**
 * - Tested Functions:
 *   - write_istate_info()
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
    // preconditions
    GlobalV::KPAR = 1;
    PARAM.input.nbands = 4;
    PARAM.sys.nbands_l = 4;
    PARAM.input.nspin = 1;
    PARAM.sys.global_out_dir = "./";
    // mpi setting
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
    kv->set_nkstot(100);
    int nkstot = kv->get_nkstot();
    kv->para_k.kinfo(nkstot, GlobalV::KPAR, GlobalV::MY_POOL, GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL, PARAM.input.nspin);
    // std::cout<<"my_rank "<<GlobalV::MY_RANK<<" pool rank/size: "
    //	<<GlobalV::RANK_IN_POOL<<"/"<<GlobalV::NPROC_IN_POOL<<std::endl;
    // std::cout<<"MY_POOL "<<GlobalV::MY_POOL<<std::endl;
    kv->set_nks(kv->para_k.nks_pool[GlobalV::MY_POOL]);
    // std::cout<<"nks "<<kv->get_nks()<<std::endl;
    ekb.create(kv->get_nks(), PARAM.input.nbands);
    wg.create(kv->get_nks(), PARAM.input.nbands);
    ekb.fill_out(0.15);
    wg.fill_out(0.0);
    kv->kvec_d.resize(kv->get_nkstot());
    int i = 0;
    for (auto& kd: kv->kvec_d)
    {
        kd.set(0.01 * i, 0.01 * i, 0.01 * i);
        ++i;
    }
    ModuleIO::write_istate_info(ekb, wg, *kv);
    std::ifstream ifs;
    ifs.open("istate.info");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(
        str,
        testing::HasSubstr("BAND               Energy(ev)               Occupation                Kpoint = 100"));
    EXPECT_THAT(str, testing::HasSubstr("(0.99 0.99 0.99)"));
    EXPECT_THAT(str, testing::HasSubstr("4                2.0408547                        0"));
    ifs.close();
    remove("istate.info");
}

TEST_F(IstateInfoTest, OutIstateInfoS2)
{
    // preconditions
    GlobalV::KPAR = 1;
    PARAM.input.nbands = 4;
    PARAM.sys.nbands_l = 4;
    PARAM.input.nspin = 2;
    PARAM.sys.global_out_dir = "./";
    // mpi setting
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
    kv->set_nkstot(100);
    int nkstot = kv->get_nkstot();
    kv->para_k.kinfo(nkstot, GlobalV::KPAR, GlobalV::MY_POOL, GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL, PARAM.input.nspin);
    // std::cout<<"my_rank "<<GlobalV::MY_RANK<<" pool rank/size: "
    //	<<GlobalV::RANK_IN_POOL<<"/"<<GlobalV::NPROC_IN_POOL<<std::endl;
    // std::cout<<"MY_POOL "<<GlobalV::MY_POOL<<std::endl;
    kv->set_nks(kv->para_k.nks_pool[GlobalV::MY_POOL]);
    // std::cout<<"nks "<<kv->get_nks()<<std::endl;
    ekb.create(kv->get_nks(), PARAM.input.nbands);
    wg.create(kv->get_nks(), PARAM.input.nbands);
    ekb.fill_out(0.15);
    wg.fill_out(0.0);
    kv->kvec_d.resize(kv->get_nkstot());
    int i = 0;
    for (auto& kd: kv->kvec_d)
    {
        kd.set(0.01 * i, 0.01 * i, 0.01 * i);
        ++i;
    }
    ModuleIO::write_istate_info(ekb, wg, *kv);
    std::ifstream ifs;
    ifs.open("istate.info");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str,
                testing::HasSubstr("BAND       Spin up Energy(ev)               Occupation     Spin down Energy(ev)    "
                                   "           Occupation                Kpoint = 50"));
    EXPECT_THAT(
        str,
        testing::HasSubstr(
            "4                  2.04085                        0                  2.04085                        0"));
    EXPECT_THAT(str, testing::HasSubstr("(0.49 0.49 0.49)"));
    ifs.close();
    remove("istate.info");
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
