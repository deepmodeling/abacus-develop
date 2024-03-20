#include "gtest/gtest.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "time.h"

// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 100;

/**
 * Unit test of HContainer for time consuming
 * HContainer is a container of hamilt::AtomPair, in this test, we test the following functions:
 * 1. insert_pair
 * 2. find_pair
 *
 */

// Unit test of HContainer with gtest framework
class HContainerTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        // set up a unitcell, with one element and three atoms, each atom has 2 orbitals
        ucell.ntype = 1;
        ucell.nat = test_size;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.iat2it = new int[ucell.nat];
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2it[iat] = 0;
        }
        ucell.atoms[0].nw = 10;

        // set up a HContainer with ucell
        HR = new hamilt::HContainer<double>(ucell);
    }

    void TearDown() override
    {
        delete HR;
        delete[] ucell.atoms;
        delete[] ucell.iat2it;
    }

    UnitCell ucell;
    hamilt::HContainer<double>* HR;
};

// using TEST to test HContainer::insert_pair

TEST(single_test_hcontainer, insert_pair)
{
    // print current used memory of system
    auto memory_start = ModuleBase::GlobalFunc::MemAvailable() / 1024.0;
    // get random number between (0, 100000000) by rand()
    srand((unsigned)time(NULL));
    hamilt::HContainer<double> HR_test(test_size);
    clock_t start, end;
    start = clock();
    for (int i = 0; i < test_size * 100; i++)
    {
        int iat1 = rand() % test_size;
        int iat2 = rand() % test_size;
        hamilt::AtomPair<double> tmp(iat1, iat2);
        tmp.set_size(10, 10);
        auto p = tmp.get_HR_values(0, 0, 0).get_pointer();
        p[0] = i * 1.0;
        HR_test.insert_pair(tmp);
    }
    end = clock();
    std::cout << "total atom pairs: " << HR_test.size_atom_pairs() << std::endl;
    int count = 0;
    for (int i = 0; i < HR_test.size_atom_pairs(); i++)
    {
        auto tmp = HR_test.get_atom_pair(i);
        count += tmp.get_R_size();
    }
    std::cout << "total I-J-R: " << count << std::endl;
    std::cout << "pure data size of I-J-R: " << count * 10 * 10 * sizeof(double) / 1024.0 / 1024.0 << " MB"
              << std::endl;
    auto memory_end = ModuleBase::GlobalFunc::MemAvailable() / 1024.0;
    std::cout << "memory: " << double(HR_test.get_memory_size()) / 1024.0 / 1024.0 << " MB" << std::endl;
    std::cout << "time: " << (double)(end - start) / CLOCKS_PER_SEC << std::endl;
}

TEST(single_test_IJR, insert_pair)
{
    // get random number between (0, 100000000) by rand()
    srand((unsigned)time(NULL));
    hamilt::HContainer<double> HR_test(test_size);
    clock_t start, end;
    start = clock();
    for (int i = 0; i < test_size; i++)
    {
        int iat1 = rand() % test_size;
        int iat2 = rand() % test_size;
        hamilt::AtomPair<double> tmp(iat1, iat2);
        tmp.set_size(10, 10);
        for (int j = 0; j < 100; j++)
        {
            auto p = tmp.get_HR_values(0, 0, j).get_pointer();
            p[0] = i * 1.0;
        }
        HR_test.insert_pair(tmp);
    }
    end = clock();
    std::cout << "total atom pairs: " << HR_test.size_atom_pairs() << std::endl;
    int count = 0;
    for (int i = 0; i < HR_test.size_atom_pairs(); i++)
    {
        auto tmp = HR_test.get_atom_pair(i);
        count += tmp.get_R_size();
    }
    std::cout << "total I-J-R: " << count << std::endl;
    std::cout << "pure data size of I-J-R: " << count * 10 * 10 * sizeof(double) / 1024.0 / 1024.0 << " MB"
              << std::endl;
    std::cout << "memory: " << double(HR_test.get_memory_size()) / 1024.0 / 1024.0 << " MB" << std::endl;
    std::cout << "time: " << (double)(end - start) / CLOCKS_PER_SEC << std::endl;
}

// using TEST_F to test HContainer::find_pair
TEST_F(HContainerTest, find_pair)
{
    // find_pair 1000000 times with random iat1 and iat2
    srand((unsigned)time(NULL));
    clock_t start, end;
    start = clock();
    for (int i = 0; i < test_size * 100; i++)
    {
        int iat1 = rand() % test_size;
        int iat2 = rand() % test_size;
        auto tmp = HR->find_pair(iat1, iat2);
        if (tmp == nullptr)
        {
            std::cout << "not found" << std::endl;
        }
    }
    end = clock();
    std::cout << "time: " << (double)(end - start) / CLOCKS_PER_SEC << std::endl;
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}