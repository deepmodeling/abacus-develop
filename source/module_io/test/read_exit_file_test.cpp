#include "module_io/read_exit_file.h"
#include "module_io/read_input.h"
#include "mpi.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <fstream>

/************************************************
 *  unit test of read_exit_file.cpp
 ***********************************************/

class ReadStopFileTest : public testing::Test
{
  protected:
    virtual void SetUp()
    {
    }
    virtual void TearDown()
    {
    }
};

TEST_F(ReadStopFileTest, read_exit_file)
{
    std::string filename = "EXIT";
    std::string output = "running.txt";
    std::ofstream ofs_running(output.c_str(), std::ios::out);

    // case 1: no such file
    int stop = ModuleIO::read_exit_file(GlobalV::MY_RANK, filename, ofs_running);
    EXPECT_EQ(stop, 0);
    MPI_Barrier(MPI_COMM_WORLD);

    // case 2: no keywords in the file
    std::ofstream ofs;
    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(filename.c_str(), std::ios::out);
        ofs << "no keywords\n\ntest" << std::endl;
        ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    stop = ModuleIO::read_exit_file(GlobalV::MY_RANK, filename, ofs_running);
    EXPECT_EQ(stop, 0);
    MPI_Barrier(MPI_COMM_WORLD);

    // case 3:  stop_ion = false    stop_elec = false
    if (GlobalV::MY_RANK == 0)
    {
        std::remove(filename.c_str());
        ofs.open(filename.c_str(), std::ios::out);
        ofs << "stop_ion    false\nstop_elec    0" << std::endl;
        ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    stop = ModuleIO::read_exit_file(GlobalV::MY_RANK, filename, ofs_running);
    EXPECT_EQ(stop, 0);
    MPI_Barrier(MPI_COMM_WORLD);

    // case 4:  stop_ion = true    stop_elec = false
    if (GlobalV::MY_RANK == 0)
    {
        std::remove(filename.c_str());
        ofs.open(filename.c_str(), std::ios::out);
        ofs << "stop_ion    true\nstop_elec    f" << std::endl;
        ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    stop = ModuleIO::read_exit_file(GlobalV::MY_RANK, filename, ofs_running);
    EXPECT_EQ(stop, 1);
    MPI_Barrier(MPI_COMM_WORLD);

    // case 5:  stop_ion = false    stop_elec = true
    if (GlobalV::MY_RANK == 0)
    {
        std::remove(filename.c_str());
        ofs.open(filename.c_str(), std::ios::out);
        ofs << "stop_ion    F\nstop_elec    1" << std::endl;
        ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    stop = ModuleIO::read_exit_file(GlobalV::MY_RANK, filename, ofs_running);
    EXPECT_EQ(stop, 2);
    MPI_Barrier(MPI_COMM_WORLD);

    // case 6:  stop_ion = true    stop_elec = true
    if (GlobalV::MY_RANK == 0)
    {
        std::remove(filename.c_str());
        ofs.open(filename.c_str(), std::ios::out);
        ofs << "stop_ion    T\nstop_elec    t" << std::endl;
        ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    stop = ModuleIO::read_exit_file(GlobalV::MY_RANK, filename, ofs_running);
    EXPECT_EQ(stop, 2);
    MPI_Barrier(MPI_COMM_WORLD);

    ofs_running.close();
    std::remove(output.c_str());
    std::remove(filename.c_str());
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