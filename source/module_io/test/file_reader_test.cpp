#include "module_io/file_reader.h"

#include <fstream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

/************************************************
 *  unit test of file_reader.cpp
 ***********************************************/

/**
 * - Tested Functions:
 * - FileReader()
 * - Constructor
 * - isOpen()
 * - Check if file is open
 * - readLine()
 * - Read a line to string stream
 */

class FileReaderTest : public testing::Test
{
  protected:
    std::ofstream ofs;
    std::string filename = "test.txt";
    void SetUp() override
    {
        ofs.open(filename.c_str());
        // write some info to the file
        ofs << "step 0" << std::endl;
        ofs << "matrix dimension: 4" << std::endl;
        ofs << "number of R: 2" << std::endl;
        ofs << "R: 0 0 0 2" << std::endl;
        ofs << "1 2" << std::endl;
        ofs.close();
    }
    void TearDown() override
    {
        remove("test.txt");
    }
    std::string output;
};

TEST_F(FileReaderTest, Constructor)
{
    ModuleIO::FileReader fr(filename);
    // Check if file is open
    EXPECT_TRUE(fr.isOpen());
}

TEST_F(FileReaderTest, ReadLine)
{
    ModuleIO::FileReader fr(filename);
    // Check if file is open
    EXPECT_TRUE(fr.isOpen());
    // read a line to string stream
    fr.readLine();
    EXPECT_EQ(fr.ss.str(), "step 0");
    fr.readLine();
    EXPECT_EQ(fr.ss.str(), "matrix dimension: 4");
    fr.readLine();
    EXPECT_EQ(fr.ss.str(), "number of R: 2");
    fr.readLine();
    EXPECT_EQ(fr.ss.str(), "R: 0 0 0 2");
    fr.readLine();
    EXPECT_EQ(fr.ss.str(), "1 2");
    fr.readLine();
    testing::internal::CaptureStdout();
    EXPECT_EXIT(fr.readLine(), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("End of file"));
}

#include "mpi.h"
int main(int argc, char **argv)
{
        MPI_Init(&argc, &argv);
        testing::InitGoogleTest(&argc, argv);

        int nproc;
        int my_rank;

        MPI_Comm_size(MPI_COMM_WORLD,&nproc);
        MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

        int result = RUN_ALL_TESTS();

        MPI_Finalize();

        return result;
}