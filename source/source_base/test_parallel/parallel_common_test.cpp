#ifdef __MPI
#include "source_base/parallel_common.h"

#include "mpi.h"

#include "gtest/gtest.h"
#include <complex>
#include <cstring>
#include <string>

/************************************************
 *  unit test of functions in parallel_common.cpp
 ***********************************************/

/**
 * The tested functions are wrappers of MPI_Bcast
 * in ABACUS, as defined in source_base/parallel_common.h.
 * The source is process 0 in all MPI_Bcast
 * wrappers.
 */

class MPIContext
{
  public:
    MPIContext()
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &_size);
    }

    int GetRank() const
    {
        return _rank;
    }
    int GetSize() const
    {
        return _size;
    }

  private:
    int _rank;
    int _size;
};

class ParaCommon : public testing::Test
{
  protected:
    bool boo = true;
    int is = 0;
    double fs = 0.0;
    std::complex<double> imgs{0.0, 0.0};
    std::string chs = "abacus";
    char cha[7] = "abacus";
    int iv[10] = {0};
    double fv[10] = {0};
    std::complex<double> imgv[10] = {0.0, 0.0};
    std::string chv[10] = {""};
    MPIContext mpiContext;
};

TEST_F(ParaCommon, Bcast)
{
    // reset data in the first process
    int MY_RANK = mpiContext.GetRank();
    if (MY_RANK == 0)
    {
        boo = false;
        is = 1;
        fs = 1.0;
        imgs = std::complex<double>(1.0, -1.0);
        chs = "ABACUS";
        strcpy(cha, chs.c_str());
        for (int i = 0; i < 10; i++)
        {
            double ii = static_cast<double>(i);
            iv[i] = i;
            fv[i] = ii;
            imgv[i] = std::complex<double>(ii, -ii);
            std::stringstream ss;
            ss << chs << i;
            chv[i] = ss.str();
        }
    }
    // call bcast wrappers
    Parallel_Common::bcast_bool(boo);
    Parallel_Common::bcast_int(is);
    Parallel_Common::bcast_double(fs);
    Parallel_Common::bcast_complex_double(imgs);
    Parallel_Common::bcast_string(chs);
    Parallel_Common::bcast_char(cha, 7);
    Parallel_Common::bcast_int(iv, 10);
    Parallel_Common::bcast_double(fv, 10);
    Parallel_Common::bcast_complex_double(imgv, 10);
    Parallel_Common::bcast_string(chv, 10);
    // make comparisons
    EXPECT_FALSE(boo);
    EXPECT_EQ(is, 1);
    EXPECT_EQ(fs, 1.0);
    EXPECT_NEAR(imgs.real(), 1.0, 1E-15);
    EXPECT_NEAR(imgs.imag(), -1.0, 1E-15);
    EXPECT_EQ(chs, "ABACUS");
    EXPECT_STREQ(cha, "ABACUS");
    for (int i = 0; i < 10; i++)
    {
        double ii = static_cast<double>(i);
        EXPECT_EQ(iv[i], i);
        EXPECT_NEAR(fv[i], ii, 1E-15);
        EXPECT_NEAR(imgv[i].real(), ii, 1E-15);
        EXPECT_NEAR(imgv[i].imag(), -ii, 1E-15);
        std::stringstream ss;
        ss << chs << i;
        EXPECT_EQ(chv[i], ss.str());
    }
}

/**
 * bcast_mpi_offset broadcasts from the root of the communicator it is given,
 * not from world rank 0. The MDCell STRU writer relies on this to share the
 * write position within its own (possibly non-world) communicator.
 */
TEST_F(ParaCommon, BcastMPIOffsetOnCommunicator)
{
    const int rank = mpiContext.GetRank();
    const int color = rank % 2;
    MPI_Comm split_comm = MPI_COMM_NULL;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &split_comm);

    int sub_rank = 0;
    MPI_Comm_rank(split_comm, &sub_rank);
    const MPI_Offset expected = static_cast<MPI_Offset>(100 * (color + 1));
    MPI_Offset offset = 0;
    if (sub_rank == 0)
    {
        offset = expected;
    }
    Parallel_Common::bcast_mpi_offset(offset, split_comm);
    EXPECT_EQ(offset, expected);

    MPI_Comm_free(&split_comm);
}

int main(int argc, char** argv)
{

    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    int result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}
#endif
