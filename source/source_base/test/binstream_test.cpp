#include "gtest/gtest.h"
#include "gmock/gmock.h"
#ifdef __unix__
#include <signal.h>
#include <sys/resource.h>
#endif
/************************************************
 *  unit test of binstream.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Binstream()
 *     - Open a binary file
 *   - close()
 *     - Close a binary file
 */

#include "source_base/module_out/binstream.h"

class BinstreamTest : public testing::Test
{
protected:
};

TEST_F(BinstreamTest, variable)
{
    int a1, a2, a3;
    Binstream wfc("wfc", "w");
    EXPECT_EQ(!wfc, false);
    wfc << 10;
    wfc.close();
    EXPECT_EQ(bool(wfc), false);
    wfc.open("wfc", "a");
    EXPECT_EQ(bool(wfc), true);
    wfc << 100;
    wfc.close();
    EXPECT_EQ(!wfc, true);

    wfc.open("wfc", "r");
    EXPECT_EQ(bool(wfc), true);
    wfc >> a1;
    wfc >> a2;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(wfc >> a3, ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("!NOTICE!"));
    EXPECT_EQ(bool(wfc), true);
    wfc.close();
    EXPECT_EQ(bool(wfc), false);

    EXPECT_EQ(a1, 10);
    EXPECT_EQ(a2, 100);

    wfc.open("wfc", "w");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(wfc >> a3, ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("!NOTICE!"));
    remove("wfc");

    wfc.open("wfc", "r");
    EXPECT_EQ(bool(wfc), false); // If file is not open, return false.

    Binstream *p = new Binstream("wfc", "r");
    delete p;

    remove("wfc"); // mohan add 2025-06-22
}

TEST_F(BinstreamTest, array)
{
    int a[10], b[11];
    for(int i = 0; i < 10; ++i)
    {
        a[i] = i;
    }
    Binstream wwfc("wfc", "w");
    wwfc.write(a, 10);
    wwfc.close();

    Binstream rwfc("wfc", "r");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(rwfc.read(b, 11);, ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("!NOTICE!"));
    rwfc.close();
    rwfc.open("wfc", "r");
    rwfc.read(b, 10);
    rwfc.close();
    remove("wfc");

    for(int i = 0; i < 10; ++i)
    {
        EXPECT_EQ(a[i], b[i]);
    }

    wwfc.open("wfc", "w");
}

TEST_F(BinstreamTest, WriteToUnopenedFile)
{
    Binstream ofs;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ofs << 42, ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("!NOTICE!"));

    int a[3] = {1, 2, 3};
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ofs.write(a, 3), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("!NOTICE!"));
}

TEST_F(BinstreamTest, ReadFromUnopenedFile)
{
    Binstream ifs;
    int v = 0;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ifs >> v, ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("!NOTICE!"));

    int a[3] = {0, 0, 0};
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ifs.read(a, 3), ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("!NOTICE!"));
}

// fwrite() may report success while the data is still buffered; with
// RLIMIT_FSIZE=1 the write() call must detect the failure (via fflush)
// instead of returning normally with only one byte on disk.
TEST_F(BinstreamTest, DelayedWriteFailureDetected)
{
#ifdef __unix__
    // Ignore SIGXFSZ so the file-size limit surfaces as EFBIG from
    // fflush/fwrite instead of killing the process.
    struct sigaction old_act;
    ASSERT_EQ(sigaction(SIGXFSZ, NULL, &old_act), 0);
    struct sigaction ign_act = old_act;
    ign_act.sa_handler = SIG_IGN;
    ASSERT_EQ(sigaction(SIGXFSZ, &ign_act, NULL), 0);

    rlimit lim;
    lim.rlim_cur = 1; // allow at most 1 byte per file
    lim.rlim_max = RLIM_INFINITY;
    ASSERT_EQ(setrlimit(RLIMIT_FSIZE, &lim), 0);

    int a[4] = {1, 2, 3, 4};
    EXPECT_EXIT(
        {
            Binstream wfc("wfc_rlimit", "w");
            wfc.write(a, 4);
            wfc.close();
        },
        ::testing::ExitedWithCode(1), "");

    lim.rlim_cur = RLIM_INFINITY;
    ASSERT_EQ(setrlimit(RLIMIT_FSIZE, &lim), 0);
    ASSERT_EQ(sigaction(SIGXFSZ, &old_act, NULL), 0);
    remove("wfc_rlimit");
#endif
}

// A close() on a stream whose earlier buffered writes failed at flush time
// must also be caught (fclose returns EOF).
TEST_F(BinstreamTest, CloseFailureDetected)
{
    // close() on an unopened stream must be a safe no-op.
    Binstream empty;
    empty.close();
    EXPECT_EQ(bool(empty), false);
}
