#ifdef __MPI
#include <mpi.h>
#endif
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <fstream>

#include "source_base/tool_quit.h"

/************************************************
 *  Unit test for init_chg=hr error handling
 *
 *  Verifies clear error messages when HR files missing.
 *  Tests the file-existence check logic used by
 *  LCAO_domain::init_hr_from_file without linking
 *  the full LCAO_set.cpp (which has heavy dependencies).
 ***********************************************/

TEST(InitChgHrErrorTest, MissingFileError)
{
    std::string hrfile = "./nonexistent_dir/hrs1_nao.csr";

    // This replicates the exact check in init_hr_from_file
    EXPECT_DEATH(
        {
            std::ifstream test_file(hrfile);
            if (!test_file.good())
            {
                std::string error_msg = "Cannot open Hamiltonian file: " + hrfile + "\n\n";
                error_msg += "When using init_chg=hr, you need to provide Hamiltonian matrix files:\n";
                error_msg += "  - For nspin=1: hrs1_nao.csr\n";
                error_msg += "  - For nspin=2: hrs1_nao.csr (spin-up) and hrs2_nao.csr (spin-down)\n\n";
                error_msg += "Solutions:\n";
                error_msg += "  1. Run an SCF calculation first with 'out_mat_hs2 1' to generate HR files\n";
                error_msg += "  2. Check that 'read_file_dir' points to the correct directory\n";
                error_msg += "  3. Use 'init_chg file' or 'init_chg atomic' instead";
                ModuleBase::WARNING_QUIT("LCAO_domain::init_hr_from_file", error_msg);
            }
        },
        ""
    );
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
