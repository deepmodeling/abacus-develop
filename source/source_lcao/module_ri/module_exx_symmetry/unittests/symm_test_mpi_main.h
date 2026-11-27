#ifndef SYMM_TEST_MPI_MAIN_H
#define SYMM_TEST_MPI_MAIN_H

// Shared entry point for the symmetry unit tests.
//
// Rule 9 asks code to go through the guarded wrappers rather than calling MPI
// directly. There is no wrapper for process start-up and shutdown, so this
// keeps the two calls in one place here instead of repeating them in every
// test main, and keeps them behind __MPI like the rest of the code base.

#include "gtest/gtest.h"

#ifdef __MPI
#include <mpi.h>
#endif

namespace ModuleSymmetryTest
{
/// Run the registered tests, initialising MPI when it is available.
inline int run_all_tests(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
} // namespace ModuleSymmetryTest

#endif
