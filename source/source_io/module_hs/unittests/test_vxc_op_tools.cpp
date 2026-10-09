/**
 * @file test_vxc_op_tools.cpp
 * @brief Unit tests for vxc_op_tools.cpp helper functions
 *
 * Tests the functions that do not require a full Parallel_Orbitals
 * / ScaLAPACK context:
 * - get_real (double and complex overloads)
 * - set_para2d_MO (non-MPI set_serial path)
 * - orbital_energy (serial Parallel_2D, double and complex)
 * - all_band_energy (serial Parallel_2D, double and complex)
 *
 * cVc is not tested here because it requires pv.desc / pv.desc_wfc
 * (full Parallel_Orbitals with BLACS context).
 * write_orb_energy is not tested here because it requires K_Vectors
 * (vtable pulls in klist.cpp -> symmetry.cpp dependency chain).
 */
#include <gtest/gtest.h>

#include "source_io/module_hs/vxc_op_tools.h"
#include "source_base/parallel_2d.h"
#include "source_base/parallel_comm.h"
#include "source_base/matrix.h"

#include <complex>
#include <vector>
#ifdef __MPI
#include <mpi.h>
#endif

namespace
{
/**
 * @brief Create a serial Parallel_2D that owns all rows/cols on rank 0.
 */
Parallel_2D make_serial_pv(int dim)
{
    Parallel_2D pv;
    pv.set_serial(dim, dim);
    return pv;
}
} // namespace

// ---------------------------------------------------------------------------
// get_real: extract real part
// ---------------------------------------------------------------------------

TEST(VxcOpTools, GetRealDouble)
{
    const double d = 3.14;
    EXPECT_DOUBLE_EQ(3.14, ModuleIO::get_real(d));
}

TEST(VxcOpTools, GetRealComplex)
{
    const std::complex<double> c(1.5, 2.5);
    EXPECT_DOUBLE_EQ(1.5, ModuleIO::get_real(c));
}

// ---------------------------------------------------------------------------
// set_para2d_MO: set up MO-space 2D distribution
// ---------------------------------------------------------------------------

TEST(VxcOpTools, SetPara2dMO)
{
    const int nbands = 4;
    Parallel_Orbitals pv;
#ifdef __MPI
    ASSERT_EQ(pv.init(nbands, nbands, 1, DIAG_WORLD), 0);
#else
    pv.set_serial(nbands, nbands);
#endif
    Parallel_2D p2d;
    ModuleIO::set_para2d_MO(pv, nbands, p2d);
    EXPECT_EQ(p2d.get_global_row_size(), nbands);
    EXPECT_EQ(p2d.get_global_col_size(), nbands);
}

// ---------------------------------------------------------------------------
// orbital_energy: extract diagonal from MO matrix
// ---------------------------------------------------------------------------

TEST(VxcOpTools, OrbitalEnergyDouble)
{
    const int nbands = 3;
    Parallel_2D p2d = make_serial_pv(nbands);
    // MO matrix (row-major, p2d local): 3x3 diagonal
    // mat_mo[i * row_size + j] where row_size = get_row_size()
    const int row_size = p2d.get_row_size();
    std::vector<double> mat_mo(row_size * p2d.get_col_size(), 0.0);
    // Set diagonal: 1.0, 2.0, 3.0
    mat_mo[0 * row_size + 0] = 1.0;
    mat_mo[1 * row_size + 1] = 2.0;
    mat_mo[2 * row_size + 2] = 3.0;

    const std::vector<double> e = ModuleIO::orbital_energy(0, nbands, mat_mo, p2d);
    ASSERT_EQ(e.size(), static_cast<size_t>(nbands));
    EXPECT_DOUBLE_EQ(1.0, e[0]);
    EXPECT_DOUBLE_EQ(2.0, e[1]);
    EXPECT_DOUBLE_EQ(3.0, e[2]);
}

TEST(VxcOpTools, OrbitalEnergyComplex)
{
    const int nbands = 2;
    Parallel_2D p2d = make_serial_pv(nbands);
    const int row_size = p2d.get_row_size();
    std::vector<std::complex<double>> mat_mo(row_size * p2d.get_col_size());
    mat_mo[0 * row_size + 0] = std::complex<double>(1.5, 0.0);
    mat_mo[1 * row_size + 1] = std::complex<double>(3.5, 0.0);

    const std::vector<double> e = ModuleIO::orbital_energy(0, nbands, mat_mo, p2d);
    ASSERT_EQ(e.size(), static_cast<size_t>(nbands));
    EXPECT_DOUBLE_EQ(1.5, e[0]);
    EXPECT_DOUBLE_EQ(3.5, e[1]);
}

// ---------------------------------------------------------------------------
// all_band_energy: weighted sum of diagonal
// ---------------------------------------------------------------------------

TEST(VxcOpTools, AllBandEnergyDouble)
{
    const int nbands = 3;
    const int ik = 0;
    Parallel_2D p2d = make_serial_pv(nbands);
    const int row_size = p2d.get_row_size();
    std::vector<double> mat_mo(row_size * p2d.get_col_size(), 0.0);
    mat_mo[0 * row_size + 0] = 1.0;
    mat_mo[1 * row_size + 1] = 2.0;
    mat_mo[2 * row_size + 2] = 3.0;

    // wg: occupation weights (nks x nbands)
    ModuleBase::matrix wg(1, nbands);
    wg(0, 0) = 1.0;
    wg(0, 1) = 0.5;
    wg(0, 2) = 0.0;

    const double e = ModuleIO::all_band_energy(ik, mat_mo, p2d, wg);
    // 1.0*1.0 + 2.0*0.5 + 3.0*0.0 = 2.0
    EXPECT_NEAR(2.0, e, 1e-12);
}

TEST(VxcOpTools, AllBandEnergyComplex)
{
    const int nbands = 2;
    const int ik = 0;
    Parallel_2D p2d = make_serial_pv(nbands);
    const int row_size = p2d.get_row_size();
    std::vector<std::complex<double>> mat_mo(row_size * p2d.get_col_size());
    mat_mo[0 * row_size + 0] = std::complex<double>(2.0, 0.0);
    mat_mo[1 * row_size + 1] = std::complex<double>(4.0, 0.0);

    ModuleBase::matrix wg(1, nbands);
    wg(0, 0) = 1.0;
    wg(0, 1) = 0.25;

    const double e = ModuleIO::all_band_energy(ik, mat_mo, p2d, wg);
    // 2.0*1.0 + 4.0*0.25 = 3.0
    EXPECT_NEAR(3.0, e, 1e-12);
}

// ---------------------------------------------------------------------------
// TODO: cVc — requires full Parallel_Orbitals with BLACS context (pv.desc,
// pv.desc_wfc) for the ScaLAPACK gemm path. To test, either:
// (a) refactor cVc to accept raw matrix pointers + Parallel_2D instead of
//     Parallel_Orbitals (governance rule 10: pass explicit arguments), or
// (b) construct a real Parallel_Orbitals via init(dim, dim, 1, DIAG_WORLD)
//     and fill desc/desc_wfc — but this also requires valid BLACS context
//     and ScaLAPACK linkage.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// TODO: write_orb_energy — requires K_Vectors whose vtable (in klist.cpp +
// reciprocal_grid.cpp) pulls in symmetry.cpp -> module_cell/module_symmetry
// dependency chain. To test, either:
// (a) refactor write_orb_energy to accept int nks instead of const
//     K_Vectors& (governance rule 10), or
// (b) add a test-only K_Vectors lightweight subclass that stubs get_nks().
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Main: MPI_Init for Parallel_Reduce::reduce_all in orbital_energy / all_band_energy
// ---------------------------------------------------------------------------

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    DIAG_WORLD = MPI_COMM_WORLD;
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
