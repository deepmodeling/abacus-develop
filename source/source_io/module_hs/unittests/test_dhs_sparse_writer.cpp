/**
 * @file test_dhs_sparse_writer.cpp
 * @brief Unit tests for dhs_sparse_writer.cpp (ModuleIO::save_dH_sparse).
 *
 * save_dH_sparse writes the derivative sparse matrices (dHx/dHy/dHz) for
 * either H or S into separate CSR files.  The tests verify that:
 * - text mode respects the sparse threshold
 * - binary mode preserves the native field order
 * - SOC spin mode (nspin == 4) writes the complex variant for all
 *   three directions.
 *
 * dhs_sparse_writer.cpp only reads GlobalV::DRANK; no PARAM state is
 * involved, so the fixture sets DRANK directly.
 */
#include "csr_test_helpers.h"

#include "source_io/module_hs/dhs_sparse_writer.h"

#include <fstream>

// Dummy output stream for running-log writes in tests; discards output.
static std::ofstream ofs_running;

#ifdef __MPI
#include <mpi.h>
#endif

TEST(DhsSparseWriter, TextCountsOnlyValuesAboveThreshold)
{
    remove_derivative_files("h");
    init_sparse_output_globals();

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    UnitCell ucell;
    init_unitcell(ucell);
    LCAO_HS_Arrays arrays;
    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    arrays.all_R_coor.insert(r_vector);
    arrays.dHRx_sparse[0][r_vector][0][0] = 1.0;
    arrays.dHRx_sparse[0][r_vector][0][1] = 1e-12;
    arrays.dHRx_sparse[0][r_vector][1][0] = 0.0;
    arrays.dHRx_sparse[0][r_vector][1][1] = -2.0;

    ModuleIO::save_dH_sparse(5, pv, &ucell, arrays, 1e-10, false, "h", 8, "./", "scf", false, 1, 2, ofs_running);

    const std::string output = read_file("dhrxs1_nao.csr");
    // istep = 5 is the sixth ionic step, printed one-based as 6
    EXPECT_THAT(output, testing::HasSubstr("Ionic_Step 6  # ABACUS version:"));
    EXPECT_THAT(output, testing::HasSubstr(" # print dH/dR matrix in real space dHx(R)\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1 # number of spin directions\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1 # spin index\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 2 # number of localized basis\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1 # number of Bravais lattice vector R\n"));
    EXPECT_THAT(output, testing::HasSubstr("user_defined_lattice\n"));
    EXPECT_THAT(output, testing::HasSubstr("Direct\n 0 0.25 0.5\n"));
    // Values below the threshold are dropped; printToCSR wraps every six entries.
    EXPECT_THAT(output, testing::HasSubstr("0 0 0 2\n # CSR values\n 1.00000000e+00 -2.00000000e+00\n"
                                           " # CSR column indices\n 0 1\n"
                                           " # CSR row pointers\n 0 1 2\n"));

    // Empty blocks still carry the three CSR sections and zero row pointers.
    const std::string y_output = read_file("dhrys1_nao.csr");
    EXPECT_THAT(y_output, testing::HasSubstr("0 0 0 0\n # CSR values\n # CSR column indices\n"
                                             " # CSR row pointers\n 0 0 0\n"));

    remove_derivative_files("h");
}

TEST(DhsSparseWriter, BinaryCountsOnlyValuesAboveThreshold)
{
    remove_derivative_files("h");
    init_sparse_output_globals();

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    UnitCell ucell;
    init_unitcell(ucell);
    LCAO_HS_Arrays arrays;
    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    arrays.all_R_coor.insert(r_vector);
    arrays.dHRx_sparse[0][r_vector][0][0] = 1.0;
    arrays.dHRx_sparse[0][r_vector][0][1] = 1e-12;
    arrays.dHRx_sparse[0][r_vector][1][0] = 0.0;
    arrays.dHRx_sparse[0][r_vector][1][1] = -2.0;

    ModuleIO::save_dH_sparse(6, pv, &ucell, arrays, 1e-10, true, "h", 8, "./", "scf", false, 1, 2, ofs_running);

    std::ifstream ifs("dhrxs1_nao.csr", std::ios::binary);
    ASSERT_TRUE(ifs.is_open());
    // istep = 6 -> one-based step 7 in the binary header
    EXPECT_EQ(read_binary_value<int>(ifs), 7);
    EXPECT_EQ(read_binary_value<int>(ifs), 2);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 2);
    EXPECT_DOUBLE_EQ(read_binary_value<double>(ifs), 1.0);
    EXPECT_DOUBLE_EQ(read_binary_value<double>(ifs), -2.0);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_EQ(read_binary_value<long long>(ifs), 0);
    EXPECT_EQ(read_binary_value<long long>(ifs), 1);
    EXPECT_EQ(read_binary_value<long long>(ifs), 2);

    remove_derivative_files("h");
}

TEST(DhsSparseWriter, SocWritesAllDirections)
{
    remove_derivative_files("s");
    init_sparse_output_globals();

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    UnitCell ucell;
    init_unitcell(ucell);
    LCAO_HS_Arrays arrays;
    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    arrays.all_R_coor.insert(r_vector);
    arrays.dHRx_soc_sparse[r_vector][0][0] = std::complex<double>(1.0, 0.0);
    arrays.dHRy_soc_sparse[r_vector][0][1] = std::complex<double>(2.0, -1.0);
    arrays.dHRz_soc_sparse[r_vector][1][1] = std::complex<double>(-3.0, 0.5);

    ModuleIO::save_dH_sparse(7, pv, &ucell, arrays, 1e-10, false, "s", 8, "./", "scf", false, 4, 2, ofs_running);

    const std::string x_output = read_file("dsrx_nao.csr");
    const std::string y_output = read_file("dsry_nao.csr");
    const std::string z_output = read_file("dsrz_nao.csr");
    EXPECT_THAT(x_output, testing::HasSubstr(" # print dS/dR matrix in real space dSx(R)\n"));
    EXPECT_THAT(y_output, testing::HasSubstr(" # print dS/dR matrix in real space dSy(R)\n"));
    EXPECT_THAT(z_output, testing::HasSubstr(" # print dS/dR matrix in real space dSz(R)\n"));
    EXPECT_THAT(x_output, testing::HasSubstr(" 1 # number of Bravais lattice vector R\n"));
    EXPECT_THAT(y_output, testing::HasSubstr(" 1 # number of Bravais lattice vector R\n"));
    EXPECT_THAT(z_output, testing::HasSubstr(" 1 # number of Bravais lattice vector R\n"));
    EXPECT_THAT(x_output, testing::HasSubstr("0 0 0 1\n # CSR values\n"));
    EXPECT_THAT(y_output, testing::HasSubstr("0 0 0 1\n # CSR values\n"));
    EXPECT_THAT(z_output, testing::HasSubstr("0 0 0 1\n # CSR values\n"));
    EXPECT_THAT(x_output, testing::HasSubstr("(1.00000000e+00,0.00000000e+00)"));
    EXPECT_THAT(y_output, testing::HasSubstr("(2.00000000e+00,-1.00000000e+00)"));
    EXPECT_THAT(z_output, testing::HasSubstr("(-3.00000000e+00,5.00000000e-01)"));

    remove_derivative_files("s");
}

TEST(DhsSparseWriter, OverlapDerivativeWritesSingleSpinChannel)
{
    remove_derivative_files("s");
    remove_derivative_files("h");
    init_sparse_output_globals();

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    UnitCell ucell;
    init_unitcell(ucell);
    LCAO_HS_Arrays arrays;
    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    arrays.all_R_coor.insert(r_vector);
    arrays.dHRx_sparse[0][r_vector][0][0] = 1.0;

    // dS/dR is spin-independent, so nspin = 2 must write a single file whose
    // name carries no spin index.
    ModuleIO::save_dH_sparse(4, pv, &ucell, arrays, 1e-10, false, "s", 8, "./", "scf", false, 2, 2, ofs_running);

    std::ifstream spin_file("dsrx_nao.csr");
    EXPECT_TRUE(spin_file.is_open());
    spin_file.close();
    std::ifstream legacy_spin1_file("dsrxs1_nao.csr");
    EXPECT_FALSE(legacy_spin1_file.is_open());
    legacy_spin1_file.close();

    const std::string output = read_file("dsrx_nao.csr");
    EXPECT_THAT(output, testing::HasSubstr("Ionic_Step 5  # ABACUS version:"));
    // The header reports a single spin direction even though nspin = 2.
    EXPECT_THAT(output, testing::HasSubstr(" 1 # number of spin directions\n"));
    EXPECT_THAT(output, testing::HasSubstr("0 0 0 1\n # CSR values\n 1.00000000e+00"));

    // dH/dR stays spin-dependent and still writes both spin channels.
    arrays.dHRx_sparse[1][r_vector][1][1] = -2.0;
    ModuleIO::save_dH_sparse(4, pv, &ucell, arrays, 1e-10, false, "h", 8, "./", "scf", false, 2, 2, ofs_running);
    std::ifstream h_spin2_file("dhrxs2_nao.csr");
    EXPECT_TRUE(h_spin2_file.is_open());
    h_spin2_file.close();

    remove_derivative_files("s");
    remove_derivative_files("h");
}

TEST(DhsSparseWriter, AppendModeReusesFixedFilename)
{
    remove_derivative_files("s");
    init_sparse_output_globals();

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    UnitCell ucell;
    init_unitcell(ucell);
    LCAO_HS_Arrays arrays;
    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    arrays.all_R_coor.insert(r_vector);
    arrays.dHRx_sparse[0][r_vector][0][0] = 1.0;

    // out_app_flag = 1 reuses the fixed filename without a g{N} suffix; the
    // first written ionic step (istep == 0) truncates it.
    ModuleIO::save_dH_sparse(0, pv, &ucell, arrays, 1e-10, false, "s", 8, "./", "md", 1, 1, 2, ofs_running);
    const std::string first_frame = read_file("dsrx_nao.csr");
    ASSERT_FALSE(first_frame.empty());
    EXPECT_EQ(count_substr(first_frame, "Ionic_Step "), 1);
    EXPECT_EQ(count_substr(first_frame, "Ionic_Step 1"), 1);

    // The next ionic step appends a second frame to the same file.
    ModuleIO::save_dH_sparse(1, pv, &ucell, arrays, 1e-10, false, "s", 8, "./", "md", 1, 1, 2, ofs_running);
    const std::string two_frames = read_file("dsrx_nao.csr");
    EXPECT_EQ(count_substr(two_frames, "Ionic_Step "), 2);
    EXPECT_EQ(count_substr(two_frames, "Ionic_Step 1"), 1);
    EXPECT_EQ(count_substr(two_frames, "Ionic_Step 2"), 1);

    // Mode 2 keeps one suffixed file per ionic step and never appends; mode 0
    // overwrites the fixed filename instead.
    ModuleIO::save_dH_sparse(0, pv, &ucell, arrays, 1e-10, false, "s", 8, "./", "md", 2, 1, 2, ofs_running);
    const std::vector<std::string> mode2_lines = read_lines("dsrxg1_nao.csr");
    ASSERT_FALSE(mode2_lines.empty());
    EXPECT_EQ(mode2_lines[0].compare(0, 12, "Ionic_Step 1"), 0);

    remove_derivative_files("s");
    remove_derivative_files("s", 0);
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
#endif

    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
