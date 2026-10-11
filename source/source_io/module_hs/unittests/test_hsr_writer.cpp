/**
 * @file test_hsr_writer.cpp
 * @brief Unit tests for hsr_writer.cpp: filename-generation helpers and the
 *        CSR text/binary H(R), S(R) writer paths (write_hcontainer_csr,
 *        write_hcontainer_csr_binary, write_hsr).
 *
 * The filename helpers are pure string functions with no physics state.
 * The writer tests pin down the exact CSR text/binary layout on disk so a
 * format change is caught by a failing test instead of a broken reader.
 */
#include "csr_test_helpers.h"

#include "source_base/module_out/csr_reader.h"
#include "source_base/constants.h"
#include "source_estate/fp_energy.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
#include "source_io/module_hs/hsr_writer.h"

#include <algorithm>
#include <complex>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

namespace
{
/**
 * @brief Helper: compare two filenames with a clear failure message.
 */
void expect_fname_eq(const std::string& expected, const std::string& actual)
{
    EXPECT_EQ(expected, actual) << "Expected \"" << expected << "\", got \"" << actual << "\"";
}
} // namespace

// ---------------------------------------------------------------------------
// hsr_gen_fname: H(R) spin-dependent filename
// ---------------------------------------------------------------------------

TEST(HsrWriterFname, HsrOverwriteWithStep)
{
    // prefix="data-", ispin=0, append=false, istep=0
    // -> "data-1g1_nao.csr"
    const std::string fname = ModuleIO::hsr_gen_fname("data-", 0, false, 0);
    expect_fname_eq("data-1g1_nao.csr", fname);
}

TEST(HsrWriterFname, HsrOverwriteWithStepTwo)
{
    // istep=2 -> g3
    const std::string fname = ModuleIO::hsr_gen_fname("data-", 1, false, 2);
    expect_fname_eq("data-2g3_nao.csr", fname);
}

TEST(HsrWriterFname, HsrAppendSuppressesStep)
{
    // append=true -> step number is omitted
    const std::string fname = ModuleIO::hsr_gen_fname("data-", 0, true, 5);
    expect_fname_eq("data-1_nao.csr", fname);
}

TEST(HsrWriterFname, HsrNegativeStepSuppressesStep)
{
    // istep=-1 -> step number is omitted (no "g" suffix)
    const std::string fname = ModuleIO::hsr_gen_fname("data-", 0, false, -1);
    expect_fname_eq("data-1_nao.csr", fname);
}

TEST(HsrWriterFname, HsrDatExtension)
{
    // out_type=2 -> ".dat" extension
    const std::string fname = ModuleIO::hsr_gen_fname("data-", 0, false, 0, 2);
    expect_fname_eq("data-1g1_nao.dat", fname);
}

TEST(HsrWriterFname, HsrCsrExtensionDefault)
{
    // out_type=1 (default) -> ".csr" extension
    const std::string fname = ModuleIO::hsr_gen_fname("data-", 0, false, 0, 1);
    expect_fname_eq("data-1g1_nao.csr", fname);
}

// ---------------------------------------------------------------------------
// sr_gen_fname: S(R) spin-independent filename
// ---------------------------------------------------------------------------

TEST(HsrWriterFname, SrOverwriteWithStep)
{
    // append=false, istep=0 -> "srg1_nao.csr"
    const std::string fname = ModuleIO::sr_gen_fname(false, 0);
    expect_fname_eq("srg1_nao.csr", fname);
}

TEST(HsrWriterFname, SrAppendSuppressesStep)
{
    // append=true -> "sr_nao.csr" (no step number)
    const std::string fname = ModuleIO::sr_gen_fname(true, 5);
    expect_fname_eq("sr_nao.csr", fname);
}

TEST(HsrWriterFname, SrNegativeStepSuppressesStep)
{
    // istep=-1 -> "sr_nao.csr"
    const std::string fname = ModuleIO::sr_gen_fname(false, -1);
    expect_fname_eq("sr_nao.csr", fname);
}

TEST(HsrWriterFname, SrDatExtension)
{
    // out_type=2 -> ".dat"
    const std::string fname = ModuleIO::sr_gen_fname(false, 0, 2);
    expect_fname_eq("srg1_nao.dat", fname);
}

// ---------------------------------------------------------------------------
// dhr_gen_fname: dH/dR spin-dependent filename
// ---------------------------------------------------------------------------

TEST(HsrWriterFname, DhrOverwriteWithStep)
{
    // prefix="H", ispin=0, append=false, istep=0
    // -> "Hrs1g1_nao.csr"
    const std::string fname = ModuleIO::dhr_gen_fname("H", 0, false, 0);
    expect_fname_eq("Hrs1g1_nao.csr", fname);
}

TEST(HsrWriterFname, DhrAppendSuppressesStep)
{
    // append=true -> step omitted
    const std::string fname = ModuleIO::dhr_gen_fname("H", 1, true, 3);
    expect_fname_eq("Hrs2_nao.csr", fname);
}

TEST(HsrWriterFname, DhrNegativeStepSuppressesStep)
{
    // istep=-1 -> step omitted
    const std::string fname = ModuleIO::dhr_gen_fname("dH", 0, false, -1);
    expect_fname_eq("dHrs1_nao.csr", fname);
}

// ---------------------------------------------------------------------------
// write_hcontainer_csr: CSR text format layout
// ---------------------------------------------------------------------------

TEST(HsrWriterIo, HContainerCsrHeaderKeepsCurrentFormat)
{
    const std::string filename = "write_hs_r_header_h.csr";
    std::remove(filename.c_str());

    UnitCell ucell;
    init_unitcell(ucell);
    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    hamilt::HContainer<double> matrix(&pv);
    double values[4] = {1.0, 0.0, 0.5, 2.0};
    fill_matrix(matrix, pv, values);

    // label "H": header carries the Fermi energy of this spin channel
    const double efermi_eV = 5.4321;
    ModuleIO::write_hcontainer_csr(filename, &ucell, 5, &matrix, 0, 0, 1, "H", "", efermi_eV, true);

    const std::string output = read_file(filename);
    EXPECT_THAT(output, testing::HasSubstr(" --- Ionic Step 1 ---\n"));
    EXPECT_THAT(output, testing::HasSubstr(" # print H matrix in real space H(R)\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1 # number of spin directions\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1 # spin index, E_Fermi = 5.4321 eV\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 2 # number of localized basis\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1 # number of Bravais lattice vector R\n"));
    EXPECT_THAT(output, testing::HasSubstr(" user_defined_lattice\n"));
    EXPECT_THAT(output, testing::HasSubstr(" Si\n"));
    EXPECT_THAT(output, testing::HasSubstr(" Direct\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 0 0.25 0.5\n"));
    EXPECT_THAT(output, testing::HasSubstr(" #                               CSR Format                             #\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 0 0 0 3\n"));
    EXPECT_THAT(output, testing::HasSubstr(" # CSR values\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1.00000e+00 5.00000e-01 2.00000e+00"));
    EXPECT_THAT(output, testing::Not(testing::HasSubstr("# representation:")));

    std::remove(filename.c_str());
}

TEST(HsrWriterIo, GammaFoldedHeaderKeepsCsrReadable)
{
    const std::string filename = "write_hs_r_gamma_folded.csr";
    const std::string representation_note
        = "gamma-only folded matrix; stored R-space contributions are summed into R = (0, 0, 0)";
    std::remove(filename.c_str());

    UnitCell ucell;
    init_unitcell(ucell);
    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    hamilt::HContainer<double> matrix(&pv);
    double values[4] = {1.0, 0.0, 0.5, 2.0};
    fill_matrix(matrix, pv, values);

    // label "H" without Fermi energy available: pass has_efermi = false
    const double no_efermi = 0.0;
    ModuleIO::write_hcontainer_csr(
        filename, &ucell, 5, &matrix, 0, 0, 1, "H", representation_note, no_efermi, false);

    const std::string output = read_file(filename);
    EXPECT_THAT(output, testing::HasSubstr("# representation: " + representation_note + "\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1 # number of Bravais lattice vector R\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 0 0 0 3\n"));
    // has_efermi = false: no Fermi annotation even though label is "H"
    EXPECT_THAT(output, testing::Not(testing::HasSubstr("E_Fermi")));

    ModuleIO::csrFileReader<double> reader(filename);
    ASSERT_EQ(reader.getNumberOfR(), 1);
    EXPECT_EQ(reader.getMatrixDimension(), 2);
    EXPECT_EQ(reader.getRCoordinate(0), std::vector<int>({0, 0, 0}));

    std::remove(filename.c_str());
}

TEST(HsrWriterIo, HContainerCsrAppendKeepsCurrentStepSections)
{
    const std::string filename = "write_hs_r_append_s.csr";
    std::remove(filename.c_str());

    UnitCell ucell;
    init_unitcell(ucell);
    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    hamilt::HContainer<double> matrix(&pv);
    double values[4] = {1.0, 0.0, 0.0, 1.0};
    fill_matrix(matrix, pv, values);

    // label "S": header carries no Fermi energy (has_efermi = false)
    const double ignored_efermi = 0.0;
    ModuleIO::write_hcontainer_csr(filename, &ucell, 4, &matrix, 0, 0, 1, "S", "", ignored_efermi, false);
    ModuleIO::write_hcontainer_csr(filename, &ucell, 4, &matrix, 1, 0, 1, "S", "", ignored_efermi, false);

    const std::string output = read_file(filename);
    EXPECT_EQ(count_substr(output, " --- Ionic Step "), 2);
    EXPECT_THAT(output, testing::HasSubstr(" --- Ionic Step 1 ---\n"));
    EXPECT_THAT(output, testing::HasSubstr(" --- Ionic Step 2 ---\n"));
    EXPECT_EQ(count_substr(output, " # print S matrix in real space S(R)\n"), 2);
    // S(R) files must not contain a Fermi energy annotation
    EXPECT_THAT(output, testing::Not(testing::HasSubstr("E_Fermi")));

    std::remove(filename.c_str());
}

TEST(HsrWriterIo, HContainerCsrHeaderCarriesPerChannelFermi)
{
    UnitCell ucell;
    init_unitcell(ucell);
    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    hamilt::HContainer<double> matrix(&pv);
    double values[4] = {1.0, 0.0, 0.5, 2.0};
    fill_matrix(matrix, pv, values);

    // Each spin channel is written to its own file (write mode, istep=0)
    const double efermi_up_eV = 5.4321;
    const double efermi_dw_eV = 3.2109;
    const std::string filename_up = "write_hs_r_header_fermi_up.csr";
    const std::string filename_dw = "write_hs_r_header_fermi_dw.csr";
    std::remove(filename_up.c_str());
    std::remove(filename_dw.c_str());
    ModuleIO::write_hcontainer_csr(filename_up, &ucell, 5, &matrix, 0, 0, 2, "H", "", efermi_up_eV, true);
    ModuleIO::write_hcontainer_csr(filename_dw, &ucell, 5, &matrix, 0, 1, 2, "H", "", efermi_dw_eV, true);

    const std::string output_up = read_file(filename_up);
    const std::string output_dw = read_file(filename_dw);
    EXPECT_THAT(output_up, testing::HasSubstr(" 1 # spin index, E_Fermi = 5.4321 eV\n"));
    EXPECT_THAT(output_dw, testing::HasSubstr(" 2 # spin index, E_Fermi = 3.2109 eV\n"));

    std::remove(filename_up.c_str());
    std::remove(filename_dw.c_str());
}

// ---------------------------------------------------------------------------
// write_hcontainer_csr_binary: native binary layout
// ---------------------------------------------------------------------------

TEST(HsrWriterIo, HContainerBinaryWritesSortedAndEmptyRBlocks)
{
    const std::string filename = "write_hs_r_native_binary_real.dat";
    std::remove(filename.c_str());

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    hamilt::HContainer<double> matrix(&pv);
    double empty_values[4] = {1e-12, 0.0, 0.0, -1e-12};
    double values[4] = {1.0, 1e-12, 0.0, -2.0};
    fill_matrix_at_R(matrix, pv, 1, 0, 0, empty_values);
    fill_matrix_at_R(matrix, pv, -1, 0, 0, values);

    ModuleIO::write_hcontainer_csr_binary(filename, &matrix, -1, false);

    std::ifstream ifs(filename.c_str(), std::ios::binary);
    ASSERT_TRUE(ifs.is_open());
    const NativeDoubleRecord record = read_native_double_record(ifs);
    EXPECT_EQ(record.step, 0);
    EXPECT_EQ(record.nbasis, 2);
    ASSERT_EQ(record.blocks.size(), 2);
    EXPECT_EQ(record.blocks[0].rx, -1);
    EXPECT_EQ(record.blocks[0].ry, 0);
    EXPECT_EQ(record.blocks[0].rz, 0);
    EXPECT_THAT(record.blocks[0].values, testing::ElementsAre(1.0, -2.0));
    EXPECT_THAT(record.blocks[0].columns, testing::ElementsAre(0, 1));
    EXPECT_THAT(record.blocks[0].row_ptr, testing::ElementsAre(0, 1, 2));
    EXPECT_EQ(record.blocks[1].rx, 1);
    EXPECT_EQ(record.blocks[1].ry, 0);
    EXPECT_EQ(record.blocks[1].rz, 0);
    EXPECT_TRUE(record.blocks[1].values.empty());
    EXPECT_TRUE(record.blocks[1].columns.empty());
    EXPECT_THAT(record.blocks[1].row_ptr, testing::ElementsAre(0, 0, 0));
    EXPECT_EQ(ifs.peek(), std::ifstream::traits_type::eof());

    std::remove(filename.c_str());
}

TEST(HsrWriterIo, HContainerBinaryWritesComplexValuesAsDoublePairs)
{
    const std::string filename = "write_hs_r_native_binary_complex.dat";
    std::remove(filename.c_str());

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    hamilt::HContainer<std::complex<double>> matrix(&pv);
    std::complex<double> values[4] = {
        std::complex<double>(1.0, 2.0),
        std::complex<double>(0.0, 0.0),
        std::complex<double>(0.0, 0.0),
        std::complex<double>(-3.0, 4.0),
    };
    fill_matrix_at_R(matrix, pv, 0, 0, 0, values);

    ModuleIO::write_hcontainer_csr_binary(filename, &matrix, 4, false);

    std::ifstream ifs(filename.c_str(), std::ios::binary);
    ASSERT_TRUE(ifs.is_open());
    EXPECT_EQ(read_binary_value<int>(ifs), 4);
    EXPECT_EQ(read_binary_value<int>(ifs), 2);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 2);
    EXPECT_DOUBLE_EQ(read_binary_value<double>(ifs), 1.0);
    EXPECT_DOUBLE_EQ(read_binary_value<double>(ifs), 2.0);
    EXPECT_DOUBLE_EQ(read_binary_value<double>(ifs), -3.0);
    EXPECT_DOUBLE_EQ(read_binary_value<double>(ifs), 4.0);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_EQ(read_binary_value<long long>(ifs), 0);
    EXPECT_EQ(read_binary_value<long long>(ifs), 1);
    EXPECT_EQ(read_binary_value<long long>(ifs), 2);
    EXPECT_EQ(ifs.peek(), std::ifstream::traits_type::eof());

    std::remove(filename.c_str());
}

TEST(HsrWriterIo, HContainerBinaryAppendsAndCanOverwrite)
{
    const std::string filename = "write_hs_r_native_binary_append.dat";
    std::remove(filename.c_str());

    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    hamilt::HContainer<double> matrix(&pv);
    double values[4] = {1.0, 0.0, 0.0, 2.0};
    fill_matrix(matrix, pv, values);

    ModuleIO::write_hcontainer_csr_binary(filename, &matrix, 0, true);
    ModuleIO::write_hcontainer_csr_binary(filename, &matrix, 1, true);

    std::ifstream appended(filename.c_str(), std::ios::binary);
    ASSERT_TRUE(appended.is_open());
    EXPECT_EQ(read_native_double_record(appended).step, 0);
    EXPECT_EQ(read_native_double_record(appended).step, 1);
    EXPECT_EQ(appended.peek(), std::ifstream::traits_type::eof());
    appended.close();

    ModuleIO::write_hcontainer_csr_binary(filename, &matrix, 2, false);
    std::ifstream replaced(filename.c_str(), std::ios::binary);
    ASSERT_TRUE(replaced.is_open());
    EXPECT_EQ(read_native_double_record(replaced).step, 2);
    EXPECT_EQ(replaced.peek(), std::ifstream::traits_type::eof());

    std::remove(filename.c_str());
}

TEST(HsrWriterIo, HContainerBinaryMpiGatherWritesCompleteFiles)
{
#ifndef __MPI
    GTEST_SKIP() << "MPI support is required.";
#else
    int mpi_size = 0;
    int mpi_rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_size != 2)
    {
        GTEST_SKIP() << "This test requires exactly two MPI ranks.";
    }

    const std::string hr_filename = "hrs1_nao.dat";
    const std::string sr_filename = "sr_nao.dat";
    if (mpi_rank == 0)
    {
        std::remove(hr_filename.c_str());
        std::remove(sr_filename.c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);

    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = 1;
    ucell.atoms = new Atom[1];
    ucell.set_atom_flag = true;
    ucell.atoms[0].na = 1;
    ucell.atoms[0].nw = 2;
    ucell.iat2it.resize(1);
    ucell.iat2it[0] = 0;
    ucell.set_iat2iwt(1);
    const int* iat2iwt = ucell.get_iat2iwt();
    Parallel_Orbitals serial_pv;
    serial_pv.set_serial(2, 2);
    serial_pv.set_atomic_trace(iat2iwt, 1, 2);
    Parallel_Orbitals parallel_pv;
    ASSERT_EQ(parallel_pv.init(2, 2, 1, MPI_COMM_WORLD), 0);
    parallel_pv.set_atomic_trace(iat2iwt, 1, 2);

    hamilt::HContainer<double> hr_serial(ucell, &serial_pv);
    hamilt::HContainer<double> sr_serial(ucell, &serial_pv);
    double hr_values[4] = {1.0, 0.5, 0.0, 2.0};
    double sr_values[4] = {1.0, 0.0, 0.0, 1.0};
    if (mpi_rank == 0)
    {
        std::copy(hr_values, hr_values + 4, hr_serial.get_atom_pair(0).get_pointer(0));
        std::copy(sr_values, sr_values + 4, sr_serial.get_atom_pair(0).get_pointer(0));
    }

    hamilt::HContainer<double> hr_parallel(ucell, &parallel_pv);
    hamilt::HContainer<double> sr_parallel(ucell, &parallel_pv);
    hamilt::transferSerial2Parallels(hr_serial, &hr_parallel, 0);
    hamilt::transferSerial2Parallels(sr_serial, &sr_parallel, 0);

    init_sparse_output_globals();
    std::vector<hamilt::HContainer<double>*> hr_vec(1, &hr_parallel);
    elecstate::Efermi eferm;
    eferm.two_efermi = false;
    eferm.ef = 0.4; // Ry; used only for the H(R) header
    std::ofstream ofs_running_null; // not opened; test does not inspect the running log
    ModuleIO::write_hsr(
        hr_vec, &sr_parallel, &ucell, 2, 8, parallel_pv, true, true, iat2iwt, 1, 0, "./", eferm,
        ofs_running_null);
    MPI_Barrier(MPI_COMM_WORLD);

    if (mpi_rank == 0)
    {
        std::ifstream hr_stream(hr_filename.c_str(), std::ios::binary);
        ASSERT_TRUE(hr_stream.is_open());
        const NativeDoubleRecord hr_record = read_native_double_record(hr_stream);
        ASSERT_EQ(hr_record.blocks.size(), 1);
        EXPECT_THAT(hr_record.blocks[0].values, testing::ElementsAre(1.0, 0.5, 2.0));
        EXPECT_THAT(hr_record.blocks[0].columns, testing::ElementsAre(0, 1, 1));
        EXPECT_THAT(hr_record.blocks[0].row_ptr, testing::ElementsAre(0, 2, 3));

        std::ifstream sr_stream(sr_filename.c_str(), std::ios::binary);
        ASSERT_TRUE(sr_stream.is_open());
        const NativeDoubleRecord sr_record = read_native_double_record(sr_stream);
        ASSERT_EQ(sr_record.blocks.size(), 1);
        EXPECT_THAT(sr_record.blocks[0].values, testing::ElementsAre(1.0, 1.0));
        EXPECT_THAT(sr_record.blocks[0].columns, testing::ElementsAre(0, 1));
        EXPECT_THAT(sr_record.blocks[0].row_ptr, testing::ElementsAre(0, 1, 2));

        std::remove(hr_filename.c_str());
        std::remove(sr_filename.c_str());
    }
    delete[] ucell.atoms;
    ucell.atoms = nullptr;
    ucell.set_atom_flag = false;
#endif
}

// ---------------------------------------------------------------------------
// write_hsr: text CSR (out_type=1) threads per-spin Fermi energy from Efermi
// ---------------------------------------------------------------------------

TEST(HsrWriterIo, WriteHsrTextCsrCarriesPerSpinFermiFromEfermi)
{
    // Exercise the fix path inside write_hsr: eferm.get_efval(ispin) * Ry_to_eV
    // for a two-Fermi (nspin=2) case. The only existing write_hsr test uses
    // out_type=2 (binary), which skips the Fermi branch entirely.
    // append=true suppresses the ionic-step suffix in the generated names.
    const std::string hr_up_filename = "hrs1_nao.csr";
    const std::string hr_dw_filename = "hrs2_nao.csr";
    const std::string sr_filename = "sr_nao.csr";
    std::remove(hr_up_filename.c_str());
    std::remove(hr_dw_filename.c_str());
    std::remove(sr_filename.c_str());

    UnitCell ucell;
    init_unitcell(ucell);
    // The __MPI path in write_hsr calls set_atomic_trace and gatherParallels,
    // which require a real atom-to-orbital map; nullptr with nat=0 throws on
    // that path.
    ucell.iat2it.resize(1);
    ucell.iat2it[0] = 0;
    ucell.set_iat2iwt(1);
    const int* iat2iwt_ptr = ucell.get_iat2iwt();
    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    // gatherParallels reaches get_indexes_row(iat) on this source pv, which
    // dereferences iat2iwt_; init_serial_orbitals does not set it, so the
    // atomic trace must be installed explicitly.
    pv.set_atomic_trace(iat2iwt_ptr, 1, 2);

    hamilt::HContainer<double> hr_up(&pv);
    hamilt::HContainer<double> hr_dw(&pv);
    hamilt::HContainer<double> sr(&pv);
    double values[4] = {1.0, 0.0, 0.5, 2.0};
    double sr_values[4] = {1.0, 0.0, 0.0, 1.0};
    fill_matrix(hr_up, pv, values);
    fill_matrix(hr_dw, pv, values);
    fill_matrix(sr, pv, sr_values);

    init_sparse_output_globals();

    elecstate::Efermi eferm;
    eferm.two_efermi = true;
    eferm.ef_up = 0.5;   // Ry
    eferm.ef_dw = 0.3;   // Ry
    const double ef_up_eV = eferm.ef_up * ModuleBase::Ry_to_eV;   // 6.802849
    const double ef_dw_eV = eferm.ef_dw * ModuleBase::Ry_to_eV;   // 4.0817094

    std::vector<hamilt::HContainer<double>*> hr_vec;
    hr_vec.push_back(&hr_up);
    hr_vec.push_back(&hr_dw);

    std::ofstream ofs_running_null; // not opened; test does not inspect the running log
    ModuleIO::write_hsr(
        hr_vec, &sr, &ucell, 1, 8, pv, true, true, iat2iwt_ptr, 1, 0, "./", eferm,
        ofs_running_null);

    const std::string output_up = read_file(hr_up_filename);
    const std::string output_dw = read_file(hr_dw_filename);

    // write_hsr multiplies each spin's Fermi by Ry_to_eV and passes has_efermi=true
    std::ostringstream expected_up;
    expected_up << " 1 # spin index, E_Fermi = " << std::setprecision(6) << ef_up_eV << " eV\n";
    std::ostringstream expected_dw;
    expected_dw << " 2 # spin index, E_Fermi = " << std::setprecision(6) << ef_dw_eV << " eV\n";
    EXPECT_THAT(output_up, testing::HasSubstr(expected_up.str()));
    EXPECT_THAT(output_dw, testing::HasSubstr(expected_dw.str()));

    std::remove(hr_up_filename.c_str());
    std::remove(hr_dw_filename.c_str());
    std::remove(sr_filename.c_str());
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
