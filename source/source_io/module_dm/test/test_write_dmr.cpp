/**
 * @file test_write_dmr.cpp
 * @brief Unit tests for write_dmr.cpp: dmr_gen_fname and write_dmr_csr.
 *
 * dmr_gen_fname builds the spin-dependent DM(R) filename; write_dmr_csr
 * writes the DM(R) matrix in the legacy CSR text layout.  Both tests pin
 * the current contract so format changes are caught by a failing test.
 * Neither function reads PARAM (the high-level write_dmr driver does, but
 * it is not exercised here), so no PARAM fixture is required.
 */
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "source_io/module_dm/write_dmr.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/module_hcontainer/atom_pair.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

namespace
{
std::string read_file(const std::string& filename)
{
    std::ifstream ifs(filename.c_str());
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return oss.str();
}

void init_unitcell(UnitCell& ucell)
{
    ucell.latName = "user_defined_lattice";
    ucell.lat0 = 10.0;
    ucell.latvec.e11 = 1.0;
    ucell.latvec.e22 = 1.0;
    ucell.latvec.e33 = 1.0;
    ucell.ntype = 1;
    ucell.nat = 1;
    ucell.atoms = new Atom[1];
    ucell.set_atom_flag = true;
    ucell.atoms[0].label = "Si";
    ucell.atoms[0].na = 1;
    ucell.atoms[0].nw = 2;
    ucell.atoms[0].taud.resize(1);
    ucell.atoms[0].taud[0] = ModuleBase::Vector3<double>(0.0, 0.25, 0.5);
}

void init_serial_orbitals(Parallel_Orbitals& pv)
{
    pv.atom_begin_row.resize(2);
    pv.atom_begin_col.resize(2);
    pv.atom_begin_row[0] = 0;
    pv.atom_begin_row[1] = 2;
    pv.atom_begin_col[0] = 0;
    pv.atom_begin_col[1] = 2;
    pv.nrow = 2;
    pv.ncol = 2;
    pv.set_serial(2, 2);
}

void fill_matrix(hamilt::HContainer<double>& matrix, Parallel_Orbitals& pv, double* values)
{
    hamilt::AtomPair<double> pair(0, 0, 0, 0, 0, &pv, values);
    matrix.insert_pair(pair);
}
} // namespace

TEST(WriteDmr, GenFnameKeepsCurrentContract)
{
    EXPECT_EQ(ModuleIO::dmr_gen_fname(1, 0, true, -1), "dmrs1_nao.csr");
    EXPECT_EQ(ModuleIO::dmr_gen_fname(1, 1, false, 2), "dmrs2g3_nao.csr");
    EXPECT_EQ(ModuleIO::dmr_gen_fname(2, 0, true, 5), "dmrs1_nao.npz");
}

TEST(WriteDmr, CsrHeaderKeepsCurrentFormat)
{
    std::string filename = "write_hs_r_header_dmr.csr";
    std::remove(filename.c_str());

    UnitCell ucell;
    init_unitcell(ucell);
    Parallel_Orbitals pv;
    init_serial_orbitals(pv);
    hamilt::HContainer<double> matrix(&pv);
    double values[4] = {0.25, 0.0, 0.0, 0.75};
    fill_matrix(matrix, pv, values);

    ModuleIO::write_dmr_csr(filename, &ucell, 4, &matrix, 0, 1, 2);

    const std::string output = read_file(filename);
    EXPECT_THAT(output, testing::HasSubstr(" --- Ionic Step 1 ---\n"));
    EXPECT_THAT(output, testing::HasSubstr(" # print density matrix in real space DM(R)\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 2 # number of spin directions\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 2 # spin index\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 2 # number of localized basis\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1 # number of Bravais lattice vector R\n"));

    std::remove(filename.c_str());
}
