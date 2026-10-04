#include <fstream>
#include <sstream>
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../write_pdos_text.h"
#include "source_base/global_variable.h"
#include "source_cell/unitcell.h"

// Build a minimal UnitCell for a single Si atom with one s and one p
// zeta: nwl=1, nw=4, iw2l={0,1,1,1}, iw2m={0,0,1,2} (p: m=0,+1,-1).
// Index arrays (iat2it/iat2ia/iat2iwt/itia2iat) are filled by hand so the
// test only needs the UnitCell/Atom constructor/destructor, not the full
// source_cell library.
static UnitCell* make_si_ucell()
{
    UnitCell* ucell = new UnitCell;
    ucell->ntype = 1;
    ucell->nat = 1;
    ucell->atoms = new Atom[1];
    ucell->set_atom_flag = true;

    Atom& atom = ucell->atoms[0];
    atom.label = "Si";
    atom.na = 1;
    atom.nwl = 1;
    atom.iw2l = {0, 1, 1, 1};
    atom.iw2m = {0, 0, 1, 2};
    atom.nw = 4;
    atom.l_nchi = {1, 1};

    // iat2it / iat2ia
    ucell->iat2it = new int[1];
    ucell->iat2ia = new int[1];
    ucell->iat2it[0] = 0;
    ucell->iat2ia[0] = 0;

    // itia2iat(it, ia) -> iat ; IntArray is (ntype, na)-indexed
    ucell->itia2iat.create(ucell->ntype, atom.na);
    ucell->itia2iat(0, 0) = 0;

    // iat2iwt / npol: first global orbital index of the atom
    ucell->set_iat2iwt_for_test({0}, 1);
    return ucell;
}

class WritePdosTextTest : public ::testing::Test
{
protected:
    UnitCell* ucell = nullptr;

    void SetUp() override
    {
        ucell = make_si_ucell();
    }
    void TearDown() override
    {
        // iat2it/iat2ia are owned by UnitCell's internal Statistics member,
        // whose destructor releases them; do not delete them here.
        delete ucell;
    }

    std::string read_file(const std::string& fname)
    {
        std::ifstream ifs(fname);
        std::stringstream buffer;
        buffer << ifs.rdbuf();
        return buffer.str();
    }
};

// Spin-unpolarized: writes pdoss1_<basis>.txt with s and 3 p columns.
TEST_F(WritePdosTextTest, BasicOutput)
{
    GlobalV::ofs_running.open("test_write_pdos_text.log");

    const int nspin = 1;
    const int nlocal = 4;
    const int npoints = 5;
    const double emin = -1.0;
    const double dos_edelta_ev = 0.5;
    const std::string out_dir = "./";
    const std::string basis = "lcao";
    const int istep = -1;

    // pdos matrix: nlocal x npoints; value = orbital * 100 + point
    ModuleBase::matrix pdos(nlocal, npoints);
    for (int w = 0; w < nlocal; ++w)
    {
        for (int n = 0; n < npoints; ++n)
        {
            pdos(w, n) = w * 100.0 + n;
        }
    }

    ModuleIO::write_pdos_text(*ucell, &pdos, nspin, nlocal, npoints,
                              emin, dos_edelta_ev, out_dir, basis, istep);

    const std::string str = read_file("./pdoss1_lcao.txt");

    // header lines
    EXPECT_THAT(str, testing::HasSubstr("# npoints: 5"));
    EXPECT_THAT(str, testing::HasSubstr("# energy(eV)  atom  species"));
    EXPECT_THAT(str, testing::HasSubstr("s(m=0) p(m=0,+1,-1)"));

    // first data row: energy=-1.0, atom 1, Si, s=0*100+0=0, p columns
    // w=1,2,3 -> 100,200,300
    EXPECT_THAT(str, testing::HasSubstr("   -1.000000    1    Si  0.000000  100.000000  200.000000  300.000000"));
    // last data row: energy=1.0, n=4 -> s=4, p=104,204,304
    EXPECT_THAT(str, testing::HasSubstr("    1.000000    1    Si  4.000000  104.000000  204.000000  304.000000"));

    remove("./pdoss1_lcao.txt");
    GlobalV::ofs_running.close();
    remove("test_write_pdos_text.log");
}

// istep >= 0 appends g{istep+1} to the file name.
TEST_F(WritePdosTextTest, IstepFileName)
{
    GlobalV::ofs_running.open("test_write_pdos_text.log");

    const int nspin = 1;
    const int nlocal = 4;
    const int npoints = 2;
    const double emin = 0.0;
    const double dos_edelta_ev = 1.0;
    ModuleBase::matrix pdos(nlocal, npoints);
    pdos.zero_out();

    ModuleIO::write_pdos_text(*ucell, &pdos, nspin, nlocal, npoints,
                              emin, dos_edelta_ev, "./", "lcao", 2);

    const std::string str = read_file("./pdoss1g3_lcao.txt");
    EXPECT_THAT(str, testing::HasSubstr("# istep: 3"));

    remove("./pdoss1g3_lcao.txt");
    GlobalV::ofs_running.close();
    remove("test_write_pdos_text.log");
}

// Spin-polarized (nspin=2): writes two files, pdos[is] per spin channel.
TEST_F(WritePdosTextTest, SpinPolarized)
{
    GlobalV::ofs_running.open("test_write_pdos_text.log");

    const int nspin = 2;
    const int nlocal = 4;
    const int npoints = 2;
    const double emin = 0.0;
    const double dos_edelta_ev = 1.0;

    ModuleBase::matrix pdos_arr[2];
    pdos_arr[0].create(nlocal, npoints);
    pdos_arr[1].create(nlocal, npoints);
    for (int is = 0; is < 2; ++is)
    {
        for (int w = 0; w < nlocal; ++w)
        {
            for (int n = 0; n < npoints; ++n)
            {
                pdos_arr[is](w, n) = (is + 1) * 1000.0 + w * 100.0 + n;
            }
        }
    }

    ModuleIO::write_pdos_text(*ucell, pdos_arr, nspin, nlocal, npoints,
                              emin, dos_edelta_ev, "./", "lcao", -1);

    const std::string s1 = read_file("./pdoss1_lcao.txt");
    const std::string s2 = read_file("./pdoss2_lcao.txt");
    // spin 1: s orbital w=0, n=0 -> 1000
    EXPECT_THAT(s1, testing::HasSubstr("1000.000000"));
    // spin 2: s orbital w=0, n=0 -> 2000
    EXPECT_THAT(s2, testing::HasSubstr("2000.000000"));

    remove("./pdoss1_lcao.txt");
    remove("./pdoss2_lcao.txt");
    GlobalV::ofs_running.close();
    remove("test_write_pdos_text.log");
}
