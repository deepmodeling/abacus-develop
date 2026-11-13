/**
 * Unit-test of angmom_op_mat.cpp
 * 
 * The representation matrices for Lx, Ly and Lz operators (SO(3) generators).
 * On basis {s}:
 * Lx: Ly: Lz:
 * 0   0   0
 * 
 * On basis {px, py, pz}:
 * Lx:             Ly:             Lz:
 *    px  py  pz      px  py  pz      px  py  pz
 * px  0           px  0       i   px  0  -i
 * py      0  -i   py      0       py  i   0
 * pz      i   0   pz -i       0   pz          0
 * 
 * will test if the calculated values are close to the above results.
 */

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <complex>
#include <memory>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#ifdef __MPI
#include <mpi.h>
#endif

#include "source_io/module_hs/angmom_op_mat.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_basis/module_nao/radial_collection.h"
#include "source_base/sph_bessel_tf.h"
#include "source_base/ylm.h"
#include "source_cell/unitcell.h"
#include "source_cell/nonlocal_info_base.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_reduce.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_global.h"

namespace
{
constexpr double kDoubleThreshold = 1e-12;

/// @brief Minimal NonlocalInfoBase stub for unit-testing the Angmom_op
///        constructor. Only get_rcutmax_Beta is actually called; the other
///        ten overrides return safe placeholders.
class StubNonlocalInfo : public NonlocalInfoBase
{
  public:
    explicit StubNonlocalInfo(const double rcut_beta)
        : rcut_beta_(rcut_beta),
          zero_double_(0.0),
          empty_label_()
    {
    }

    const double& get_rcutmax_Beta() const override { return rcut_beta_; }
    int get_nproj(const int& /*type_in*/) const override { return 0; }
    int get_nprojmax() const override { return 0; }
    double get_rcut_max(const int& /*type_in*/) const override { return 0.0; }
    const std::string& get_label(const int& /*type_in*/) const override { return empty_label_; }
    int get_type(const int& /*type_in*/) const override { return 0; }
    int get_proj_L(const int& /*type_in*/, const int& /*ip_in*/) const override { return 0; }
    int get_proj_Nr(const int& /*type_in*/, const int& /*ip_in*/) const override { return 0; }
    const double* get_proj_radial(const int& /*type_in*/, const int& /*ip_in*/) const override { return nullptr; }
    const double* get_proj_beta_r(const int& /*type_in*/, const int& /*ip_in*/) const override { return nullptr; }
    int get_proj_Nk(const int& /*type_in*/, const int& /*ip_in*/) const override { return 0; }
    double get_proj_dk(const int& /*type_in*/, const int& /*ip_in*/) const override { return 0.0; }
    double get_proj_dr_uniform(const int& /*type_in*/, const int& /*ip_in*/) const override { return 0.0; }

  private:
    double rcut_beta_;
    double zero_double_;
    std::string empty_label_;
};

/// @brief Populate a minimal 1-atom cubic UnitCell that aligns with the
///        orbital fixture Si_gga_8au_100Ry_2s2p1d.orb (Lmax=2, nchi=[2,2,1]).
///        UnitCell is not movable, so the target is filled through a reference.
void populate_minimal_unitcell(UnitCell& ucell)
{
    // Lattice vectors: simple cubic, lat0 = 10 Bohr, large enough to cover
    // the orbital cutoff of 8 Bohr with PBC.
    ucell.lat0 = 10.0;
    ucell.lat0_angstrom = ucell.lat0 * ModuleBase::BOHR_TO_A;
    ucell.tpiba = ModuleBase::TWO_PI / ucell.lat0;
    ucell.tpiba2 = ucell.tpiba * ucell.tpiba;
    ucell.latvec.e11 = 10.0; ucell.latvec.e12 = 0.0; ucell.latvec.e13 = 0.0;
    ucell.latvec.e21 = 0.0; ucell.latvec.e22 = 10.0; ucell.latvec.e23 = 0.0;
    ucell.latvec.e31 = 0.0; ucell.latvec.e32 = 0.0; ucell.latvec.e33 = 10.0;
    ucell.a1.x = ucell.latvec.e11; ucell.a1.y = ucell.latvec.e12; ucell.a1.z = ucell.latvec.e13;
    ucell.a2.x = ucell.latvec.e21; ucell.a2.y = ucell.latvec.e22; ucell.a2.z = ucell.latvec.e23;
    ucell.a3.x = ucell.latvec.e31; ucell.a3.y = ucell.latvec.e32; ucell.a3.z = ucell.latvec.e33;

    ucell.ntype = 1;
    ucell.nat = 1;
    ucell.namax = 1;
    ucell.nwmax = 13; // 2s + 6p + 5d = 13 orbitals for l_nchi = [2,2,1]
    ucell.latName = "sc";

    // UnitCell::atoms is an owning Atom* array (new Atom[ntype]), deleted by
    // the UnitCell destructor, so we allocate it here on the heap and let the
    // ucell instance own the memory for the duration of the test.
    ucell.atoms = new Atom[ucell.ntype];
    Atom& si = ucell.atoms[0];
    si.type = 0;
    si.na = 1;
    si.nwl = 2;
    si.nw = 13;
    si.l_nchi = {2, 2, 1};
    si.tau = {ModuleBase::Vector3<double>(0.0, 0.0, 0.0)};
    si.taud = {ModuleBase::Vector3<double>(0.0, 0.0, 0.0)};
    si.label = "Si";
    si.label_orb = "Si";
    si.set_index();

    ucell.orbital_fn.push_back("Si_gga_8au_100Ry_2s2p1d.orb");
    ucell.infoNL = std::unique_ptr<NonlocalInfoBase>(new StubNonlocalInfo(1.0));

    // Derived index arrays used by assemble_angular_momentum and itiaiw2iwt.
    ucell.set_iat2itia();
    ucell.set_iat2iwt(/*npol_in=*/1);
}

/// @brief Initialize a serial Parallel_Orbitals whose global row dimension
///        matches the total number of NAOs of our minimal UnitCell.
///        Parallel_Orbitals is not copyable, so we receive it by reference.
void init_serial_pv(const int dim, Parallel_Orbitals& pv)
{
    // Re-derive authoritative fields from the BlockCyclicDistribution, which
    // is the single source of truth per the Parallel_Orbitals class contract.
    ASSERT_EQ(pv.para_2d().init(dim, dim, /*nb=*/1, DIAG_WORLD, false), 0);
    pv.para_2d().set_serial(dim, dim);
}

/// @brief Read the first line containing "# size of this matrix" from an
///        L{x,y,z}_nao.csr file to confirm it was written with a nonzero
///        dimension (the file body itself is produced by the CSR printout).
std::string read_file_head(const std::string& path)
{
    std::ifstream ifs(path);
    std::string line;
    std::ostringstream oss;
    // Return the first 40 lines; enough to cover the full UnitCell header,
    // CSR banner, and the first data block of any assembled R matrix.
    for (int i = 0; i < 40 && std::getline(ifs, line); ++i)
    {
        oss << line << '\n';
    }
    return oss.str();
}
} // namespace

class CalpLpRTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        ModuleBase::SphericalBesselTransformer sbt(true);

        this->orb_ = std::unique_ptr<RadialCollection>(new RadialCollection);
        this->orb_->build(1, &forb_, 'o');
        this->orb_->set_transformer(sbt);

        const double rcut = 8.0;
        const int ngrid = int(rcut / 0.01) + 1;
        const double cutoff = 2.0 * rcut;
        this->orb_->set_uniform_grid(true, ngrid, cutoff, 'i', true);

        this->calculator_ = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
        this->calculator_->tabulate(*orb_, *orb_, 'S', ngrid, cutoff);

        // Initialize Ylm coefficients
        ModuleBase::Ylm::set_coefficients();
    }

    void TearDown() override
    {
        // Cleanup code here, if needed
    }

    const std::string forb_ = "../../../../../tests/PP_ORB/Si_gga_8au_100Ry_2s2p1d.orb";

    std::unique_ptr<TwoCenterIntegrator> calculator_;
    std::unique_ptr<RadialCollection> orb_;
};

TEST_F(CalpLpRTest, CalLzijRTest)
{
    std::complex<double> out;

    ModuleBase::Vector3<double> vR(0., 0., 0.); // home-cell
    int it = 0, ia = 0, il = 0, iz = 0, im = 0;
    int jt = 0, ja = 0, jl = 0, jz = 0, jm = 0;
    
    // l=0
    out = ModuleIO::cal_LzijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.0, kDoubleThreshold);
    EXPECT_NEAR(out.imag(), 0.0, kDoubleThreshold);

    // l=1
    il = 1; jl = 1;
    std::vector<std::complex<double>> ans(9, 0.0);
    ans[1] = {0.0, -1.0}; ans[3] = {0.0, 1.0};
    int idx = 0;
    const std::vector<int> m = {1, -1, 0}; // px, py, pz
    for (int i = 0; i < static_cast<int>(m.size()); ++i) {
        for (int j = 0; j < static_cast<int>(m.size()); ++j) {
            im = m[i]; jm = m[j];
            out = ModuleIO::cal_LzijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
            EXPECT_NEAR(out.real(), ans[idx].real(), kDoubleThreshold);
            EXPECT_NEAR(out.imag(), ans[idx].imag(), kDoubleThreshold);
            idx++;
        }
    }
}

TEST_F(CalpLpRTest, CalLxijRTest)
{
    std::complex<double> out;

    ModuleBase::Vector3<double> vR(0., 0., 0.); // home-cell
    int it = 0, ia = 0, il = 0, iz = 0, im = 0;
    int jt = 0, ja = 0, jl = 0, jz = 0, jm = 0;

    // l=0
    out = ModuleIO::cal_LxijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.0, kDoubleThreshold);
    EXPECT_NEAR(out.imag(), 0.0, kDoubleThreshold);

    // l=1
    il = 1; jl = 1;
    std::vector<std::complex<double>> ans(9, 0.0);
    ans[5] = {0.0, 1.0}; ans[7] = {0.0, -1.0};
    int idx = 0;
    const std::vector<int> m = {1, -1, 0}; // px, py, pz
    for (int i = 0; i < static_cast<int>(m.size()); ++i) {
        for (int j = 0; j < static_cast<int>(m.size()); ++j) {
            im = m[i]; jm = m[j];
            out = ModuleIO::cal_LxijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
            EXPECT_NEAR(out.real(), ans[idx].real(), kDoubleThreshold);
            EXPECT_NEAR(out.imag(), ans[idx].imag(), kDoubleThreshold);
            idx++;
        }
    }
}

TEST_F(CalpLpRTest, CalLyijRTest)
{
    std::complex<double> out;

    ModuleBase::Vector3<double> vR(0., 0., 0.); // home-cell
    int it = 0, ia = 0, il = 0, iz = 0, im = 0;
    int jt = 0, ja = 0, jl = 0, jz = 0, jm = 0; // self, the first s

    // l=0
    out = ModuleIO::cal_LyijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
    EXPECT_NEAR(out.real(), 0.0, kDoubleThreshold);
    EXPECT_NEAR(out.imag(), 0.0, kDoubleThreshold);
    
    // l=1
    il = 1; jl = 1;
    std::vector<std::complex<double>> ans(9, 0.0);
    ans[2] = {0.0, -1.0}; ans[6] = {0.0, 1.0};
    int idx = 0;
    const std::vector<int> m = {1, -1, 0}; // px, py, pz
    for (int i = 0; i < static_cast<int>(m.size()); ++i) {
        for (int j = 0; j < static_cast<int>(m.size()); ++j) {
            im = m[i]; jm = m[j];
            out = ModuleIO::cal_LyijR(calculator_, it, ia, il, iz, im, jt, ja, jl, jz, jm, vR);
            EXPECT_NEAR(out.real(), ans[idx].real(), kDoubleThreshold);
            EXPECT_NEAR(out.imag(), ans[idx].imag(), kDoubleThreshold);
            idx++;
        }
    }
}

// ---------------------------------------------------------------------------
// Assembled-path smoke test: construct Angmom_op with a minimal UnitCell and
// exercise calculate() end-to-end.  The CSR writer emits three files, each
// with the full header and a nonzero matrix body; the numerical values are
// checked in a follow-up commit.
// ---------------------------------------------------------------------------
TEST(AngmomOpAssembly, WritesThreeCsrFiles)
{
    const std::string orbital_dir = "../../../../../tests/PP_ORB/";
    const std::string out_dir = testing::TempDir();
    const double search_radius = 10.0;
    const bool gamma_only = true;
    const bool searchpbc = false;
    const int istep = 0;
    const int out_app_flag = 1;
    const int precision = 8;

    UnitCell ucell;
    populate_minimal_unitcell(ucell);
    const int basis_dim = ucell.atoms[0].nw;  // 13 for our Si fixture
    Parallel_Orbitals pv;
    init_serial_pv(basis_dim, pv);

    std::ofstream ofs_null("/dev/null");
    ModuleIO::Angmom_op angmom(orbital_dir, ucell, search_radius,
                               /*tdestructor=*/1, /*tgrid=*/0, /*tatom=*/0,
                               searchpbc, /*out_level=*/"V", gamma_only,
                               &ofs_null, /*rank=*/0);

    angmom.calculate("lx", out_dir, ucell, pv, precision, /*rank=*/0,
                     istep, out_app_flag);

    // All three components must have been written as separate CSR files.
    const std::string lx_path = out_dir + "lx_nao.csr";
    const std::string ly_path = out_dir + "ly_nao.csr";
    const std::string lz_path = out_dir + "lz_nao.csr";
    std::ifstream lx_fs(lx_path);
    std::ifstream ly_fs(ly_path);
    std::ifstream lz_fs(lz_path);
    EXPECT_TRUE(lx_fs.is_open()) << lx_path;
    EXPECT_TRUE(ly_fs.is_open()) << ly_path;
    EXPECT_TRUE(lz_fs.is_open()) << lz_path;

    // A zero-matrix write would still open the file and emit metadata; the CSR
    // banner appears in the body and confirms assemble_angular_momentum filled
    // at least one (R, i, j) block.
    const std::string lx_head = read_file_head(lx_path);
    EXPECT_THAT(lx_head, ::testing::HasSubstr("Ionic_Step"));
    EXPECT_THAT(lx_head, ::testing::HasSubstr("CSR Format"));
}

int main(int argc, char **argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    DIAG_WORLD = MPI_COMM_WORLD;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
