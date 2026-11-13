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

/// Elements with a magnitude below this value are dropped by the assembly
/// layer; the expected matrix below applies the same cutoff.
constexpr double kLThreshold = 1e-10;

/// Tolerance for comparing against CSR values, which the writer prints with
/// precision 8 (eight digits after the decimal point).
constexpr double kCsrValueTolerance = 1e-6;

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

/// @brief Read the leading lines of an L{x,y,z}_nao.csr file so the test can
///        check the provenance line and the CSR banner.
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

/// @brief One Bravais-vector block read back from a *_nao.csr file.
struct CsrBlock
{
    int rx = 0;
    int ry = 0;
    int rz = 0;
    int nnz = 0;
    std::vector<double> values;
    std::vector<int> columns;
    std::vector<int> row_ptrs;
};

/// @brief Accept only a line made of exactly four integers, which is the
///        "Rx Ry Rz nnz" header of one R block.  Every other non-comment line
///        is either shorter (UnitCell header) or longer (data).
bool parse_block_header(const std::string& line, CsrBlock& block)
{
    std::istringstream iss(line);
    int rx = 0;
    int ry = 0;
    int rz = 0;
    int nnz = 0;
    if (!(iss >> rx >> ry >> rz >> nnz))
    {
        return false;
    }
    std::string trailing;
    if (iss >> trailing)
    {
        return false;
    }
    block.rx = rx;
    block.ry = ry;
    block.rz = rz;
    block.nnz = nnz;
    return true;
}

/// @brief Parse every R block of an L{x,y,z}_nao.csr file.  The writer emits
///        "Rx Ry Rz nnz", then " # CSR values", " # CSR column indices" and
///        " # CSR row pointers" sections, each of which may wrap over lines.
std::vector<CsrBlock> read_csr_blocks(const std::string& path)
{
    enum Section
    {
        kHeader,
        kValues,
        kColumns,
        kRowPtrs
    };

    std::ifstream ifs(path);
    std::vector<CsrBlock> blocks;
    std::string line;
    Section section = kHeader;
    CsrBlock current;
    int dim = 0;

    while (std::getline(ifs, line))
    {
        if (line.find("# number of localized basis") != std::string::npos)
        {
            std::istringstream iss(line);
            iss >> dim;
            continue;
        }
        if (line.find("# CSR values") != std::string::npos)
        {
            section = kValues;
            continue;
        }
        if (line.find("# CSR column indices") != std::string::npos)
        {
            section = kColumns;
            continue;
        }
        if (line.find("# CSR row pointers") != std::string::npos)
        {
            section = kRowPtrs;
            continue;
        }
        if (line.find('#') != std::string::npos)
        {
            continue; // provenance, banner and count comments
        }
        if (line.find_first_not_of(" \t\r") == std::string::npos)
        {
            continue; // blank separator between sections
        }

        if (section == kHeader)
        {
            parse_block_header(line, current);
            continue;
        }

        std::istringstream iss(line);
        if (section == kValues)
        {
            double value = 0.0;
            while (iss >> value)
            {
                current.values.push_back(value);
            }
            if (static_cast<int>(current.values.size()) == current.nnz)
            {
                section = kHeader;
            }
        }
        else if (section == kColumns)
        {
            int column = 0;
            while (iss >> column)
            {
                current.columns.push_back(column);
            }
            if (static_cast<int>(current.columns.size()) == current.nnz)
            {
                section = kHeader;
            }
        }
        else
        {
            int row_ptr = 0;
            while (iss >> row_ptr)
            {
                current.row_ptrs.push_back(row_ptr);
            }
            if (static_cast<int>(current.row_ptrs.size()) == dim + 1)
            {
                blocks.push_back(current);
                current = CsrBlock();
                section = kHeader;
            }
        }
    }
    return blocks;
}

/// @brief Dense lookup inside one CSR block; missing entries read as zero.
double csr_value(const CsrBlock& block, const int row, const int col)
{
    const int begin = block.row_ptrs[row];
    const int end = block.row_ptrs[row + 1];
    for (int k = begin; k < end; ++k)
    {
        if (block.columns[k] == col)
        {
            return block.values[k];
        }
    }
    return 0.0;
}

/// @brief Position of a signed magnetic quantum number inside its (2l+1)
///        block, using the ABACUS real-harmonic order 0, +1, -1, +2, -2, ...
int m_position(const int m)
{
    if (m > 0)
    {
        return 2 * m - 1;
    }
    if (m < 0)
    {
        return -2 * m;
    }
    return 0;
}

/// @brief First local orbital index of every angular-momentum block of one
///        atom.  Re-derived here from l_nchi instead of sharing the
///        production helper, so a change of the production mapping is caught.
std::vector<int> block_offsets(const Atom& atom)
{
    std::vector<int> offsets(atom.nwl + 1, 0);
    for (int l = 0; l < atom.nwl; ++l)
    {
        const int block_size = atom.l_nchi[l] * (2 * l + 1);
        offsets[l + 1] = offsets[l] + block_size;
    }
    return offsets;
}

/// @brief Home-cell L matrix evaluated directly with the Lx/Ly/Lz kernels.
///        The one-atom cell makes every pair a self pair, so the displacement
///        handed to the kernels is zero.
std::vector<std::vector<double>> expected_l_matrix(
    const std::unique_ptr<TwoCenterIntegrator>& calculator,
    const Atom& atom,
    const char dir)
{
    const std::vector<int> offsets = block_offsets(atom);
    std::vector<std::vector<double>> mat(atom.nw, std::vector<double>(atom.nw, 0.0));
    const ModuleBase::Vector3<double> dtau(0.0, 0.0, 0.0);

    for (int li = 0; li <= atom.nwl; ++li)
    {
        const int block_size_i = 2 * li + 1;
        for (int iz = 0; iz < atom.l_nchi[li]; ++iz)
        {
            for (int mi = -li; mi <= li; ++mi)
            {
                const int iw_i = offsets[li] + iz * block_size_i + m_position(mi);
                for (int lj = 0; lj <= atom.nwl; ++lj)
                {
                    const int block_size_j = 2 * lj + 1;
                    for (int jz = 0; jz < atom.l_nchi[lj]; ++jz)
                    {
                        for (int mj = -lj; mj <= lj; ++mj)
                        {
                            std::complex<double> val(0.0, 0.0);
                            if (dir == 'x')
                            {
                                val = ModuleIO::cal_LxijR(calculator,
                                    0, 0, li, iz, mi, 0, 0, lj, jz, mj, dtau);
                            }
                            else if (dir == 'y')
                            {
                                val = ModuleIO::cal_LyijR(calculator,
                                    0, 0, li, iz, mi, 0, 0, lj, jz, mj, dtau);
                            }
                            else
                            {
                                val = ModuleIO::cal_LzijR(calculator,
                                    0, 0, li, iz, mi, 0, 0, lj, jz, mj, dtau);
                            }

                            if (std::fabs(val.imag()) > kLThreshold)
                            {
                                const int iw_j = offsets[lj] + jz * block_size_j + m_position(mj);
                                mat[iw_i][iw_j] = val.imag();
                            }
                        }
                    }
                }
            }
        }
    }
    return mat;
}

/// @brief Compare a written CSR block against the expected dense matrix.
void expect_block_matches(const CsrBlock& block, const std::vector<std::vector<double>>& expected)
{
    const int dim = static_cast<int>(expected.size());
    ASSERT_EQ(static_cast<int>(block.row_ptrs.size()), dim + 1);
    ASSERT_EQ(static_cast<int>(block.columns.size()), block.nnz);
    ASSERT_EQ(static_cast<int>(block.values.size()), block.nnz);

    for (int row = 0; row < dim; ++row)
    {
        const int begin = block.row_ptrs[row];
        const int end = block.row_ptrs[row + 1];
        ASSERT_GE(begin, 0);
        ASSERT_LE(begin, end);
        ASSERT_LE(end, block.nnz);

        std::vector<double> dense(dim, 0.0);
        for (int k = begin; k < end; ++k)
        {
            const int col = block.columns[k];
            ASSERT_GE(col, 0);
            ASSERT_LT(col, dim);
            dense[col] = block.values[k];
        }
        for (int col = 0; col < dim; ++col)
        {
            EXPECT_NEAR(dense[col], expected[row][col], kCsrValueTolerance)
                << "CSR row " << row << ", column " << col;
        }
    }
}

/// @brief Anchor the p block of Lz to the analytic 3x3 matrix documented at
///        the top of this file: px and py couple as -i and +i while pz is
///        inert.  The writer stores the imaginary part, so -1 and +1 appear.
void expect_lz_p_block(const CsrBlock& block, const Atom& atom)
{
    const std::vector<int> offsets = block_offsets(atom);
    const int p_start = offsets[1];
    const int p_block_size = 3;
    const int px_position = m_position(1);
    const int py_position = m_position(-1);

    for (int shell = 0; shell < atom.l_nchi[1]; ++shell)
    {
        const int px = p_start + shell * p_block_size + px_position;
        const int py = p_start + shell * p_block_size + py_position;
        EXPECT_NEAR(csr_value(block, px, py), -1.0, kCsrValueTolerance);
        EXPECT_NEAR(csr_value(block, py, px), 1.0, kCsrValueTolerance);
        // A purely imaginary Hermitian block has a vanishing diagonal.
        EXPECT_NEAR(csr_value(block, px, px), 0.0, kCsrValueTolerance);
        EXPECT_NEAR(csr_value(block, py, py), 0.0, kCsrValueTolerance);
    }
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
// Assembled-path test: construct Angmom_op with a minimal UnitCell and
// exercise calculate() end-to-end.  The writer emits one CSR block per
// Bravais vector for each component; every block is read back and compared
// with the L kernels evaluated on an independently tabulated two-center
// table, so index mapping, thresholding and the CSR layout are all covered.
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

    // The provenance line is one-based (Ionic_Step 1 for istep = 0) and the
    // body must carry the CSR banner that precedes the matrix blocks.
    const std::string lx_head = read_file_head(lx_path);
    EXPECT_THAT(lx_head, ::testing::HasSubstr("Ionic_Step"));
    EXPECT_THAT(lx_head, ::testing::HasSubstr("CSR Format"));

    // Independent oracle: tabulate the same two-center table the constructor
    // builds, then evaluate the kernels directly.  The assembled values are
    // checked against the kernels, not against hard-coded numbers.
    ModuleBase::SphericalBesselTransformer sbt(true);
    const std::string oracle_orbital = orbital_dir + ucell.orbital_fn[0];
    RadialCollection oracle_orbs;
    oracle_orbs.build(1, &oracle_orbital, 'o');
    oracle_orbs.set_transformer(sbt);
    const double orbital_rcut = oracle_orbs.rcut_max();
    const int oracle_ngrid = int(orbital_rcut / 0.01) + 1;
    const double oracle_cutoff = 2.0 * orbital_rcut;
    oracle_orbs.set_uniform_grid(true, oracle_ngrid, oracle_cutoff, 'i', true);
    std::unique_ptr<TwoCenterIntegrator> oracle_calculator(new TwoCenterIntegrator);
    oracle_calculator->tabulate(oracle_orbs, oracle_orbs, 'S', oracle_ngrid, oracle_cutoff);

    const std::vector<std::string> paths = {lx_path, ly_path, lz_path};
    const std::vector<char> components = {'x', 'y', 'z'};
    std::vector<std::vector<CsrBlock>> all_blocks;
    for (size_t ic = 0; ic < paths.size(); ++ic)
    {
        const std::vector<CsrBlock> blocks = read_csr_blocks(paths[ic]);
        // The single atom has no periodic images inside the search radius,
        // so exactly the home-cell block must be present.
        ASSERT_EQ(blocks.size(), 1u) << paths[ic];
        EXPECT_EQ(blocks[0].rx, 0);
        EXPECT_EQ(blocks[0].ry, 0);
        EXPECT_EQ(blocks[0].rz, 0);
        EXPECT_GT(blocks[0].nnz, 0) << paths[ic];

        const std::vector<std::vector<double>> expected =
            expected_l_matrix(oracle_calculator, ucell.atoms[0], components[ic]);
        expect_block_matches(blocks[0], expected);
        all_blocks.push_back(blocks);
    }

    // Analytic anchor for the p block of Lz, as documented at the top of this
    // file: <px|Lz|py> = -i and <py|Lz|px> = +i, independent of the oracle.
    expect_lz_p_block(all_blocks[2][0], ucell.atoms[0]);
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
