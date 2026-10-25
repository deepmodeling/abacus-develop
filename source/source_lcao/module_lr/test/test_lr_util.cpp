#include "source_lcao/module_lr/utils/lr_util_hcontainer.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <complex>
#include <map>
#include <vector>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include "source_base/parallel_common.h"
#include "source_estate/module_pot/td_field_manager.h"

namespace
{
const int test_size = 10;
const int test_nw = 10;
class LRDMTest : public testing::Test
{
protected:
    std::vector<Atom> atoms;
    UnitCell ucell;
    Parallel_Orbitals distribution;
    Parallel_Orbitals* paraV = nullptr;

    void SetUp() override
    {
        atoms.resize(1);
        ucell.atoms = atoms.data(); // Borrowed by the mock UnitCell; vector owns the atoms.
        ucell.ntype = 1;
        ucell.nat = test_size;
        ucell.iat2it.assign(test_size, 0);
        ucell.iat2ia.resize(test_size);
        atoms[0].na = test_size;
        atoms[0].nw = test_nw;
        atoms[0].tau.resize(test_size);
        atoms[0].iw2l.assign(test_nw, 0);
        atoms[0].iw2m.assign(test_nw, 0);
        atoms[0].iw2n.assign(test_nw, 0);
        for (int atom = 0; atom < test_size; ++atom)
        {
            ucell.iat2ia[atom] = atom;
        }
        ucell.set_iat2iwt(1);
        const int dimension = test_size * test_nw;
#ifdef __MPI
        distribution.init(dimension, dimension, 2, MPI_COMM_WORLD);
#else
        distribution.set_serial(dimension, dimension);
#endif
        distribution.set_atomic_trace(ucell.get_iat2iwt(), test_size, dimension);
        paraV = &distribution;
    }
};
}

TEST_F(LRDMTest, LRDotRCountsEachImageOnce)
{
    hamilt::HContainer<double> h1(paraV);
    hamilt::HContainer<double> h2(paraV);
    for (int r = -1; r <= 1; ++r)
    {
        hamilt::AtomPair<double> a(0, 1, r, 0, 0, paraV);
        const int opposite_r = -r;
        hamilt::AtomPair<double> b(0, 1, opposite_r, 0, 0, paraV);
        h1.insert_pair(a);
        h2.insert_pair(b);
    }
    h1.allocate(nullptr, true);
    h2.allocate(nullptr, true);
    double expected = 0.0;
    for (int r = -1; r <= 1; ++r)
    {
        auto& a = h1.find_pair(0, 1)->get_HR_values(r, 0, 0);
        auto& b = h2.find_pair(0, 1)->get_HR_values(r, 0, 0);
        const int size = a.get_row_size() * a.get_col_size();
        for (int i = 0; i < size; ++i)
        {
            const double x = r + 2.0;
            const double y = i + 1.0;
            a.get_pointer()[i] = x;
            b.get_pointer()[i] = y;
            expected += x * y;
        }
    }
    EXPECT_DOUBLE_EQ(LR_Util::dot_R_matrix(h1, h2), expected);
}

TEST_F(LRDMTest, LRTransposePreservesAtomIndices)
{
    const int n = test_size * test_nw;
    const std::vector<ModuleBase::Vector3<double>> kvec(1);
    module_dm::DensityMatrix<double, double> dm(paraV, 1, kvec, 1);
    hamilt::HContainer<double> hr(paraV);
    for (int ia = 0; ia < test_size; ++ia)
    {
        for (int ja = 0; ja < test_size; ++ja)
        {
            hamilt::AtomPair<double> pair(ia, ja, paraV);
            hr.insert_pair(pair);
        }
    }
    hr.allocate(nullptr, true);
    hr.fix_gamma();
    dm.init_dmr(hr);
    auto& dk = dm.get_dmk_vec()[0];
    for (int j = 0; j < paraV->ncol; ++j)
    {
        for (int i = 0; i < paraV->nrow; ++i)
        {
            const int mu = paraV->local2global_row(i);
            const int nu = paraV->local2global_col(j);
            dk[j * paraV->nrow + i] = mu * n + nu;
        }
    }
    dm.cal_dmr(-1);
    const int passes = 2;
    for (int pass = 0; pass < passes; ++pass)
    {
        for (int ip = 0; ip < dm.get_dmr_ptr(1)->size_atom_pairs(); ++ip)
        {
            const auto& pair = dm.get_dmr_ptr(1)->get_atom_pair(ip);
            const auto* values = pair.get_HR_values(0).get_pointer();
            for (int i = 0; i < pair.get_row_size(); ++i)
            {
                for (int j = 0; j < pair.get_col_size(); ++j)
                {
                    const int mu = paraV->local2global_row(pair.get_begin_row() + i);
                    const int nu = paraV->local2global_col(pair.get_begin_col() + j);
                    const double expected = pass == 0 ? mu * n + nu : nu * n + mu;
                    EXPECT_DOUBLE_EQ(values[i * pair.get_col_size() + j], expected);
                }
            }
        }
        if (pass == 0)
        {
            LR_Util::transpose_DMR(dm, *paraV);
        }
    }
}

#ifdef __EXX
namespace RI_2D_Comm
{
// Observe the adapter's spin/k input before the expensive RI redistribution.
template <>
std::vector<std::map<TA, std::map<TAC, RI::Tensor<double>>>>
split_m2D_ktoR<double, std::vector<double>>(
    const UnitCell&, const K_Vectors&, const std::vector<const std::vector<double>*>& mks,
    const Parallel_2D&, const int nspin, const bool)
{
    EXPECT_EQ(nspin, 2);
    EXPECT_EQ(mks.size(), 4U);
    for (std::size_t ik = 0; ik < mks.size(); ++ik)
    {
        EXPECT_DOUBLE_EQ(mks[ik]->front(), ik + 1.0);
    }
    return std::vector<std::map<TA, std::map<TAC, RI::Tensor<double>>>>(nspin);
}
}

TEST_F(LRDMTest, LRExxAdapterIncludesBothSpinChannels)
{
    K_Vectors kv;
    kv.set_nks(4);
    kv.kvec_d.resize(4);
    module_dm::DensityMatrix<double, double> dm(paraV, 2, kv.kvec_d, 2);
    hamilt::HContainer<double> hr(paraV);
    hr.allocate(nullptr, true);
    dm.init_dmr(hr);
    auto& dmk = dm.get_dmk_vec();
    ASSERT_EQ(dmk.size(), 4U);
    for (std::size_t ik = 0; ik < dmk.size(); ++ik)
    {
        std::fill(dmk[ik].begin(), dmk[ik].end(), ik + 1.0);
    }
    const auto ds = LR_Util::get_exx_Ds_gs(dm, ucell, kv, *paraV);
    EXPECT_EQ(ds.size(), 2U);
}
#endif

// The shared R-coordinate implementation also contains an unused TD cleanup path.
int elecstate::TDFieldManager::gauge() const { return 0; }

TEST_F(LRDMTest, SparseWriterKeepsSpinFilesAndDisjointRankKeys)
{
    const int rank = paraV->get_coord_row() * paraV->get_dim1() + paraV->get_coord_col();
    const int process_count = paraV->get_dim0() * paraV->get_dim1();
    const int active_ranks = std::max(1, process_count - 1);
    GlobalV::DRANK = rank; // The common writer uses the existing diagonal-grid output rank.
    hamilt::HContainer<double> structure(paraV);
    const int zero_image = 99;
    hamilt::AtomPair<double> zero(0, 0, zero_image, 0, 0, paraV);
    structure.insert_pair(zero);
    const int image = rank + 1;
    if (rank < active_ranks)
    {
        hamilt::AtomPair<double> pair(0, 0, image, 0, 0, paraV);
        structure.insert_pair(pair);
    }
    structure.allocate(nullptr, true);
    const std::vector<ModuleBase::Vector3<double>> kvec(1);
    module_dm::DensityMatrix<double, double> dm(paraV, 2, kvec, 1);
    dm.init_dmr(structure);
    ASSERT_EQ(dm.get_dmr_vec().size(), 2);
    for (int spin = 0; spin < 2; ++spin)
    {
        if (rank < active_ranks)
        {
            auto& block = dm.get_dmr_vec()[spin]->find_pair(0, 0)->get_HR_values(image, 0, 0);
            const int count = block.get_row_size() * block.get_col_size();
            const double value = (spin + 1) * image;
            std::fill(block.get_pointer(), block.get_pointer() + count, value);
        }
    }
    std::string filename;
    if (rank == 0)
    {
        filename = testing::TempDir() + "lr_sparse_spins_" + std::to_string(getpid());
    }
#ifdef __MPI
    Parallel_Common::bcast_string(filename);
#endif
    const std::string out_dir;
    const double threshold = 1e-10;
    std::ostringstream running_log;
    LR_Util::save_DMR(dm, filename, *paraV, out_dir, threshold, running_log);
    if (rank == 0)
    {
        const int dimension = test_size * test_nw;
        for (int spin = 0; spin < 2; ++spin)
        {
            const std::string spin_file = filename + "_s" + std::to_string(spin);
            std::ifstream input(spin_file);
            ASSERT_TRUE(input.is_open());
            std::string line;
            std::getline(input, line);
            EXPECT_EQ(line, "STEP: 0");
            std::getline(input, line);
            EXPECT_EQ(line, "Matrix Dimension of EDM(R): " + std::to_string(dimension));
            std::getline(input, line);
            EXPECT_EQ(line, "Matrix number of EDM(R): " + std::to_string(active_ranks));
            for (int expected_image = 1; expected_image <= active_ranks; ++expected_image)
            {
                int rx = 0;
                int ry = 0;
                int rz = 0;
                int count = 0;
                input >> rx >> ry >> rz >> count;
                ASSERT_FALSE(input.fail());
                EXPECT_EQ(rx, expected_image);
                EXPECT_EQ(ry, 0);
                EXPECT_EQ(rz, 0);
                ASSERT_GT(count, 0);
                for (int i = 0; i < count; ++i)
                {
                    double value = 0.0;
                    input >> value;
                    const double expected_value = (spin + 1) * expected_image;
                    EXPECT_DOUBLE_EQ(value, expected_value);
                }
                int index = 0;
                for (int i = 0; i < count + dimension + 1; ++i)
                {
                    input >> index;
                }
                ASSERT_FALSE(input.fail());
            }
            input >> std::ws;
            EXPECT_TRUE(input.eof());
            input.close();
            std::remove(spin_file.c_str());
            EXPECT_NE(running_log.str().find(spin_file), std::string::npos);
        }
    }
}
