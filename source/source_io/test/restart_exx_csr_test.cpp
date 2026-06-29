#include "gtest/gtest.h"

#include "source_io/module_restart/restart_exx_csr.h"

#include <map>

namespace
{
void init_unitcell_for_ri(UnitCell& ucell)
{
    ucell.ntype = 1;
    ucell.nat = 1;
    ucell.atoms = new Atom[1];
    ucell.set_atom_flag = true;
    ucell.atoms[0].na = 1;
    ucell.atoms[0].nw = 2;
    ucell.atoms[0].stapos_wf = 0;
    ucell.iat2it = new int[1]{0};
    ucell.iat2ia = new int[1]{0};
}
} // namespace

TEST(RestartExxCsr, CalculateRITensorSparseDropsBelowThresholdEntries)
{
    UnitCell ucell;
    init_unitcell_for_ri(ucell);

    RI::Tensor<double> matrix({2, 2});
    matrix(0, 0) = 1.0;
    matrix(0, 1) = 1e-12;
    matrix(1, 0) = 0.0;
    matrix(1, 1) = -2.0;

    std::map<int, std::map<ModuleIO::TAC, RI::Tensor<double>>> Hexxs;
    Hexxs[0][ModuleIO::TAC{0, ModuleIO::TC{0, 0, 0}}] = matrix;

    const auto sparse = ModuleIO::calculate_RI_Tensor_sparse(1e-10, Hexxs, ucell);

    const Abfs::Vector3_Order<int> r_vector(0, 0, 0);
    ASSERT_EQ(sparse.count(r_vector), 1);
    const auto& block = sparse.at(r_vector);
    ASSERT_EQ(block.count(0), 1);
    EXPECT_EQ(block.at(0).count(0), 1);
    EXPECT_EQ(block.at(0).at(0), 1.0);
    EXPECT_EQ(block.at(0).count(1), 0);
    ASSERT_EQ(block.count(1), 1);
    EXPECT_EQ(block.at(1).count(0), 0);
    EXPECT_EQ(block.at(1).count(1), 1);
    EXPECT_EQ(block.at(1).at(1), -2.0);
}
