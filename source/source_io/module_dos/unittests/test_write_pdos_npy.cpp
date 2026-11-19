#include "../write_pdos_npy.h"
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"

#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <stdexcept>

namespace
{
class WritePdosNpy : public testing::Test
{
protected:
    UnitCell cell;
    std::vector<ModuleBase::matrix> pdos{1};
    std::string prefix;

    void SetUp() override
    {
        prefix = testing::TempDir() + "/abacus-" + testing::UnitTest::GetInstance()->current_test_info()->name() + "-";
        cell.nat = 1;
        cell.ntype = 1;
        cell.atoms = new Atom[1];
        cell.set_atom_flag = true;
        cell.iat2it = {0};
        cell.atoms[0].label = "Si";
        cell.atoms[0].nw = 2;
        cell.atoms[0].nwl = 0;
        cell.atoms[0].iw2l = {0, 0};
        cell.atoms[0].iw2m = {0, 0};
        cell.set_iat2iwt_for_test({0}, 1);
        pdos.front().create(2, 2);
        pdos.front()(0, 0) = 4.25;
        pdos.front()(1, 0) = -1.25;
        pdos.front()(0, 1) = 1e-8;
        pdos.front()(1, 1) = 0;
    }

    std::string write(int step, int mode)
    {
        return ModuleIO::write_pdos_npy(cell, pdos, 1, -1.0, 0.25, prefix, step, mode);
    }

    std::string contents(const std::string& path)
    {
        std::ifstream in(path, std::ios::binary);
        return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }

    void TearDown() override
    {
        std::remove((prefix + "pdos.npy").c_str());
        std::remove((prefix + "pdosg2.npy").c_str());
    }
};

TEST_F(WritePdosNpy, FramesAppendOverwriteAndRecoverShortTail)
{
    const auto path = write(0, 1);
    const auto first = contents(path);
    ASSERT_EQ(first.substr(0, 8), std::string("\x93NUMPY\x01\x00", 8));
    const double expected[] = {3.0, 0.0};
    const auto tail = first.substr(first.size() - sizeof(expected));
    EXPECT_EQ(tail, std::string(reinterpret_cast<const char*>(expected), sizeof(expected)));
    const auto other = write(1, 2);
    ASSERT_EQ(other, prefix + "pdosg2.npy");
    const auto second = contents(other);
    {
        std::ofstream interrupted(path, std::ios::app | std::ios::binary);
        interrupted.write(second.data(), second.size() - 1);
    }
    write(1, 1);
    EXPECT_EQ(contents(path), first + second);
    write(0, 1);
    EXPECT_EQ(contents(path), first);
    write(1, 0);
    EXPECT_EQ(contents(path), second);
}

TEST_F(WritePdosNpy, IncompatibleAppendPreservesFile)
{
    const auto path = write(0, 0);
    {
        std::ofstream existing(path, std::ios::binary);
        existing << "unrelated file";
    }
    EXPECT_THROW(write(1, 1), std::runtime_error);
    EXPECT_EQ(contents(path), "unrelated file");
}
}
