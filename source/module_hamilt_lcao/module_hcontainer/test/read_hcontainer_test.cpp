#include <fstream>

#include "../hcontainer.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_io/csr_reader.h"
#include "prepare_unitcell.h"

// mock functions
#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
LCAO_Orbitals::LCAO_Orbitals()
{
}
LCAO_Orbitals::~LCAO_Orbitals()
{
}
#endif
Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}
// mocke functions

class ReadHContainerTest : public testing::Test
{
  protected:
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    std::vector<int> nw = {13};
    int nlocal;
    void SetUp() override
    {
        ucell = utp.SetUcellInfo(nw, nlocal);
    }
    void TearDown() override
    {
        delete ucell;
    }
};

TEST_F(ReadHContainerTest, Foo)
{
    std::string filename = "./support/SR.csr";
    ModuleIO::csrFileReader<double> csr(filename);
    std::cout << "csr.getStep " << csr.getStep() << std::endl;
    std::cout << "csr.getMatrixDimension " << csr.getMatrixDimension() << std::endl;
    std::cout << "csr.getNumberOfR " << csr.getNumberOfR() << std::endl;
    std::cout << "nlocal " << nlocal << std::endl;
    Parallel_Orbitals paraV;
    std::ofstream ofs("test.log");
    paraV.set_global2local(nlocal, nlocal, false, ofs);
    ofs.close();
    remove("test.log");
    paraV.set_atomic_trace(ucell->iat2iwt.data(), ucell->nat, nlocal);
    std::cout << paraV.atom_begin_col[0] << " " << paraV.atom_begin_col[1] << std::endl;
    std::cout << paraV.atom_begin_row[0] << " " << paraV.atom_begin_row[1] << std::endl;

}