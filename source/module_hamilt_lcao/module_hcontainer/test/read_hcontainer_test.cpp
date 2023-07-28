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

class TEST_Parallel_Orbitals : public Parallel_Orbitals
{
    public:
    TEST_Parallel_Orbitals() : Parallel_Orbitals()
    {
    }
    ~TEST_Parallel_Orbitals()
    {
    }
    void set_serial(int& nrow_in, const int& ncol_in)
    {
        this->nrow = nrow_in;
        this->ncol = ncol_in;
        this->nloc = this->nrow * this->ncol;
    }
}

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
    Parallel_Orbitals ParaV;
}