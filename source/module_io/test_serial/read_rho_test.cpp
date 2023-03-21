#include "gtest/gtest.h"
#include "module_io/rho_io.h"
#include "module_base/global_variable.h"
#include "prepare_unitcell.h"

#ifdef __LCAO
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
LCAO_Orbitals::LCAO_Orbitals(){}
LCAO_Orbitals::~LCAO_Orbitals(){}
#endif
Magnetism::Magnetism(){}
Magnetism::~Magnetism(){}

/************************************************
 *  unit test of read_rho
 ***********************************************/

/**
 * - Tested Functions:
 *   - read_rho()
 *     - the function to read_rho from file
 *     - the serial version without MPI
 */

class ReadRhoTest : public ::testing::Test
{
protected:
	int nspin = 1;
	int nrxx = 36*36*36;
	int prenspin = 1;
	double **rho;
	std::unique_ptr<UnitCell> ucell{new UnitCell};
	void SetUp()
	{
		rho = new double*[nspin];
		for(int is=0; is<nspin; ++is)
		{
			rho[is] = new double[nrxx];
		}
	}
	void TearDown()
	{
		for(int is=0; is<nspin; ++is)
		{
			delete[] rho[is];
		}
		delete[] rho;
	}
};

TEST_F(ReadRhoTest,Read)
{
	int is = 0;
	std::string fn = "./support/SPIN1_CHG.cube";
	int nx = 36;
	int ny = 36;
	int nz = 36;
	double ef;
	UcellTestPrepare utp = UcellTestLib["Si"];
	ucell = utp.SetUcellInfo();
	ModuleIO::read_rho(is,nspin,fn,rho[is],nx,ny,nz,ef,*ucell,prenspin);
	EXPECT_DOUBLE_EQ(ef,0.461002);
	EXPECT_DOUBLE_EQ(rho[0][0],1.27020863940e-03);
	EXPECT_DOUBLE_EQ(rho[0][46655],1.33581335706e-02);
}
