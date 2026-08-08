#include "gtest/gtest.h"
#include "gmock/gmock.h"
#define private public
#include "source_io/module_parameter/parameter.h"
#undef private
#include "source_io/module_output/write_orb_info.h"
#include "source_cell/unitcell.h"
#include "prepare_unitcell.h"
#include "source_cell/read_pseudo.h"

Magnetism::Magnetism()
{
	this->tot_mag = 0.0;
	this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

/************************************************
 *  unit test of write_orb_info
 ***********************************************/

/**
 * - Tested Functions:
 *   - write_orb_info()
 */


TEST(OrbInfo,WriteOrbInfo)
{
    UnitCell* ucell = new UnitCell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    ucell = utp.SetUcellInfo();
    ucell->set_iat2itia();
    std::string pp_dir = "./support/";
    std::ofstream ofs;
    ofs.open("running.log");
    const std::string global_out_dir = "./";
    const std::string dft_functional = "default";
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    const int nspin = 1;
    const int nlocal = 18;
    const int npol = 1;
    const std::string basis_type = "pw";
    const std::string esolver_type = "ksdft";
    const std::string init_wfc = "";
    const int nbands = 6;
    unitcell::read_cell_pseudopots(pp_dir, ofs, *ucell, global_out_dir, dft_functional, lspinorb, pseudo_rcut, soc_lambda);
    unitcell::cal_nwfc(ofs,*ucell,ucell->atoms, nspin, nlocal, npol, basis_type, esolver_type, init_wfc, nbands);
    ModuleIO::write_orb_info(ucell);
    ofs.close();
    std::ifstream ifs("Orbital");
    std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("#io    spec    l    m    z  sym"));
    EXPECT_THAT(str, testing::HasSubstr("0      Si    0    0    1              s"));
    EXPECT_THAT(str, testing::HasSubstr("0      Si    2    4    1            dxy"));
    EXPECT_THAT(str, testing::HasSubstr("#sym =Symmetry name of real orbital"));
    ifs.close();
    delete ucell;
    std::remove("Orbital");
    std::remove("running.log");
}

#ifdef __MPI
#include "mpi.h"
int main(int argc, char **argv)
{


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}
#endif


