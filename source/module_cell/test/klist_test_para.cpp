#include "module_base/mathzone.h"
#include "module_base/parallel_global.h"
#define private public
#include "module_parameter/parameter.h"
#undef private
#include "module_cell/parallel_kpoints.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <iostream>
#include <streambuf>
#define private public
#include "../klist.h"
#include "module_basis/module_ao/ORB_gaunt_table.h"
#include "module_cell/atom_pseudo.h"
#include "module_cell/atom_spec.h"
#include "module_cell/parallel_kpoints.h"
#include "module_cell/pseudo.h"
#include "module_cell/setup_nonlocal.h"
#include "module_cell/unitcell.h"
#include "module_elecstate/magnetism.h"
#include "module_hamilt_pw/hamilt_pwdft/VL_in_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"
#include "module_io/berryphase.h"
#undef private
bool berryphase::berry_phase_flag = false;

pseudo::pseudo()
{
}
pseudo::~pseudo()
{
}
Atom::Atom()
{
}
Atom::~Atom()
{
}
Atom_pseudo::Atom_pseudo()
{
}
Atom_pseudo::~Atom_pseudo()
{
}
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
ORB_gaunt_table::ORB_gaunt_table()
{
}
ORB_gaunt_table::~ORB_gaunt_table()
{
}
pseudopot_cell_vl::pseudopot_cell_vl()
{
}
pseudopot_cell_vl::~pseudopot_cell_vl()
{
}
pseudopot_cell_vnl::pseudopot_cell_vnl()
{
}
pseudopot_cell_vnl::~pseudopot_cell_vnl()
{
}
Soc::~Soc()
{
}
Fcoef::~Fcoef()
{
}

namespace GlobalC
{
Parallel_Kpoints Pkpoints;
UnitCell ucell;
} // namespace GlobalC

/************************************************
 *  unit test of class K_Vectors
 ***********************************************/

/**
 * - Tested Functions:
 *   - Set
 *     - this is a "kind of" integerated test
 *       for set() and mpi_k()
 *   - SetAfterVC
 *     - this is a "kind of" integerated test
 *       for set_after_vc() and mpi_k_after_vc()
 *     - a bug is found from here, that is,
 *       KPAR > 1 is not support yet in vc-relax calculation
 *       due to the size of kvec_d, kvec_c being nks, rather
 *       than nkstot in set_both_kvec_after_vc
 */

// abbriviated from module_symmetry/test/symmetry_test.cpp
struct atomtype_
{
    std::string atomname;
    std::vector<std::vector<double>> coordinate;
};

struct stru_
{
    int ibrav;
    std::string point_group;    // Schoenflies symbol
    std::string point_group_hm; // Hermann-Mauguin notation.
    std::string space_group;
    std::vector<double> cell;
    std::vector<atomtype_> all_type;
};

std::vector<stru_> stru_lib{stru_{1,
                                  "O_h",
                                  "m-3m",
                                  "Pm-3m",
                                  std::vector<double>{1., 0., 0., 0., 1., 0., 0., 0., 1.},
                                  std::vector<atomtype_>{atomtype_{"C",
                                                                   std::vector<std::vector<double>>{
                                                                       {0., 0., 0.},
                                                                   }}}}};
// used to construct cell and analyse its symmetry

class KlistParaTest : public testing::Test
{
  protected:
    std::unique_ptr<K_Vectors> kv{new K_Vectors};
    std::ifstream ifs;
    std::ofstream ofs;
    std::ofstream ofs_running;
    std::string output;
    UnitCell ucell;
    // used to construct cell and analyse its symmetry
    void construct_ucell(stru_& stru)
    {
        std::vector<atomtype_> coord = stru.all_type;
        ucell.a1 = ModuleBase::Vector3<double>(stru.cell[0], stru.cell[1], stru.cell[2]);
        ucell.a2 = ModuleBase::Vector3<double>(stru.cell[3], stru.cell[4], stru.cell[5]);
        ucell.a3 = ModuleBase::Vector3<double>(stru.cell[6], stru.cell[7], stru.cell[8]);
        ucell.ntype = stru.all_type.size();
        ucell.atoms = new Atom[ucell.ntype];
        ucell.nat = 0;
        ucell.latvec.e11 = ucell.a1.x;
        ucell.latvec.e12 = ucell.a1.y;
        ucell.latvec.e13 = ucell.a1.z;
        ucell.latvec.e21 = ucell.a2.x;
        ucell.latvec.e22 = ucell.a2.y;
        ucell.latvec.e23 = ucell.a2.z;
        ucell.latvec.e31 = ucell.a3.x;
        ucell.latvec.e32 = ucell.a3.y;
        ucell.latvec.e33 = ucell.a3.z;
        ucell.GT = ucell.latvec.Inverse();
        ucell.G = ucell.GT.Transpose();
        ucell.lat0 = 1.8897261254578281;
        for (int i = 0; i < coord.size(); i++)
        {
            ucell.atoms[i].label = coord[i].atomname;
            ucell.atoms[i].na = coord[i].coordinate.size();
            ucell.atoms[i].tau.resize(ucell.atoms[i].na);
            ucell.atoms[i].taud.resize(ucell.atoms[i].na);
            for (int j = 0; j < ucell.atoms[i].na; j++)
            {
                std::vector<double> this_atom = coord[i].coordinate[j];
                ucell.atoms[i].tau[j] = ModuleBase::Vector3<double>(this_atom[0], this_atom[1], this_atom[2]);
                ModuleBase::Mathzone::Cartesian_to_Direct(ucell.atoms[i].tau[j].x,
                                                          ucell.atoms[i].tau[j].y,
                                                          ucell.atoms[i].tau[j].z,
                                                          ucell.a1.x,
                                                          ucell.a1.y,
                                                          ucell.a1.z,
                                                          ucell.a2.x,
                                                          ucell.a2.y,
                                                          ucell.a2.z,
                                                          ucell.a3.x,
                                                          ucell.a3.y,
                                                          ucell.a3.z,
                                                          ucell.atoms[i].taud[j].x,
                                                          ucell.atoms[i].taud[j].y,
                                                          ucell.atoms[i].taud[j].z);
            }
            ucell.nat += ucell.atoms[i].na;
        }
    }
    // clear ucell
    void ClearUcell()
    {
        delete[] ucell.atoms;
    }
};

#ifdef __MPI
TEST_F(KlistParaTest, Set)
{
    // construct cell and symmetry
    ModuleSymmetry::Symmetry symm;
    construct_ucell(stru_lib[0]);
    if (GlobalV::MY_RANK == 0) {
        GlobalV::ofs_running.open("tmp_klist_5");
}
    symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running);
    // read KPT
    std::string k_file = "./support/KPT1";
    // set klist
    kv->nspin = 1;
    PARAM.input.nspin = 1;
    if (GlobalV::NPROC == 4)
    {
        GlobalV::KPAR = 2;
    }
    Parallel_Global::init_pools(GlobalV::NPROC,
                                GlobalV::MY_RANK,
                                PARAM.input.bndpar,
                                GlobalV::KPAR,
                                GlobalV::NPROC_IN_BNDGROUP,
                                GlobalV::RANK_IN_BPGROUP,
                                GlobalV::MY_BNDGROUP,
                                GlobalV::NPROC_IN_POOL,
                                GlobalV::RANK_IN_POOL,
                                GlobalV::MY_POOL);
    ModuleSymmetry::Symmetry::symm_flag = 1;
    kv->set(ucell,symm, k_file, kv->nspin, ucell.G, ucell.latvec,  GlobalV::ofs_running);
    EXPECT_EQ(kv->get_nkstot(), 35);
    EXPECT_TRUE(kv->kc_done);
    EXPECT_TRUE(kv->kd_done);
    if (GlobalV::NPROC == 4)
    {
        if (GlobalV::MY_RANK == 0) {
            EXPECT_EQ(kv->get_nks(), 18);
}
        if (GlobalV::MY_RANK == 1) {
            EXPECT_EQ(kv->get_nks(), 18);
}
        if (GlobalV::MY_RANK == 2) {
            EXPECT_EQ(kv->get_nks(), 17);
}
        if (GlobalV::MY_RANK == 3) {
            EXPECT_EQ(kv->get_nks(), 17);
}
    }
    ClearUcell();
    if (GlobalV::MY_RANK == 0)
    {
        GlobalV::ofs_running.close();
        remove("tmp_klist_5");
        remove("kpoints");
    }
}

TEST_F(KlistParaTest, SetAfterVC)
{
    // construct cell and symmetry
    ModuleSymmetry::Symmetry symm;
    construct_ucell(stru_lib[0]);
    if (GlobalV::MY_RANK == 0) {
        GlobalV::ofs_running.open("tmp_klist_6");
}
    symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running);
    // read KPT
    std::string k_file = "./support/KPT1";
    // set klist
    kv->nspin = 1;
    PARAM.input.nspin = 1;
    if (GlobalV::NPROC == 4)
    {
        GlobalV::KPAR = 1;
    }
    Parallel_Global::init_pools(GlobalV::NPROC,
                                GlobalV::MY_RANK,
                                PARAM.input.bndpar,
                                GlobalV::KPAR,
                                GlobalV::NPROC_IN_BNDGROUP,
                                GlobalV::RANK_IN_BPGROUP,
                                GlobalV::MY_BNDGROUP,
                                GlobalV::NPROC_IN_POOL,
                                GlobalV::RANK_IN_POOL,
                                GlobalV::MY_POOL);
    ModuleSymmetry::Symmetry::symm_flag = 1;
    kv->set(ucell,symm, k_file, kv->nspin, ucell.G, ucell.latvec, GlobalV::ofs_running);
    EXPECT_EQ(kv->get_nkstot(), 35);
    EXPECT_TRUE(kv->kc_done);
    EXPECT_TRUE(kv->kd_done);
    if (GlobalV::NPROC == 4)
    {
        if (GlobalV::MY_RANK == 0) {
            EXPECT_EQ(kv->get_nks(), 35);
}
        if (GlobalV::MY_RANK == 1) {
            EXPECT_EQ(kv->get_nks(), 35);
}
        if (GlobalV::MY_RANK == 2) {
            EXPECT_EQ(kv->get_nks(), 35);
}
        if (GlobalV::MY_RANK == 3) {
            EXPECT_EQ(kv->get_nks(), 35);
}
    }
    // call set_after_vc here
    kv->kc_done = false;
    kv->set_after_vc(kv->nspin, ucell.G, ucell.latvec);
    EXPECT_TRUE(kv->kc_done);
    EXPECT_TRUE(kv->kd_done);
    // clear
    ClearUcell();
    if (GlobalV::MY_RANK == 0)
    {
        GlobalV::ofs_running.close();
        remove("tmp_klist_6");
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
    int result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}
#endif
