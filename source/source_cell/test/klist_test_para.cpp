#include "source_base/mathzone.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_global.h"
#define private public
#include "source_io/module_parameter/parameter.h"
#undef private
#include "source_cell/parallel_kpoints.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <iostream>
#include <streambuf>
#define private public
#include "../klist.h"
#include "source_basis/module_ao/ORB_gaunt_table.h"
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/parallel_kpoints.h"
#include "source_cell/pseudo.h"
#include "source_cell/setup_nonlocal.h"
#include "source_cell/unitcell.h"
#include "source_estate/magnetism.h"
#include "source_pw/module_pwdft/vl_pw.h"
#include "source_pw/module_pwdft/vnl_pw.h"
#include "source_pw/module_pwdft/parallel_grid.h"
#include "source_io/module_unk/berryphase.h"
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
SepPot::SepPot(){}
SepPot::~SepPot(){}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}


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
    EXPECT_EQ(kv->get_nkstot_full(), 512);
    EXPECT_GT(kv->get_nkstot_full(), kv->get_nkstot());
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
    std::vector<double> local_kvec_c_full(kv->kvec_c_full.size() * 3);
    for (size_t ik = 0; ik < kv->kvec_c_full.size(); ++ik)
    {
        local_kvec_c_full[3 * ik] = kv->kvec_c_full[ik].x;
        local_kvec_c_full[3 * ik + 1] = kv->kvec_c_full[ik].y;
        local_kvec_c_full[3 * ik + 2] = kv->kvec_c_full[ik].z;
    }
    const int local_count = static_cast<int>(local_kvec_c_full.size());
    std::vector<int> counts;
    std::vector<int> displs;
    std::vector<int> pools;
    if (GlobalV::MY_RANK == 0)
    {
        counts.resize(GlobalV::NPROC);
        pools.resize(GlobalV::NPROC);
    }
    MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&GlobalV::MY_POOL, 1, MPI_INT, pools.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<double> gathered_kvec_c_full;
    if (GlobalV::MY_RANK == 0)
    {
        displs.resize(GlobalV::NPROC, 0);
        for (int irank = 1; irank < GlobalV::NPROC; ++irank)
        {
            displs[irank] = displs[irank - 1] + counts[irank - 1];
        }
        gathered_kvec_c_full.resize(displs.back() + counts.back());
    }
    MPI_Gatherv(local_kvec_c_full.data(),
                local_count,
                MPI_DOUBLE,
                gathered_kvec_c_full.data(),
                counts.data(),
                displs.data(),
                MPI_DOUBLE,
                0,
                MPI_COMM_WORLD);
    if (GlobalV::MY_RANK == 0)
    {
        for (int irank = 0; irank < GlobalV::NPROC; ++irank)
        {
            for (int jrank = irank + 1; jrank < GlobalV::NPROC; ++jrank)
            {
                if (pools[irank] != pools[jrank])
                {
                    continue;
                }
                ASSERT_EQ(counts[irank], counts[jrank]);
                for (int i = 0; i < counts[irank]; ++i)
                {
                    EXPECT_NEAR(gathered_kvec_c_full[displs[irank] + i],
                                gathered_kvec_c_full[displs[jrank] + i],
                                1e-12);
                }
            }
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
//    kv->set_after_vc(kv->nspin, ucell.G, ucell.latvec);
    KVectorUtils::set_after_vc(*kv, kv->nspin, ucell.G);
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
