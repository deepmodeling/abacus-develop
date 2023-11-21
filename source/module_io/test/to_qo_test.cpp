#include <gtest/gtest.h>
#include "module_io/to_qo.h"

Atom_pseudo::Atom_pseudo() {}
Atom_pseudo::~Atom_pseudo() {}
#ifdef __MPI
void Atom_pseudo::bcast_atom_pseudo() {}
#endif
pseudo::pseudo() {}
pseudo::~pseudo() {}

Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}
#ifdef __LCAO
InfoNonlocal::InfoNonlocal() {}
InfoNonlocal::~InfoNonlocal() {}
#endif
void output::printM3(std::ofstream &ofs, const std::string &description, const ModuleBase::Matrix3 &m) {}

class toQOTest : public testing::Test
{
  protected:
    void SetUp() override
    {    
        ucell.atoms = new Atom[2];
        ucell.set_atom_flag = true;
        ucell.ntype = 2;
        ucell.atoms[0].tau = new ModuleBase::Vector3<double>[1];
        ucell.atoms[1].tau = new ModuleBase::Vector3<double>[1];
        ucell.atoms[0].tau[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
        ucell.atoms[1].tau[0] = ModuleBase::Vector3<double>(4.0, 4.0, 4.0);
        ucell.atoms[0].taud = new ModuleBase::Vector3<double>[1];
        ucell.atoms[1].taud = new ModuleBase::Vector3<double>[1];
        ucell.atoms[0].taud[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
        ucell.atoms[1].taud[0] = ModuleBase::Vector3<double>(0.5, 0.5, 0.5);
        ucell.atoms[0].na = 1;
        ucell.atoms[1].na = 1;
        ucell.a1 = ModuleBase::Vector3<double>(8.0, 8.0, 0.0);
        ucell.a2 = ModuleBase::Vector3<double>(8.0, 0.0, 8.0);
        ucell.a3 = ModuleBase::Vector3<double>(0.0, 8.0, 8.0);
        ucell.atoms[0].ncpp.zv = 4.0;
        ucell.atoms[1].ncpp.zv = 4.0;
        ucell.atoms[0].ncpp.psd = "Si";
        ucell.atoms[1].ncpp.psd = "C";
        ucell.latvec.e11 = 8.0; ucell.latvec.e12 = 8.0; ucell.latvec.e13 = 0.0;
        ucell.latvec.e21 = 8.0; ucell.latvec.e22 = 0.0; ucell.latvec.e23 = 8.0;
        ucell.latvec.e31 = 0.0; ucell.latvec.e32 = 8.0; ucell.latvec.e33 = 8.0;
        ucell.GT = ucell.latvec.Inverse();
        ucell.G = ucell.GT.Transpose();
        ucell.GGT = ucell.G * ucell.GT;
        ucell.orbital_fn = new std::string[2];
        ucell.orbital_fn[0] = "../../../../tests/PP_ORB/Si_gga_8au_100Ry_2s2p1d.orb";
        ucell.orbital_fn[1] = "../../../../tests/PP_ORB/C_gga_8au_100Ry_2s2p1d.orb";
    }

    void TearDown() override
    {
    }
    UnitCell ucell;
};

TEST_F(toQOTest, Constructor)
{
    toQO tqo("hydrogen");
    EXPECT_EQ(tqo.qo_basis(), "hydrogen");
    EXPECT_EQ(tqo.strategy(), "minimal");
    EXPECT_EQ(tqo.nkpts(), 0);
    EXPECT_EQ(tqo.p_ucell(), nullptr);
}

TEST_F(toQOTest, Initialize)
{
    toQO tqo("hydrogen");
    UnitCell ucell;
    tqo.initialize(&ucell, 1);
    EXPECT_EQ(tqo.nkpts(), 1);
    EXPECT_EQ(tqo.p_ucell(), &ucell);
}

TEST_F(toQOTest, EliminateDuplicateVector3)
{
    std::vector<ModuleBase::Vector3<int>> v;
    v.push_back(ModuleBase::Vector3<int>(0, 0, 0));
    v.push_back(ModuleBase::Vector3<int>(0, 0, 0));
    v.push_back(ModuleBase::Vector3<int>(0, 0, 0));
    v.push_back(ModuleBase::Vector3<int>(1, 0, 0));
    v.push_back(ModuleBase::Vector3<int>(1, 0, 0));
    v.push_back(ModuleBase::Vector3<int>(1, 0, 0));
    v.push_back(ModuleBase::Vector3<int>(1, 1, 0));
    v.push_back(ModuleBase::Vector3<int>(1, 1, 0));
    v.push_back(ModuleBase::Vector3<int>(1, 1, 0));
    v.push_back(ModuleBase::Vector3<int>(1, 1, 1));
    v.push_back(ModuleBase::Vector3<int>(1, 1, 1));
    v.push_back(ModuleBase::Vector3<int>(1, 1, 1));
    toQO tqo("hydrogen");
    tqo.eliminate_duplicate_vector3(v);
    EXPECT_EQ(v.size(), 4);
}

TEST_F(toQOTest, Norm2RijSupercell)
{
    toQO tqo("hydrogen");
    tqo.initialize(&ucell, 1);
    ModuleBase::Vector3<double> rij(1.0, 0.0, 0.0);
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 0, 0, 0), 2.0);
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 1, 0, 0), 146.0);
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 0, 1, 0), 146.0);
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 0, 0, 1), 130.0);
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 1, 1, 0), 418.0);
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 1, 0, 1), 402.0);
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 0, 1, 1), 402.0);
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 1, 1, 1), 802.0);
}

TEST_F(toQOTest, ScanSupercellForAtom)
{
    toQO tqo("hydrogen");
    tqo.initialize(&ucell, 1);
    std::vector<ModuleBase::Vector3<int>> n1n2n3 = tqo.scan_supercell_for_atom(0, 0);
    EXPECT_EQ(n1n2n3.size(), 25);
    for(int iR = 0; iR < n1n2n3.size(); iR++)
    {
        std::cout << n1n2n3[iR][0] << " " << n1n2n3[iR][1] << " " << n1n2n3[iR][2] << std::endl;
    }
}
int main(int argc, char** argv)
{
    // current getcwd()
    std::cout << "Current getcwd: " << getcwd(nullptr, 0) << std::endl;

#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif
    
    return result;
}
