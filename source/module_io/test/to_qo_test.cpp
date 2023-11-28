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
        ucell.lat0 = 1.889726124565062;
        ucell.atoms[0].tau = new ModuleBase::Vector3<double>[1];
        ucell.atoms[1].tau = new ModuleBase::Vector3<double>[1];
        ucell.atoms[0].tau[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
        ucell.atoms[1].tau[0] = ModuleBase::Vector3<double>(2.0, 2.0, 2.0);
        ucell.atoms[0].taud = new ModuleBase::Vector3<double>[1];
        ucell.atoms[1].taud = new ModuleBase::Vector3<double>[1];
        ucell.atoms[0].taud[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
        ucell.atoms[1].taud[0] = ModuleBase::Vector3<double>(0.25, 0.25, 0.25);
        ucell.atoms[0].na = 1;
        ucell.atoms[1].na = 1;
        ucell.atoms[0].nwl = 2;
        ucell.atoms[1].nwl = 2;
        ucell.a1 = ModuleBase::Vector3<double>(8.0, 8.0, 0.0);
        ucell.a2 = ModuleBase::Vector3<double>(8.0, 0.0, 8.0);
        ucell.a3 = ModuleBase::Vector3<double>(0.0, 8.0, 8.0);
        ucell.atoms[0].ncpp.zv = 4.0;
        ucell.atoms[1].ncpp.zv = 4.0;
        ucell.atoms[0].ncpp.psd = "Si";
        ucell.atoms[1].ncpp.psd = "C";
        ucell.atoms[0].label = "Si";
        ucell.atoms[1].label = "C";
        ucell.latvec.e11 = 8.0; ucell.latvec.e12 = 8.0; ucell.latvec.e13 = 0.0;
        ucell.latvec.e21 = 8.0; ucell.latvec.e22 = 0.0; ucell.latvec.e23 = 8.0;
        ucell.latvec.e31 = 0.0; ucell.latvec.e32 = 8.0; ucell.latvec.e33 = 8.0;
        ucell.GT = ucell.latvec.Inverse();
        ucell.G = ucell.GT.Transpose();
        ucell.GGT = ucell.G * ucell.GT;
        ucell.orbital_fn = new std::string[2];
        ucell.orbital_fn[0] = "../../../../tests/PP_ORB/Si_gga_8au_100Ry_2s2p1d.orb";
        ucell.orbital_fn[1] = "../../../../tests/PP_ORB/C_gga_8au_100Ry_2s2p1d.orb";

        GlobalV::global_out_dir = "./";
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

TEST_F(toQOTest, UnwrapUnitcell)
{
    toQO tqo("hydrogen");
    tqo.unwrap_unitcell(&ucell);
    EXPECT_EQ(tqo.ntype(), ucell.ntype);
    EXPECT_EQ(tqo.symbols().size(), ucell.ntype);
    EXPECT_EQ(tqo.charges().size(), ucell.ntype);
    EXPECT_EQ(tqo.symbols()[0], "Si");
    EXPECT_EQ(tqo.symbols()[1], "C");
    EXPECT_EQ(tqo.charges()[0], 4.0);
    EXPECT_EQ(tqo.charges()[1], 4.0);
}

TEST_F(toQOTest, BuildNao)
{
    toQO tqo("hydrogen");
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, ucell.orbital_fn);
    EXPECT_EQ(tqo.p_nao()->nchi(), 10); // not (l, m)-resoluted
    EXPECT_EQ(tqo.nphi(), 26); // (l, m)-resoluted
}

TEST_F(toQOTest, BuildAo)
{
    toQO tqo("hydrogen");
    tqo.unwrap_unitcell(&ucell);
    std::vector<int> nmax = std::vector<int>(ucell.ntype);
    for(int itype = 0; itype < ucell.ntype; itype++)
    {
        nmax[itype] = tqo.atom_database().principle_quantum_number[tqo.symbols()[itype]];
    }
    tqo.build_ao(ucell.ntype, tqo.charges().data(), nmax.data());
    EXPECT_EQ(tqo.p_ao()->nchi(), 5); // Si: 1s, 2p, 3d, C: 1s, 2p
    EXPECT_EQ(tqo.nchi(), 13); // Si: 1s, 2px, 2py, 2pz, 3dz2, 3dxz, 3dyz, 3dx2-y2, 3dxy, C: 1s, 2px, 2py, 2pz
}
// the scan_supercell_for_atom() calls
TEST_F(toQOTest, Norm2RijSupercell)
{
    toQO tqo("hydrogen");
    tqo.unwrap_unitcell(&ucell);
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
// the scan_supercell() calls
TEST_F(toQOTest, ScanSupercellForAtom)
{
    toQO tqo("hydrogen");
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, ucell.orbital_fn);
    std::vector<int> nmax = std::vector<int>(ucell.ntype);
    for(int itype = 0; itype < ucell.ntype; itype++)
    {
        nmax[itype] = tqo.atom_database().principle_quantum_number[tqo.symbols()[itype]];
    }
    tqo.build_ao(ucell.ntype, tqo.charges().data(), nmax.data());
    std::vector<ModuleBase::Vector3<int>> n1n2n3 = tqo.scan_supercell_for_atom(0, 0);
    EXPECT_EQ(n1n2n3.size(), 22);
}
// the scan_supercell() calls
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
    tqo.eliminate_duplicate_vector3<int>(v);
    EXPECT_EQ(v.size(), 4);
}

TEST_F(toQOTest, ScanSupercell)
{
    toQO tqo("hydrogen");
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, ucell.orbital_fn);
    std::vector<int> nmax = std::vector<int>(ucell.ntype);
    for(int itype = 0; itype < ucell.ntype; itype++)
    {
        nmax[itype] = tqo.atom_database().principle_quantum_number[tqo.symbols()[itype]];
    }
    tqo.build_ao(ucell.ntype, tqo.charges().data(), nmax.data());
    tqo.scan_supercell();
    EXPECT_EQ(tqo.nR(), 25);
}

TEST_F(toQOTest, AllocateOvlp)
{
    toQO tqo("hydrogen");
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, ucell.orbital_fn);
    std::vector<int> nmax = std::vector<int>(ucell.ntype);
    for(int itype = 0; itype < ucell.ntype; itype++)
    {
        nmax[itype] = tqo.atom_database().principle_quantum_number[tqo.symbols()[itype]];
    }
    tqo.build_ao(ucell.ntype, tqo.charges().data(), nmax.data());
    tqo.scan_supercell();
    tqo.allocate_ovlp(true);
    tqo.allocate_ovlp(false);
    EXPECT_EQ(tqo.ovlp_R().size(), 1); // in total 25 Rs, therefore 25 ovlp_R. but save_mem, so 1
    EXPECT_EQ(tqo.ovlp_k().size(), tqo.nchi()); // for single kpoint, ao*nao matrix
    EXPECT_EQ(tqo.ovlp_R()[0].size(), tqo.nchi()); // for single cell, ao*nao matrix
    EXPECT_EQ(tqo.ovlp_k()[0].size(), tqo.nphi()); // for each atomic orbital at single kpoint, get number of nao
    EXPECT_EQ(tqo.ovlp_R()[0][0].size(), tqo.nphi()); // similarly
    // all values in them should be zero or complex zero
    for(int iR = 0; iR < tqo.ovlp_R().size(); iR++)
    {
        for(int i = 0; i < tqo.ovlp_R()[iR].size(); i++)
        {
            for(int j = 0; j < tqo.ovlp_R()[iR][i].size(); j++)
            {
                EXPECT_EQ(tqo.ovlp_R()[iR][i][j], 0.0);
            }
        }
    }
    for(int i = 0; i < tqo.ovlp_k().size(); i++)
    {
        for(int j = 0; j < tqo.ovlp_k()[i].size(); j++)
        {
            EXPECT_EQ(tqo.ovlp_k()[i][j], std::complex<double>(0.0, 0.0));
        }
    }
}

TEST_F(toQOTest, Initialize)
{
    toQO tqo("hydrogen");
    std::vector<ModuleBase::Vector3<double>> kvecs_c;
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0)); // Gamma point
    tqo.initialize(&ucell, kvecs_c);
}

TEST_F(toQOTest, CalculateOvlpR)
{
    toQO tqo("hydrogen");
    std::vector<ModuleBase::Vector3<double>> kvecs_c;
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0)); // Gamma point
    tqo.initialize(&ucell, kvecs_c);
    tqo.calculate_ovlp_R(12);
    // not all elements are zero
    bool all_zero = true;
    for(int i = 0; i < tqo.ovlp_R()[0].size(); i++)
    {
        for(int j = 0; j < tqo.ovlp_R()[0][i].size(); j++)
        {
            if(tqo.ovlp_R()[0][i][j] != 0.0)
            {
                all_zero = false;
            }
        }
    }
    EXPECT_EQ(all_zero, false);
}

TEST_F(toQOTest, CalculateOvlpK)
{
    toQO tqo("hydrogen");
    std::vector<ModuleBase::Vector3<double>> kvecs_c;
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0)); // Gamma point
    tqo.initialize(&ucell, kvecs_c);
    tqo.calculate_ovlp_k(0);
    // all should be real numbers at Gamma point
    bool all_real = true;
    for(int i = 0; i < tqo.ovlp_k().size(); i++)
    {
        for(int j = 0; j < tqo.ovlp_k()[i].size(); j++)
        {
            if(tqo.ovlp_k()[i][j].imag() != 0.0)
            {
                all_real = false;
            }
        }
    }
}

TEST_F(toQOTest, Calculate)
{
    toQO tqo("hydrogen");
    std::vector<ModuleBase::Vector3<double>> kvecs_c;
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0)); // Gamma point
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.5, 0.0, 0.0));
    tqo.initialize(&ucell, kvecs_c);
    tqo.calculate();
    // for the latest kpoint, not all numbers are complex zero
    bool all_zero = true;
    for(int i = 0; i < tqo.ovlp_k().size(); i++)
    {
        for(int j = 0; j < tqo.ovlp_k()[i].size(); j++)
        {
            if(tqo.ovlp_k()[i][j] != std::complex<double>(0.0, 0.0))
            {
                all_zero = false;
            }
        }
    }
    EXPECT_EQ(all_zero, false);
    // delete files generated namely QO_ovlp_0.dat and QO_ovlp_1.dat
    std::remove("QO_ovlp_0.dat");
    std::remove("QO_ovlp_1.dat");
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
