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
        ucell.pseudo_fn = new std::string[2];
        ucell.pseudo_fn[0] = "../../../../tests/PP_ORB/Si_dojo_soc.upf";
        ucell.pseudo_fn[1] = "../../../../tests/PP_ORB/C.LDA.UPF";

        GlobalV::global_out_dir = "./";
        GlobalV::qo_screening_coeff = {0.1, 0.1};
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
    toQO tqo("hydrogen", "minimal");
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, ucell.orbital_fn);
    EXPECT_EQ(tqo.p_nao()->nchi(), 10); // not (l, m)-resoluted
    EXPECT_EQ(tqo.nphi(), 26); // (l, m)-resoluted
}

TEST_F(toQOTest, BuildHydrogenMinimal)
{
    toQO tqo("hydrogen", "minimal");
    tqo.unwrap_unitcell(&ucell);
    tqo.build_ao(ucell.ntype);
    EXPECT_EQ(tqo.p_ao()->nchi(), 5); // Si: 1s, 2p, 3d, C: 1s, 2p
    EXPECT_EQ(tqo.nchi(), 13); // Si: 1s, 2px, 2py, 2pz, 3dz2, 3dxz, 3dyz, 3dx2-y2, 3dxy, C: 1s, 2px, 2py, 2pz
    tqo.p_ao()->to_file("special_use_unittest");
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
    tqo.set_strategy("minimal");
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, ucell.orbital_fn);
    std::vector<int> nmax = std::vector<int>(ucell.ntype);
    for(int itype = 0; itype < ucell.ntype; itype++)
    {
        nmax[itype] = tqo.atom_database().principle_quantum_number[tqo.symbols()[itype]];
    }
    tqo.build_ao(ucell.ntype);
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
    tqo.set_strategy("minimal");
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, ucell.orbital_fn);
    std::vector<int> nmax = std::vector<int>(ucell.ntype);
    for(int itype = 0; itype < ucell.ntype; itype++)
    {
        nmax[itype] = tqo.atom_database().principle_quantum_number[tqo.symbols()[itype]];
    }
    tqo.build_ao(ucell.ntype);
    tqo.scan_supercell();
    EXPECT_EQ(tqo.nR(), 25);
}

TEST_F(toQOTest, AllocateOvlpMinimal)
{
    toQO tqo("hydrogen");
    tqo.set_strategy("minimal");
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, ucell.orbital_fn);
    std::vector<int> nmax = std::vector<int>(ucell.ntype);
    for(int itype = 0; itype < ucell.ntype; itype++)
    {
        nmax[itype] = tqo.atom_database().principle_quantum_number[tqo.symbols()[itype]];
    }
    tqo.build_ao(ucell.ntype);
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

TEST_F(toQOTest, CalculateSelfOvlpRMinimal)
{
    toQO tqo("hydrogen");
    tqo.set_strategy("minimal");
    std::vector<ModuleBase::Vector3<double>> kvecs_c;
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0)); // Gamma point
    ucell.orbital_fn[0] = "Si_special_use_unittest.orb"; // generated in unittest BuildAo
    ucell.orbital_fn[1] = "C_special_use_unittest.orb"; // generated in unittest BuildAo
    ucell.atoms[1].nwl = 1; // only s and p for C
    tqo.initialize(&ucell, kvecs_c);
    // find the R = 0,0,0
    for(int iR = 0; iR < tqo.nR(); iR++)
    {
        if(tqo.supercells()[iR].x == 0 && tqo.supercells()[iR].y == 0 && tqo.supercells()[iR].z == 0)
        {
            tqo.calculate_ovlp_R(iR);
            break;
        }
    }
    // check if diagonal elements are 1
    for(int i = 0; i < tqo.nphi(); i++)
    {
        EXPECT_NEAR(tqo.ovlp_R()[0][i][i], 1.0, 1e-4);
    }
    //std::remove("Si_special_use_unittest.orb");
    //std::remove("C_special_use_unittest.orb");
    //tqo.write_ovlp(tqo.ovlp_R()[0], "QO_self_ovlp.dat");
}

TEST_F(toQOTest, CalculateSelfOvlpKSymmetrical)
{
    toQO tqo("hydrogen");
    tqo.set_strategy("minimal");
    std::vector<ModuleBase::Vector3<double>> kvecs_c;
    kvecs_c.push_back(ModuleBase::Vector3<double>(-0.25, -0.25, -0.25)); // pair 1
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.25, 0.25, 0.25));
    kvecs_c.push_back(ModuleBase::Vector3<double>(-0.25, 0.25, 0.25)); // pair 2
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.25, -0.25, -0.25));
    kvecs_c.push_back(ModuleBase::Vector3<double>(-0.25, -0.25, 0.25)); // pair 3
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.25, 0.25, -0.25));
    kvecs_c.push_back(ModuleBase::Vector3<double>(-0.25, 0.25, -0.25)); // pair 4
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.25, -0.25, 0.25));
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0)); // Gamma
    ucell.orbital_fn[0] = "Si_special_use_unittest.orb"; // generated in unittest BuildAo
    ucell.orbital_fn[1] = "C_special_use_unittest.orb"; // generated in unittest BuildAo
    ucell.atoms[1].nwl = 1; // only s and p for C
    tqo.initialize(&ucell, kvecs_c);
    // test symmetry cancellation on pair1
    tqo.calculate_ovlp_k(0);
    std::vector<std::vector<std::complex<double>>> ovlp_k_1 = tqo.ovlp_k();
    tqo.calculate_ovlp_k(1);
    std::vector<std::vector<std::complex<double>>> ovlp_k_2 = tqo.ovlp_k();
    bool all_zero = true;
    // check if all imaginary parts are cancelled
    for(int i = 0; i < ovlp_k_1.size(); i++)
    {
        for(int j = 0; j < ovlp_k_1[i].size(); j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlp_R_ij = ovlp_k_1[i][j] + ovlp_k_2[i][j];
            EXPECT_NEAR(ovlp_R_ij.imag(), 0.0, 1e-10);
            if(ovlp_R_ij.real() > 1e-10)
            {
                all_zero = false;
            }
        }
    }
    EXPECT_FALSE(all_zero);
    // test symmetry cancellation on pair2
    tqo.calculate_ovlp_k(2);
    std::vector<std::vector<std::complex<double>>> ovlp_k_3 = tqo.ovlp_k();
    tqo.calculate_ovlp_k(3);
    std::vector<std::vector<std::complex<double>>> ovlp_k_4 = tqo.ovlp_k();
    all_zero = true;
    // check if all imaginary parts are cancelled
    for(int i = 0; i < ovlp_k_3.size(); i++)
    {
        for(int j = 0; j < ovlp_k_3[i].size(); j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlp_R_ij = ovlp_k_3[i][j] + ovlp_k_4[i][j];
            EXPECT_NEAR(ovlp_R_ij.imag(), 0.0, 1e-10);
            if(ovlp_R_ij.real() > 1e-10)
            {
                all_zero = false;
            }
        }
    }
    EXPECT_FALSE(all_zero);
    // test symmetry cancellation on pair3
    tqo.calculate_ovlp_k(4);
    std::vector<std::vector<std::complex<double>>> ovlp_k_5 = tqo.ovlp_k();
    tqo.calculate_ovlp_k(5);
    std::vector<std::vector<std::complex<double>>> ovlp_k_6 = tqo.ovlp_k();
    all_zero = true;
    // check if all imaginary parts are cancelled
    for(int i = 0; i < ovlp_k_5.size(); i++)
    {
        for(int j = 0; j < ovlp_k_5[i].size(); j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlp_R_ij = ovlp_k_5[i][j] + ovlp_k_6[i][j];
            EXPECT_NEAR(ovlp_R_ij.imag(), 0.0, 1e-10);
            if(ovlp_R_ij.real() > 1e-10)
            {
                all_zero = false;
            }
        }
    }
    EXPECT_FALSE(all_zero);
    // test symmetry cancellation on pair4
    tqo.calculate_ovlp_k(6);
    std::vector<std::vector<std::complex<double>>> ovlp_k_7 = tqo.ovlp_k();
    tqo.calculate_ovlp_k(7);
    std::vector<std::vector<std::complex<double>>> ovlp_k_8 = tqo.ovlp_k();
    all_zero = true;
    // check if all imaginary parts are cancelled
    for(int i = 0; i < ovlp_k_7.size(); i++)
    {
        for(int j = 0; j < ovlp_k_7[i].size(); j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlp_R_ij = ovlp_k_7[i][j] + ovlp_k_8[i][j];
            EXPECT_NEAR(ovlp_R_ij.imag(), 0.0, 1e-10);
            if(ovlp_R_ij.real() > 1e-10)
            {
                all_zero = false;
            }
        }
    }
    EXPECT_FALSE(all_zero);
    // test symmetry cancellation on pair5
    tqo.calculate_ovlp_k(8);
    std::vector<std::vector<std::complex<double>>> ovlp_k_9 = tqo.ovlp_k();
    all_zero = true;
    // check if all imaginary parts are cancelled
    for(int i = 0; i < ovlp_k_9.size(); i++)
    {
        for(int j = 0; j < ovlp_k_9[i].size(); j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            EXPECT_NEAR(ovlp_k_9[i][j].imag(), 0.0, 1e-10);
            if(ovlp_k_9[i][j].real() > 1e-10)
            {
                all_zero = false;
            }
        }
    }
    EXPECT_FALSE(all_zero);
    std::remove("Si_special_use_unittest.orb");
    std::remove("C_special_use_unittest.orb");
    //tqo.write_ovlp(tqo.ovlp_R()[0], "QO_self_ovlp.dat");
}

TEST_F(toQOTest, BuildHydrogenFull)
{
    toQO tqo("hydrogen", "full");
    tqo.unwrap_unitcell(&ucell);
    GlobalV::qo_thr = 1e-10;
    tqo.build_ao(ucell.ntype);
    EXPECT_EQ(tqo.p_ao()->nchi(), 9); // Si: 1s, 2s, 2p, 3s, 3p, 3d, C: 1s, 2s, 2p
    EXPECT_EQ(tqo.nchi(), 19); 
    tqo.p_ao()->to_file("special_use_unittest");
}

TEST_F(toQOTest, CalculateSelfOvlpRFull)
{
    toQO tqo("hydrogen", "full");
    std::vector<ModuleBase::Vector3<double>> kvecs_c;
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0)); // Gamma point
    ucell.orbital_fn[0] = "Si_special_use_unittest.orb"; // generated in unittest BuildAo
    ucell.orbital_fn[1] = "C_special_use_unittest.orb"; // generated in unittest BuildAo
    ucell.atoms[1].nwl = 1; // only s and p for C
    GlobalV::qo_thr = 1e-10;
    tqo.initialize(&ucell, kvecs_c);
    // find the R = 0,0,0
    for(int iR = 0; iR < tqo.nR(); iR++)
    {
        if(tqo.supercells()[iR].x == 0 && tqo.supercells()[iR].y == 0 && tqo.supercells()[iR].z == 0)
        {
            tqo.calculate_ovlp_R(iR);
            break;
        }
    }
    // check if diagonal elements are 1
    for(int i = 0; i < tqo.nphi(); i++)
    {
        EXPECT_NEAR(tqo.ovlp_R()[0][i][i], 1.0, 1e-4);
    }
    // check if symmetrical
    for(int i = 0; i < tqo.nchi(); i++)
    {
        for(int j = 0; j < tqo.nphi(); j++)
        {
            EXPECT_NEAR(tqo.ovlp_R()[0][i][j], tqo.ovlp_R()[0][j][i], 1e-4);
        }
    }
    std::remove("Si_special_use_unittest.orb");
    std::remove("C_special_use_unittest.orb");
    //tqo.write_ovlp(tqo.ovlp_R()[0], "QO_self_ovlp.dat");
}

TEST_F(toQOTest, BuildPswfc)
{
    toQO tqo("pswfc");
    tqo.unwrap_unitcell(&ucell);
    tqo.build_ao(ucell.ntype, ucell.pseudo_fn);
    EXPECT_EQ(tqo.p_ao()->nchi(), 5); 
    EXPECT_EQ(tqo.nchi(), 11);
    tqo.p_ao()->to_file("special_use_unittest");
}

TEST_F(toQOTest, CalculateSelfOvlpRPswfc)
{
    toQO tqo("pswfc");
    std::vector<ModuleBase::Vector3<double>> kvecs_c;
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0)); // Gamma point
    ucell.orbital_fn[0] = "Si_special_use_unittest.orb"; // generated in unittest BuildAo
    ucell.orbital_fn[1] = "C_special_use_unittest.orb"; // generated in unittest BuildAo
    ucell.atoms[1].nwl = 1; // only s and p for C
    //GlobalV::qo_thr = 1e-10;
    tqo.initialize(&ucell, kvecs_c);
    // find the R = 0,0,0
    for(int iR = 0; iR < tqo.nR(); iR++)
    {
        if(tqo.supercells()[iR].x == 0 && tqo.supercells()[iR].y == 0 && tqo.supercells()[iR].z == 0)
        {
            tqo.calculate_ovlp_R(iR);
            break;
        }
    }
    // check if diagonal elements are 1
    for(int i = 0; i < tqo.nphi(); i++)
    {
        EXPECT_NEAR(tqo.ovlp_R()[0][i][i], 1.0, 1e-4);
    }
    // check if symmetrical
    for(int i = 0; i < tqo.nchi(); i++)
    {
        for(int j = 0; j < tqo.nphi(); j++)
        {
            EXPECT_NEAR(tqo.ovlp_R()[0][i][j], tqo.ovlp_R()[0][j][i], 1e-4);
        }
    }
    std::remove("Si_special_use_unittest.orb");
    std::remove("C_special_use_unittest.orb");
    //tqo.write_ovlp(tqo.ovlp_R()[0], "QO_self_ovlp.dat");
}

TEST_F(toQOTest, CalculateOvlpKGamma)
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

TEST_F(toQOTest, CalculateSelfOvlpKPswfcSymmetrical)
{
    toQO tqo("pswfc");
    std::vector<ModuleBase::Vector3<double>> kvecs_c;
    kvecs_c.push_back(ModuleBase::Vector3<double>(-0.25, -0.25, -0.25)); // pair 1
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.25, 0.25, 0.25));
    kvecs_c.push_back(ModuleBase::Vector3<double>(-0.25, 0.25, 0.25)); // pair 2
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.25, -0.25, -0.25));
    kvecs_c.push_back(ModuleBase::Vector3<double>(-0.25, -0.25, 0.25)); // pair 3
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.25, 0.25, -0.25));
    kvecs_c.push_back(ModuleBase::Vector3<double>(-0.25, 0.25, -0.25)); // pair 4
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.25, -0.25, 0.25));
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0)); // Gamma
    tqo.initialize(&ucell, kvecs_c);
    std::cout << "Number of supercells: " << tqo.nR() << ", number of kpoints: " << tqo.nkpts() << std::endl;
    tqo.calculate_ovlp_k(0);
    std::vector<std::vector<std::complex<double>>> ovlp_k_1 = tqo.ovlp_k();
    tqo.calculate_ovlp_k(1);
    std::vector<std::vector<std::complex<double>>> ovlp_k_2 = tqo.ovlp_k();
    // check if all imaginary parts are cancelled
    bool all_real = true;
    for(int i = 0; i < ovlp_k_1.size(); i++)
    {
        for(int j = 0; j < ovlp_k_1[i].size(); j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlp_R_ij = ovlp_k_1[i][j] + ovlp_k_2[i][j];
            EXPECT_NEAR(ovlp_R_ij.imag(), 0.0, 1e-10);
        }
    }
    tqo.calculate_ovlp_k(2);
    std::vector<std::vector<std::complex<double>>> ovlp_k_3 = tqo.ovlp_k();
    tqo.calculate_ovlp_k(3);
    std::vector<std::vector<std::complex<double>>> ovlp_k_4 = tqo.ovlp_k();
    // check if all imaginary parts are cancelled
    for(int i = 0; i < ovlp_k_3.size(); i++)
    {
        for(int j = 0; j < ovlp_k_3[i].size(); j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlp_R_ij = ovlp_k_3[i][j] + ovlp_k_4[i][j];
            EXPECT_NEAR(ovlp_R_ij.imag(), 0.0, 1e-10);
        }
    }
    tqo.calculate_ovlp_k(4);
    std::vector<std::vector<std::complex<double>>> ovlp_k_5 = tqo.ovlp_k();
    tqo.calculate_ovlp_k(5);
    std::vector<std::vector<std::complex<double>>> ovlp_k_6 = tqo.ovlp_k();
    // check if all imaginary parts are cancelled
    for(int i = 0; i < ovlp_k_5.size(); i++)
    {
        for(int j = 0; j < ovlp_k_5[i].size(); j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlp_R_ij = ovlp_k_5[i][j] + ovlp_k_6[i][j];
            EXPECT_NEAR(ovlp_R_ij.imag(), 0.0, 1e-10);
        }
    }
    tqo.calculate_ovlp_k(6);
    std::vector<std::vector<std::complex<double>>> ovlp_k_7 = tqo.ovlp_k();
    tqo.calculate_ovlp_k(7);
    std::vector<std::vector<std::complex<double>>> ovlp_k_8 = tqo.ovlp_k();
    // check if all imaginary parts are cancelled
    for(int i = 0; i < ovlp_k_7.size(); i++)
    {
        for(int j = 0; j < ovlp_k_7[i].size(); j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlp_R_ij = ovlp_k_7[i][j] + ovlp_k_8[i][j];
            EXPECT_NEAR(ovlp_R_ij.imag(), 0.0, 1e-10);
        }
    }
    tqo.calculate_ovlp_k(8);
    std::vector<std::vector<std::complex<double>>> ovlp_k_9 = tqo.ovlp_k();
    // check if all imaginary parts are cancelled
    for(int i = 0; i < ovlp_k_9.size(); i++)
    {
        for(int j = 0; j < ovlp_k_9[i].size(); j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            EXPECT_NEAR(ovlp_k_9[i][j].imag(), 0.0, 1e-10);
        }
    }
    //tqo.write_ovlp(tqo.ovlp_R()[0], "QO_self_ovlp.dat");
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
