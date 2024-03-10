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

void define_fcc_cell(UnitCell& ucell)
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
    GlobalV::qo_thr = 1e-6;
    GlobalV::ofs_running = std::ofstream("unittest.log");
}

void define_sc_cell(UnitCell& ucell)
{    
    ucell.atoms = new Atom[1];
    ucell.set_atom_flag = true;
    ucell.ntype = 1;
    ucell.lat0 = 1.889726124565062;
    ucell.atoms[0].tau = new ModuleBase::Vector3<double>[1];
    ucell.atoms[0].tau[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    ucell.atoms[0].taud = new ModuleBase::Vector3<double>[1];
    ucell.atoms[0].taud[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    ucell.atoms[0].na = 1;
    ucell.atoms[0].nwl = 2;
    ucell.a1 = ModuleBase::Vector3<double>(8.0, 0.0, 0.0);
    ucell.a2 = ModuleBase::Vector3<double>(0.0, 8.0, 0.0);
    ucell.a3 = ModuleBase::Vector3<double>(0.0, 0.0, 8.0);
    ucell.atoms[0].ncpp.zv = 4.0;
    ucell.atoms[0].ncpp.psd = "Si";
    ucell.atoms[0].label = "Si";
    ucell.latvec.e11 = 8.0; ucell.latvec.e12 = 0.0; ucell.latvec.e13 = 0.0;
    ucell.latvec.e21 = 0.0; ucell.latvec.e22 = 8.0; ucell.latvec.e23 = 0.0;
    ucell.latvec.e31 = 0.0; ucell.latvec.e32 = 0.0; ucell.latvec.e33 = 8.0;
    ucell.GT = ucell.latvec.Inverse();
    ucell.G = ucell.GT.Transpose();
    ucell.GGT = ucell.G * ucell.GT;
    ucell.orbital_fn = new std::string[1];
    ucell.orbital_fn[0] = "../../../../tests/PP_ORB/Si_gga_8au_100Ry_2s2p1d.orb";
    ucell.pseudo_fn = new std::string[1];
    ucell.pseudo_fn[0] = "../../../../tests/PP_ORB/Si_dojo_soc.upf";

    GlobalV::global_out_dir = "./";
    GlobalV::qo_screening_coeff = {0.1};
    GlobalV::qo_thr = 1e-6;
    GlobalV::ofs_running = std::ofstream("unittest.log");
}

class toQOTest : public testing::Test
{
  protected:
    void SetUp() override
    {
    }

    void TearDown() override
    {
    }
    UnitCell ucell;
};

TEST_F(toQOTest, Constructor)
{
    define_fcc_cell(ucell);
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
    EXPECT_EQ(tqo.qo_basis(), "hydrogen");
    EXPECT_EQ(tqo.strategy(0), "minimal-nodeless");
    EXPECT_EQ(tqo.strategy(1), "minimal-nodeless");
    EXPECT_EQ(tqo.nkpts(), 0);
    EXPECT_EQ(tqo.p_ucell(), nullptr);
}

TEST_F(toQOTest, OrbitalFilter)
{
    define_fcc_cell(ucell);
    toQO tqo("pswfc", {"s", "spd", "all", "dfps"});
    EXPECT_TRUE(tqo.orbital_filter(0, tqo.strategy(0))); // whether l=0 is possible for stratgy of type 0
    EXPECT_FALSE(tqo.orbital_filter(1, tqo.strategy(0))); // whether l=1 is possible for stratgy of type 0
    EXPECT_FALSE(tqo.orbital_filter(2, tqo.strategy(0))); // whether l=2 is possible for stratgy of type 0
    EXPECT_TRUE(tqo.orbital_filter(0, tqo.strategy(1))); // whether l=0 is possible for stratgy of type 1
    EXPECT_TRUE(tqo.orbital_filter(1, tqo.strategy(1))); // whether l=1 is possible for stratgy of type 1
    EXPECT_TRUE(tqo.orbital_filter(2, tqo.strategy(1))); // whether l=2 is possible for stratgy of type 1
    EXPECT_FALSE(tqo.orbital_filter(3, tqo.strategy(1))); // whether l=3 is possible for stratgy of type 1
    EXPECT_TRUE(tqo.orbital_filter(0, tqo.strategy(2))); // whether l=0 is possible for stratgy of type 2
    EXPECT_TRUE(tqo.orbital_filter(1, tqo.strategy(2))); // whether l=1 is possible for stratgy of type 2
    EXPECT_TRUE(tqo.orbital_filter(2, tqo.strategy(2))); // whether l=2 is possible for stratgy of type 2
    EXPECT_TRUE(tqo.orbital_filter(3, tqo.strategy(2))); // whether l=3 is possible for stratgy of type 2
    EXPECT_TRUE(tqo.orbital_filter(4, tqo.strategy(2))); // whether l=4 is possible for stratgy of type 2
    EXPECT_TRUE(tqo.orbital_filter(5, tqo.strategy(2))); // whether l=5 is possible for stratgy of type 2
    EXPECT_TRUE(tqo.orbital_filter(6, tqo.strategy(2))); // whether l=6 is possible for stratgy of type 2
    EXPECT_TRUE(tqo.orbital_filter(0, tqo.strategy(3))); // whether l=0 is possible for stratgy of type 3
    EXPECT_TRUE(tqo.orbital_filter(1, tqo.strategy(3))); // whether l=1 is possible for stratgy of type 3
    EXPECT_TRUE(tqo.orbital_filter(2, tqo.strategy(3))); // whether l=2 is possible for stratgy of type 3
    EXPECT_TRUE(tqo.orbital_filter(3, tqo.strategy(3))); // whether l=3 is possible for stratgy of type 3
    EXPECT_FALSE(tqo.orbital_filter(4, tqo.strategy(3))); // whether l=4 is possible for stratgy of type 3
}

TEST_F(toQOTest, UnwrapUnitcell)
{
    define_fcc_cell(ucell);
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
    tqo.unwrap_unitcell(&ucell);
    EXPECT_EQ(tqo.ntype(), ucell.ntype);
    EXPECT_EQ(tqo.symbols().size(), ucell.ntype);
    EXPECT_EQ(tqo.charges().size(), ucell.ntype);
    EXPECT_EQ(tqo.symbols()[0], "Si");
    EXPECT_EQ(tqo.symbols()[1], "C");
    EXPECT_EQ(tqo.charges()[0], 14.0);
    EXPECT_EQ(tqo.charges()[1], 6.0);
}

TEST_F(toQOTest, BuildNao)
{
    define_fcc_cell(ucell);
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, "./", ucell.orbital_fn, 0);
    EXPECT_EQ(tqo.p_nao()->nchi(), 10); // not (l, m)-resoluted
    EXPECT_EQ(tqo.nphi(), 26); // (l, m)-resoluted
}

TEST_F(toQOTest, RadialCollectionIndexing)
{
    define_fcc_cell(ucell);
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, "./", ucell.orbital_fn, 0);
    // ucell.orbital_fn[0] = "../../../../tests/PP_ORB/Si_gga_8au_100Ry_2s2p1d.orb";
    // ucell.orbital_fn[1] = "../../../../tests/PP_ORB/C_gga_8au_100Ry_2s2p1d.orb";
    // test1 1Si, 1C
    std::vector<int> natoms = {1, 1};
    std::map<std::tuple<int,int,int,int,int>,int> index;
    std::map<int,std::tuple<int,int,int,int,int>> index_reverse;
    tqo.radialcollection_indexing(*(tqo.p_nao()), natoms, index, index_reverse);
    EXPECT_EQ(index.size(), 26); // 2*1 + 2*3 + 1*5 + 2*1 + 2*3 + 1*5 = 26
    EXPECT_EQ(index_reverse.size(), 26);
    // it, ia, l, izeta, m
    // it = 0
    // ia = 0
    // l = 0
    // izeta = 0
    EXPECT_EQ(index[std::make_tuple(0,0,0,0,0)], 0); // Si, 1st atom, 1s
    // izeta = 1
    EXPECT_EQ(index[std::make_tuple(0,0,0,1,0)], 1); // Si, 1st atom, 2s
    // l = 1
    // izeta = 0
    EXPECT_EQ(index[std::make_tuple(0,0,1,0,0)], 2); // Si, 1st atom, 1p, m = 0
    EXPECT_EQ(index[std::make_tuple(0,0,1,0,1)], 3); // Si, 1st atom, 1p, m = 1
    EXPECT_EQ(index[std::make_tuple(0,0,1,0,-1)], 4); // Si, 1st atom, 1p, m = -1
    // izeta = 1
    EXPECT_EQ(index[std::make_tuple(0,0,1,1,0)], 5); // Si, 1st atom, 2p, m = 0
    EXPECT_EQ(index[std::make_tuple(0,0,1,1,1)], 6); // Si, 1st atom, 2p, m = 1
    EXPECT_EQ(index[std::make_tuple(0,0,1,1,-1)], 7); // Si, 1st atom, 2p, m = -1
    // l = 2
    // izeta = 0
    EXPECT_EQ(index[std::make_tuple(0,0,2,0,0)], 8); // Si, 1st atom, 1d, m = 0
    EXPECT_EQ(index[std::make_tuple(0,0,2,0,1)], 9); // Si, 1st atom, 1d, m = 1
    EXPECT_EQ(index[std::make_tuple(0,0,2,0,-1)], 10); // Si, 1st atom, 1d, m = -1
    EXPECT_EQ(index[std::make_tuple(0,0,2,0,2)], 11); // Si, 1st atom, 1d, m = 2
    EXPECT_EQ(index[std::make_tuple(0,0,2,0,-2)], 12); // Si, 1st atom, 1d, m = -2
    // it = 1
    // ia = 0
    // l = 0
    // izeta = 0
    EXPECT_EQ(index[std::make_tuple(1,0,0,0,0)], 13); // C, 1st atom, 1s
    // izeta = 1
    EXPECT_EQ(index[std::make_tuple(1,0,0,1,0)], 14); // C, 1st atom, 2s
    // l = 1
    // izeta = 0
    EXPECT_EQ(index[std::make_tuple(1,0,1,0,0)], 15); // C, 1st atom, 1p, m = 0
    EXPECT_EQ(index[std::make_tuple(1,0,1,0,1)], 16); // C, 1st atom, 1p, m = 1
    EXPECT_EQ(index[std::make_tuple(1,0,1,0,-1)], 17); // C, 1st atom, 1p, m = -1
    // izeta = 1
    EXPECT_EQ(index[std::make_tuple(1,0,1,1,0)], 18); // C, 1st atom, 2p, m = 0
    EXPECT_EQ(index[std::make_tuple(1,0,1,1,1)], 19); // C, 1st atom, 2p, m = 1
    EXPECT_EQ(index[std::make_tuple(1,0,1,1,-1)], 20); // C, 1st atom, 2p, m = -1
    // l = 2
    // izeta = 0
    EXPECT_EQ(index[std::make_tuple(1,0,2,0,0)], 21); // C, 1st atom, 1d, m = 0
    EXPECT_EQ(index[std::make_tuple(1,0,2,0,1)], 22); // C, 1st atom, 1d, m = 1
    EXPECT_EQ(index[std::make_tuple(1,0,2,0,-1)], 23); // C, 1st atom, 1d, m = -1
    EXPECT_EQ(index[std::make_tuple(1,0,2,0,2)], 24); // C, 1st atom, 1d, m = 2
    EXPECT_EQ(index[std::make_tuple(1,0,2,0,-2)], 25); // C, 1st atom, 1d, m = -2
    // reverse
    EXPECT_EQ(index_reverse[0], std::make_tuple(0,0,0,0,0)); // Si, 1st atom, 1s
    EXPECT_EQ(index_reverse[1], std::make_tuple(0,0,0,1,0)); // Si, 1st atom, 2s
    EXPECT_EQ(index_reverse[2], std::make_tuple(0,0,1,0,0)); // Si, 1st atom, 1p, m = 0
    EXPECT_EQ(index_reverse[3], std::make_tuple(0,0,1,0,1)); // Si, 1st atom, 1p, m = 1
    EXPECT_EQ(index_reverse[4], std::make_tuple(0,0,1,0,-1)); // Si, 1st atom, 1p, m = -1
    EXPECT_EQ(index_reverse[5], std::make_tuple(0,0,1,1,0)); // Si, 1st atom, 2p, m = 0
    EXPECT_EQ(index_reverse[6], std::make_tuple(0,0,1,1,1)); // Si, 1st atom, 2p, m = 1
    EXPECT_EQ(index_reverse[7], std::make_tuple(0,0,1,1,-1)); // Si, 1st atom, 2p, m = -1
    EXPECT_EQ(index_reverse[8], std::make_tuple(0,0,2,0,0)); // Si, 1st atom, 1d, m = 0
    EXPECT_EQ(index_reverse[9], std::make_tuple(0,0,2,0,1)); // Si, 1st atom, 1d, m = 1
    EXPECT_EQ(index_reverse[10], std::make_tuple(0,0,2,0,-1)); // Si, 1st atom, 1d, m = -1
    EXPECT_EQ(index_reverse[11], std::make_tuple(0,0,2,0,2)); // Si, 1st atom, 1d, m = 2
    EXPECT_EQ(index_reverse[12], std::make_tuple(0,0,2,0,-2)); // Si, 1st atom, 1d, m = -2
    EXPECT_EQ(index_reverse[13], std::make_tuple(1,0,0,0,0)); // C, 1st atom, 1s
    EXPECT_EQ(index_reverse[14], std::make_tuple(1,0,0,1,0)); // C, 1st atom, 2s
    EXPECT_EQ(index_reverse[15], std::make_tuple(1,0,1,0,0)); // C, 1st atom, 1p, m = 0
    EXPECT_EQ(index_reverse[16], std::make_tuple(1,0,1,0,1)); // C, 1st atom, 1p, m = 1
    EXPECT_EQ(index_reverse[17], std::make_tuple(1,0,1,0,-1)); // C, 1st atom, 1p, m = -1
    EXPECT_EQ(index_reverse[18], std::make_tuple(1,0,1,1,0)); // C, 1st atom, 2p, m = 0
    EXPECT_EQ(index_reverse[19], std::make_tuple(1,0,1,1,1)); // C, 1st atom, 2p, m = 1
    EXPECT_EQ(index_reverse[20], std::make_tuple(1,0,1,1,-1)); // C, 1st atom, 2p, m = -1
    EXPECT_EQ(index_reverse[21], std::make_tuple(1,0,2,0,0)); // C, 1st atom, 1d, m = 0
    EXPECT_EQ(index_reverse[22], std::make_tuple(1,0,2,0,1)); // C, 1st atom, 1d, m = 1
    EXPECT_EQ(index_reverse[23], std::make_tuple(1,0,2,0,-1)); // C, 1st atom, 1d, m = -1
    EXPECT_EQ(index_reverse[24], std::make_tuple(1,0,2,0,2)); // C, 1st atom, 1d, m = 2
    EXPECT_EQ(index_reverse[25], std::make_tuple(1,0,2,0,-2)); // C, 1st atom, 1d, m = -2
    // test2 2Si, 3C
    natoms = {2, 3};
    tqo.radialcollection_indexing(*(tqo.p_nao()), natoms, index, index_reverse);
    EXPECT_EQ(index.size(), 65); // (2*1 + 2*3 + 1*5)*2 + (2*1 + 2*3 + 1*5)*3 = 65
    EXPECT_EQ(index_reverse.size(), 65);
    // it, ia, l, izeta, m
    EXPECT_EQ(index[std::make_tuple(0,0,0,0,0)], 0); // Si, 1st atom, 1s
    EXPECT_EQ(index[std::make_tuple(0,0,2,0,-2)], 12); // Si, 1st atom, 1d, m = -2
    EXPECT_EQ(index[std::make_tuple(0,1,0,0,0)], 13); // Si, 2nd atom, 1s
    EXPECT_EQ(index[std::make_tuple(0,1,2,0,-2)], 25); // Si, 2nd atom, 1d, m = -2
    EXPECT_EQ(index[std::make_tuple(1,0,0,0,0)], 26); // C, 1st atom, 1s
    EXPECT_EQ(index[std::make_tuple(1,0,2,0,-2)], 38); // C, 1st atom, 1d, m = -2
    EXPECT_EQ(index[std::make_tuple(1,1,0,0,0)], 39); // C, 2nd atom, 1s
    EXPECT_EQ(index[std::make_tuple(1,1,2,0,-2)], 51); // C, 2nd atom, 1d, m = -2
    EXPECT_EQ(index[std::make_tuple(1,2,0,0,0)], 52); // C, 3rd atom, 1s
    EXPECT_EQ(index[std::make_tuple(1,2,2,0,-2)], 64); // C, 3rd atom, 1d, m = -2
    // reverse
    EXPECT_EQ(index_reverse[0], std::make_tuple(0,0,0,0,0)); // Si, 1st atom, 1s
    EXPECT_EQ(index_reverse[12], std::make_tuple(0,0,2,0,-2)); // Si, 1st atom, 1d, m = -2
    EXPECT_EQ(index_reverse[13], std::make_tuple(0,1,0,0,0)); // Si, 2nd atom, 1s
    EXPECT_EQ(index_reverse[25], std::make_tuple(0,1,2,0,-2)); // Si, 2nd atom, 1d, m = -2
    EXPECT_EQ(index_reverse[26], std::make_tuple(1,0,0,0,0)); // C, 1st atom, 1s
    EXPECT_EQ(index_reverse[38], std::make_tuple(1,0,2,0,-2)); // C, 1st atom, 1d, m = -2
    EXPECT_EQ(index_reverse[39], std::make_tuple(1,1,0,0,0)); // C, 2nd atom, 1s
    EXPECT_EQ(index_reverse[51], std::make_tuple(1,1,2,0,-2)); // C, 2nd atom, 1d, m = -2
    EXPECT_EQ(index_reverse[52], std::make_tuple(1,2,0,0,0)); // C, 3rd atom, 1s
    EXPECT_EQ(index_reverse[64], std::make_tuple(1,2,2,0,-2)); // C, 3rd atom, 1d, m = -2
}

TEST_F(toQOTest, BuildHydrogenMinimal)
{
    define_fcc_cell(ucell);
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
    tqo.unwrap_unitcell(&ucell);
    tqo.build_ao(ucell.ntype,
                 "./",
                 ucell.pseudo_fn,
                 {},
                 GlobalV::qo_thr,
                 GlobalV::ofs_running,
                 0);
    EXPECT_EQ(tqo.p_ao()->nchi(), 5); // Si: 1s, 2p, 3d, C: 1s, 2p
    EXPECT_EQ(tqo.nchi(), 13); // Si: 1s, 2px, 2py, 2pz, 3dz2, 3dxz, 3dyz, 3dx2-y2, 3dxy, C: 1s, 2px, 2py, 2pz
    tqo.p_ao()->to_file("special_use_unittest");
}

// the scan_supercell_for_atom() calls
TEST_F(toQOTest, Norm2RijSupercell)
{
    define_fcc_cell(ucell);
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
    tqo.unwrap_unitcell(&ucell);
    ModuleBase::Vector3<double> rij(1.0, 0.0, 0.0);
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 0, 0, 0), 1.0); // R = 0, 0, 0
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 1, 0, 0), 145.0);
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 0, 1, 0), 145.0);
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 0, 0, 1), 129.0);
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 1, 1, 0), 417.0);
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 1, 0, 1), 401.0);
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 0, 1, 1), 401.0);
    EXPECT_EQ(tqo.norm2_rij_supercell(rij, 1, 1, 1), 801.0);
}
// the scan_supercell() calls
TEST_F(toQOTest, ScanSupercellForAtom)
{
    define_fcc_cell(ucell);
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype,     // ntype
                  "./",            // orbital_dir
                  ucell.orbital_fn,// orbital_fn
                  0);              // rank
    std::vector<int> nmax = std::vector<int>(ucell.ntype);
    for(int itype = 0; itype < ucell.ntype; itype++)
    {
        nmax[itype] = tqo.atom_database().principle_quantum_number[tqo.symbols()[itype]];
    }
    tqo.build_ao(ucell.ntype,
                 "./",
                 ucell.pseudo_fn,
                 {},
                 GlobalV::qo_thr,
                 GlobalV::ofs_running,
                 0);
    std::vector<ModuleBase::Vector3<int>> n1n2n3 = tqo.scan_supercell_for_atom(0, 0);
    EXPECT_EQ(n1n2n3.size(), 13); // 13 = 3*3*3 - 2 - 3*4
}
// the scan_supercell() calls
TEST_F(toQOTest, EliminateDuplicateVector3)
{
    define_fcc_cell(ucell);
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
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
    tqo.eliminate_duplicate_vector3<int>(v);
    EXPECT_EQ(v.size(), 4);
}

TEST_F(toQOTest, ScanSupercellFCC)
{
    define_fcc_cell(ucell);
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, 
                  "./",
                  ucell.orbital_fn,
                  0);
    tqo.build_ao(ucell.ntype,
                 "./",
                 ucell.pseudo_fn,
                 {},
                 GlobalV::qo_thr,
                 GlobalV::ofs_running,
                 0);
    tqo.scan_supercell();
    EXPECT_EQ(tqo.nR(), 13);
}

TEST_F(toQOTest, ScanSupercellSC1)
{
    define_sc_cell(ucell);
    toQO tqo("hydrogen", {"minimal-nodeless"});
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, 
                  "./",
                  ucell.orbital_fn,
                  0);
    GlobalV::qo_thr = 1e-6;
    tqo.build_ao(ucell.ntype,
                 "./",
                 ucell.pseudo_fn,
                 {},
                 GlobalV::qo_thr,
                 GlobalV::ofs_running,
                 0);
    tqo.scan_supercell();
    EXPECT_EQ(tqo.nR(), 19); // 3*3*3 - 8 (corner 111, -1-1-1, etc)
}

TEST_F(toQOTest, AllocateOvlpMinimal)
{
    define_fcc_cell(ucell);
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, 
                  "./",
                  ucell.orbital_fn,
                  0);
    std::vector<int> nmax = std::vector<int>(ucell.ntype);
    for(int itype = 0; itype < ucell.ntype; itype++)
    {
        nmax[itype] = tqo.atom_database().principle_quantum_number[tqo.symbols()[itype]];
    }
    tqo.build_ao(ucell.ntype,
                 "./",
                 ucell.pseudo_fn,
                 {},
                 GlobalV::qo_thr,
                 GlobalV::ofs_running,
                 0);
    tqo.scan_supercell();
    tqo.allocate_ovlp(true);
    tqo.allocate_ovlp(false);
    EXPECT_EQ(tqo.ovlpk().size(), tqo.nchi()*tqo.nphi()); // for single kpoint, ao*nao matrix
    EXPECT_EQ(tqo.ovlpR().size(), tqo.nchi()*tqo.nphi()); // for single cell, ao*nao matrix
    // all values in them should be zero or complex zero
    for(int i = 0; i < tqo.nchi(); i++)
    {
        for(int j = 0; j < tqo.nphi(); j++)
        {
            EXPECT_EQ(tqo.ovlpR(i, j), 0.0);
            EXPECT_EQ(tqo.ovlpk(i, j), std::complex<double>(0.0, 0.0));
        }
    }
}

TEST_F(toQOTest, Initialize)
{
    define_fcc_cell(ucell);
    GlobalV::qo_screening_coeff = {};
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
    std::vector<ModuleBase::Vector3<double>> kvecs_c;
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0)); // Gamma point
    tqo.initialize(&ucell, kvecs_c);
}

TEST_F(toQOTest, CalculateOvlpR)
{
    define_fcc_cell(ucell);
    GlobalV::qo_screening_coeff = {};
    GlobalV::qo_thr = 1e-10;
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
    std::vector<ModuleBase::Vector3<double>> kvecs_c;
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0)); // Gamma point
    tqo.initialize(&ucell, kvecs_c);
    // find the R = 0,0,0
    for(int iR = 0; iR < tqo.nR(); iR++)
    {
        if(tqo.supercells()[iR].x == 0 && tqo.supercells()[iR].y == 0 && tqo.supercells()[iR].z == 0)
        {
            tqo.calculate_ovlpR(iR);
            break;
        }
    }
    int nrows = tqo.nchi();
    int ncols = tqo.nphi();
    // not all elements are zero
    bool all_zero = true;
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            if(tqo.ovlpR(i, j) != 0.0)
            {
                all_zero = false;
            }
        }
    }
    EXPECT_EQ(all_zero, false);
}

TEST_F(toQOTest, CalculateSelfOvlpRMinimal)
{
    define_fcc_cell(ucell);
    GlobalV::qo_screening_coeff = {};
    GlobalV::qo_thr = 1e-10;
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
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
            tqo.calculate_ovlpR(iR);
            break;
        }
    }
    // check if diagonal elements are 1
    for(int i = 0; i < tqo.nphi(); i++)
    {
        EXPECT_NEAR(tqo.ovlpR(i, i), 1.0, 1e-3); // this is too tight for 1s orbital, which fluctuates a lot in narrow region
    }
    //std::remove("Si_special_use_unittest.orb");
    //std::remove("C_special_use_unittest.orb");
    //tqo.write_ovlp(tqo.ovlpR()[0], "QO_self_ovlp.dat");
}

TEST_F(toQOTest, CalculateSelfOvlpKSymmetrical)
{
    define_fcc_cell(ucell);
    GlobalV::qo_thr = 1e-10;
    GlobalV::qo_screening_coeff = {};
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
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
    int nrows = tqo.nchi();
    int ncols = tqo.nphi();

    tqo.calculate_ovlpk(0);
    std::vector<std::complex<double>> ovlpk_1 = tqo.ovlpk();
    tqo.calculate_ovlpk(1);
    std::vector<std::complex<double>> ovlpk_2 = tqo.ovlpk();
    bool all_zero = true;
    // check if all imaginary parts are cancelled
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlpR_ij = ovlpk_1[i*nrows+j] + ovlpk_2[i*nrows+j];
            EXPECT_NEAR(ovlpR_ij.imag(), 0.0, 1e-10);
            if(ovlpR_ij.real() > 1e-10)
            {
                all_zero = false;
            }
        }
    }
    EXPECT_FALSE(all_zero);
    // test symmetry cancellation on pair2
    tqo.calculate_ovlpk(2);
    std::vector<std::complex<double>> ovlpk_3 = tqo.ovlpk();
    tqo.calculate_ovlpk(3);
    std::vector<std::complex<double>> ovlpk_4 = tqo.ovlpk();
    all_zero = true;
    // check if all imaginary parts are cancelled
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlpR_ij = ovlpk_3[i*nrows+j] + ovlpk_4[i*nrows+j];
            EXPECT_NEAR(ovlpR_ij.imag(), 0.0, 1e-10);
            if(ovlpR_ij.real() > 1e-10)
            {
                all_zero = false;
            }
        }
    }
    EXPECT_FALSE(all_zero);
    // test symmetry cancellation on pair3
    tqo.calculate_ovlpk(4);
    std::vector<std::complex<double>> ovlpk_5 = tqo.ovlpk();
    tqo.calculate_ovlpk(5);
    std::vector<std::complex<double>> ovlpk_6 = tqo.ovlpk();
    all_zero = true;
    // check if all imaginary parts are cancelled
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlpR_ij = ovlpk_5[i*nrows+j] + ovlpk_6[i*nrows+j];
            EXPECT_NEAR(ovlpR_ij.imag(), 0.0, 1e-10);
            if(ovlpR_ij.real() > 1e-10)
            {
                all_zero = false;
            }
        }
    }
    EXPECT_FALSE(all_zero);
    // test symmetry cancellation on pair4
    tqo.calculate_ovlpk(6);
    std::vector<std::complex<double>> ovlpk_7 = tqo.ovlpk();
    tqo.calculate_ovlpk(7);
    std::vector<std::complex<double>> ovlpk_8 = tqo.ovlpk();
    all_zero = true;
    // check if all imaginary parts are cancelled
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlpR_ij = ovlpk_7[i*nrows+j] + ovlpk_8[i*nrows+j];
            EXPECT_NEAR(ovlpR_ij.imag(), 0.0, 1e-10);
            if(ovlpR_ij.real() > 1e-10)
            {
                all_zero = false;
            }
        }
    }
    EXPECT_FALSE(all_zero);
    // test symmetry cancellation on pair5
    tqo.calculate_ovlpk(8);
    std::vector<std::complex<double>> ovlpk_9 = tqo.ovlpk();
    all_zero = true;
    // check if all imaginary parts are cancelled
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            EXPECT_NEAR(ovlpk_9[i*nrows+j].imag(), 0.0, 1e-10);
            if(ovlpk_9[i*nrows+j].real() > 1e-10)
            {
                all_zero = false;
            }
        }
    }
    EXPECT_FALSE(all_zero);
    std::remove("Si_special_use_unittest.orb");
    std::remove("C_special_use_unittest.orb");
    //tqo.write_ovlp(tqo.ovlpR()[0], "QO_self_ovlp.dat");
}

TEST_F(toQOTest, BuildHydrogenFull)
{
    define_fcc_cell(ucell);
    GlobalV::qo_thr = 1e-10;
    toQO tqo("hydrogen", {"full", "full"});
    tqo.unwrap_unitcell(&ucell);
    GlobalV::qo_thr = 1e-10;
    tqo.build_ao(ucell.ntype,
                 "./",
                 ucell.pseudo_fn,
                 {},
                 GlobalV::qo_thr,
                 GlobalV::ofs_running,
                 0);
    EXPECT_EQ(tqo.p_ao()->nchi(), 9); // Si: 1s, 2s, 2p, 3s, 3p, 3d, C: 1s, 2s, 2p
    EXPECT_EQ(tqo.nchi(), 19); 
    tqo.p_ao()->to_file("special_use_unittest");
}

TEST_F(toQOTest, CalculateSelfOvlpRFull)
{
    define_fcc_cell(ucell);
    GlobalV::qo_thr = 1e-10;
    GlobalV::qo_screening_coeff = {};
    toQO tqo("hydrogen", {"full", "full"});
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
            tqo.calculate_ovlpR(iR);
            break;
        }
    }
    // check if diagonal elements are 1
    for(int i = 0; i < tqo.nphi(); i++)
    {
        EXPECT_NEAR(tqo.ovlpR(i, i), 1.0, 5e-4); // this is too tight for 1s orbital, which fluctuates a lot in narrow region
    }
    // check if symmetrical
    for(int i = 0; i < tqo.nchi(); i++)
    {
        for(int j = 0; j < tqo.nphi(); j++)
        {
            EXPECT_NEAR(tqo.ovlpR(i, j), tqo.ovlpR(j, i), 1e-8);
        }
    }
    std::remove("Si_special_use_unittest.orb");
    std::remove("C_special_use_unittest.orb");
    //tqo.write_ovlp(tqo.ovlpR()[0], "QO_self_ovlp.dat");
}

/* Si_dojo_soc.upf is special: two p orbitals, one s orbital */

TEST_F(toQOTest, BuildPswfcPartial1)
{
    define_fcc_cell(ucell);
    toQO tqo("pswfc", {"s", "s"});
    tqo.unwrap_unitcell(&ucell);
    tqo.build_ao(ucell.ntype,
                 "./",
                 ucell.pseudo_fn,
                 GlobalV::qo_screening_coeff,
                 GlobalV::qo_thr,
                 GlobalV::ofs_running,
                 0);
    EXPECT_EQ(tqo.p_ao()->nchi(), 5); // AO will always read and import all orbitals
    EXPECT_EQ(tqo.nchi(), 2);
}

TEST_F(toQOTest, BuildPswfcPartial2)
{
    define_fcc_cell(ucell);
    toQO tqo("pswfc", {"ps", "s"});
    tqo.unwrap_unitcell(&ucell);
    tqo.build_ao(ucell.ntype,
                 "./",
                 ucell.pseudo_fn,
                 GlobalV::qo_screening_coeff,
                 GlobalV::qo_thr,
                 GlobalV::ofs_running,
                 0);
    EXPECT_EQ(tqo.p_ao()->nchi(), 5); // AO will always read and import all orbitals
    EXPECT_EQ(tqo.nchi(), 8); // the first element is Si, it has two p orbitals, so 3+3+1+1
}

TEST_F(toQOTest, BuildPswfcPartial3)
{
    define_fcc_cell(ucell);
    toQO tqo("pswfc", {"all", "p"});
    tqo.unwrap_unitcell(&ucell);
    tqo.build_ao(ucell.ntype,
                 "./",
                 ucell.pseudo_fn,
                 GlobalV::qo_screening_coeff,
                 GlobalV::qo_thr,
                 GlobalV::ofs_running,
                 0);
    EXPECT_EQ(tqo.p_ao()->nchi(), 5); // AO will always read and import all orbitals
    EXPECT_EQ(tqo.nchi(), 10);
}

TEST_F(toQOTest, BuildPswfcAll)
{
    define_fcc_cell(ucell);
    GlobalV::qo_thr = 1e-10;
    toQO tqo("pswfc", {"all", "all"});
    tqo.unwrap_unitcell(&ucell);
    tqo.build_ao(ucell.ntype,
                 "./",
                 ucell.pseudo_fn,
                 GlobalV::qo_screening_coeff,
                 GlobalV::qo_thr,
                 GlobalV::ofs_running,
                 0);
    EXPECT_EQ(tqo.p_ao()->nchi(), 5); 
    EXPECT_EQ(tqo.nchi(), 11);
    tqo.p_ao()->to_file("special_use_unittest");
}

TEST_F(toQOTest, ScanSupercellSC2)
{
    define_sc_cell(ucell);
    toQO tqo("pswfc", {"all"});
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, 
                  "./",
                  ucell.orbital_fn,
                  0);
    GlobalV::qo_screening_coeff[0] = 0.1; // use this to control the tailing of radial function
    GlobalV::qo_thr = 1e-6;
    tqo.build_ao(ucell.ntype,
                 "./",
                 ucell.pseudo_fn,
                 GlobalV::qo_screening_coeff,
                 GlobalV::qo_thr,
                 GlobalV::ofs_running,
                 0); // radius = 13.6 Bohr
    tqo.scan_supercell();
    EXPECT_EQ(tqo.nR(), 81); // 5*5*5 - 12(edge center) - 8*4(corner)
}

TEST_F(toQOTest, ScanSupercellSC3)
{
    define_sc_cell(ucell);
    toQO tqo("pswfc", {"all"});
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, 
                  "./",
                  ucell.orbital_fn,
                  0);
    GlobalV::qo_screening_coeff[0] = 0.25; // use this to control the tailing of radial function
    GlobalV::qo_thr = 1e-6;
    tqo.build_ao(ucell.ntype,
                 "./",
                 ucell.pseudo_fn,
                 GlobalV::qo_screening_coeff,
                 GlobalV::qo_thr,
                 GlobalV::ofs_running,
                 0); // radius = 13.6 Bohr
    tqo.scan_supercell();
    EXPECT_EQ(tqo.nR(), 57); // 5*5*5 - 12(edge center) - 8*(8-1)(corner) = 5*5*5 - 12(edge center) - 8*(2*2*2-1)(corner)
    GlobalV::qo_screening_coeff[0] = 0.1;
}

TEST_F(toQOTest, ScanSupercellSC4)
{
    define_sc_cell(ucell);
    toQO tqo("pswfc", {"all"});
    tqo.unwrap_unitcell(&ucell);
    tqo.build_nao(ucell.ntype, 
                  "./",
                  ucell.orbital_fn,
                  0);
    GlobalV::qo_screening_coeff[0] = 0.5; // use this to control the tailing of radial function
    GlobalV::qo_thr = 1e-6;
    tqo.build_ao(ucell.ntype,
                 "./",
                 ucell.pseudo_fn,
                 GlobalV::qo_screening_coeff,
                 GlobalV::qo_thr,
                 GlobalV::ofs_running,
                 0); // radius = 13.6 Bohr
    tqo.scan_supercell();
    EXPECT_EQ(tqo.nR(), 33); // 3*3*3 + 6(face)
    GlobalV::qo_screening_coeff[0] = 0.1;
}

TEST_F(toQOTest, CalculateSelfOvlpRPswfc)
{
    define_fcc_cell(ucell);
    GlobalV::qo_thr = 1e-10;
    toQO tqo("pswfc", {"all", "all"});
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
            tqo.calculate_ovlpR(iR);
            break;
        }
    }
    // check if diagonal elements are 1
    for(int i = 0; i < tqo.nphi(); i++)
    {
        EXPECT_NEAR(tqo.ovlpR(i, i), 1.0, 1e-4);
    }
    // check if symmetrical
    for(int i = 0; i < tqo.nchi(); i++)
    {
        for(int j = 0; j < tqo.nphi(); j++)
        {
            EXPECT_NEAR(tqo.ovlpR(i, j), tqo.ovlpR(j, i), 1e-4);
        }
    }
    std::remove("Si_special_use_unittest.orb");
    std::remove("C_special_use_unittest.orb");
}

TEST_F(toQOTest, CalculateOvlpKGamma)
{
    define_fcc_cell(ucell);
    GlobalV::qo_thr = 1e-10;
    GlobalV::qo_screening_coeff = {};
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
    std::vector<ModuleBase::Vector3<double>> kvecs_c;
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0)); // Gamma point
    tqo.initialize(&ucell, kvecs_c);
    tqo.calculate_ovlpk(0);
    int nrows = tqo.nchi();
    int ncols = tqo.nphi();
    // all should be real numbers at Gamma point
    bool all_real = true;
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            if(tqo.ovlpk()[i*nrows+j].imag() != 0.0)
            {
                all_real = false;
            }
        }
    }
    EXPECT_TRUE(all_real);
}

TEST_F(toQOTest, CalculateOvlpKSlaterGamma)
{
    define_fcc_cell(ucell);
    GlobalV::qo_thr = 1e-10;
    GlobalV::qo_screening_coeff = {0.1};
    toQO tqo("hydrogen", {"energy-full", "energy-full"});
    std::vector<ModuleBase::Vector3<double>> kvecs_c;
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0)); // Gamma point
    tqo.initialize(&ucell, kvecs_c);
    tqo.calculate_ovlpk(0);
    int nrows = tqo.nchi();
    int ncols = tqo.nphi();
    // all should be real numbers at Gamma point
    bool all_real = true;
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            if(tqo.ovlpk()[i*nrows+j].imag() != 0.0)
            {
                all_real = false;
            }
        }
    }
    EXPECT_TRUE(all_real);
}

TEST_F(toQOTest, CalculateSelfOvlpKPswfcSymmetrical)
{
    define_fcc_cell(ucell);
    GlobalV::qo_thr = 1e-10;
    GlobalV::qo_screening_coeff = {2.0, 2.0};
    toQO tqo("pswfc", {"all", "all"});
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
    int nrows = tqo.nchi();
    int ncols = tqo.nphi();
    std::cout << "Number of supercells: " << tqo.nR() << ", number of kpoints: " << tqo.nkpts() << std::endl;
    tqo.calculate_ovlpk(0);
    std::vector<std::complex<double>> ovlpk_1 = tqo.ovlpk();
    tqo.calculate_ovlpk(1);
    std::vector<std::complex<double>> ovlpk_2 = tqo.ovlpk();
    // check if all imaginary parts are cancelled
    bool all_real = true;
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlpR_ij = ovlpk_1[i*nrows+j] + ovlpk_2[i*nrows+j];
            EXPECT_NEAR(ovlpR_ij.imag(), 0.0, 1e-8);
        }
    }
    tqo.calculate_ovlpk(2);
    std::vector<std::complex<double>> ovlpk_3 = tqo.ovlpk();
    tqo.calculate_ovlpk(3);
    std::vector<std::complex<double>> ovlpk_4 = tqo.ovlpk();
    // check if all imaginary parts are cancelled
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlpR_ij = ovlpk_3[i*nrows+j] + ovlpk_4[i*nrows+j];
            EXPECT_NEAR(ovlpR_ij.imag(), 0.0, 1e-8);
        }
    }
    tqo.calculate_ovlpk(4);
    std::vector<std::complex<double>> ovlpk_5 = tqo.ovlpk();
    tqo.calculate_ovlpk(5);
    std::vector<std::complex<double>> ovlpk_6 = tqo.ovlpk();
    // check if all imaginary parts are cancelled
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlpR_ij = ovlpk_5[i*nrows+j] + ovlpk_6[i*nrows+j];
            EXPECT_NEAR(ovlpR_ij.imag(), 0.0, 1e-8);
        }
    }
    tqo.calculate_ovlpk(6);
    std::vector<std::complex<double>> ovlpk_7 = tqo.ovlpk();
    tqo.calculate_ovlpk(7);
    std::vector<std::complex<double>> ovlpk_8 = tqo.ovlpk();
    // check if all imaginary parts are cancelled
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            std::complex<double> ovlpR_ij = ovlpk_7[i*nrows+j] + ovlpk_8[i*nrows+j];
            EXPECT_NEAR(ovlpR_ij.imag(), 0.0, 1e-8);
        }
    }
    tqo.calculate_ovlpk(8);
    std::vector<std::complex<double>> ovlpk_9 = tqo.ovlpk();
    // check if all imaginary parts are cancelled
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            // R = 0, 0, 0, then unfolding kphase would be e-ikR = 1,
            // becomes direct summation over kpoints
            EXPECT_NEAR(ovlpk_9[i*nrows+j].imag(), 0.0, 1e-8);
        }
    }
    //tqo.write_ovlp(tqo.ovlpR()[0], "QO_self_ovlp.dat");
}

TEST_F(toQOTest, CalculateHydrogenlike)
{
    define_fcc_cell(ucell);
    GlobalV::qo_thr = 1e-10;
    GlobalV::qo_screening_coeff = {};
    toQO tqo("hydrogen", {"minimal-nodeless", "minimal-nodeless"});
    std::vector<ModuleBase::Vector3<double>> kvecs_c;
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.0, 0.0, 0.0)); // Gamma point
    kvecs_c.push_back(ModuleBase::Vector3<double>(0.5, 0.0, 0.0));
    tqo.initialize(&ucell, kvecs_c);
    tqo.calculate();
    int nrows = tqo.nchi();
    int ncols = tqo.nphi();
    // for the latest kpoint, not all numbers are complex zero
    bool all_zero = true;
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            if(tqo.ovlpk()[i*nrows+j] != std::complex<double>(0.0, 0.0))
            {
                all_zero = false;
            }
        }
    }
    EXPECT_EQ(all_zero, false);
    // delete files generated namely QO_ovlp_0.dat and QO_ovlp_1.dat
    std::remove("QO_ovlp_0.dat");
    std::remove("QO_ovlp_1.dat");
    std::remove("QO_supercells.dat");
}

/**/
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
