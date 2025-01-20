#include <gtest/gtest.h>
#define private public
#include "module_parameter/parameter.h"
#undef private
#include "../psi_initializer.h"
#include "../psi_initializer_atomic.h"
#include "../psi_initializer_atomic_random.h"
#include "../psi_initializer_nao.h"
#include "../psi_initializer_nao_random.h"
#include "../psi_initializer_random.h"
#include "module_hamilt_pw/hamilt_pwdft/VL_in_pw.h"
#include "module_cell/klist.h"

/*
=========================
psi initializer unit test
=========================
- Tested functions:
    - psi_initializer_random::psi_initializer_random
      - constructor of psi_initializer_random
    - psi_initializer_atomic::psi_initializer_atomic
      - constructor of psi_initializer_atomic
    - psi_initializer_atomic_random::psi_initializer_atomic_random
      - constructor of psi_initializer_atomic_random
    - psi_initializer_nao::psi_initializer_nao
      - constructor of psi_initializer_nao
    - psi_initializer_nao_random::psi_initializer_nao_random
      - constructor of psi_initializer_nao_random
    - psi_initializer::cast_to_T (psi_initializer specialized as random)
      - function cast std::complex<double> to float, double, std::complex<float>, std::complex<double>
    - psi_initializer_random::allocate
      - allocate wavefunctions with random-specific method
    - psi_initializer_atomic::allocate
      - allocate wavefunctions with atomic-specific method
    - psi_initializer_atomic_random::allocate
      - allocate wavefunctions with atomic-specific method
    - psi_initializer_nao::allocate
      - allocate wavefunctions with nao-specific method
    - psi_initializer_nao_random::allocate
      - allocate wavefunctions with nao-specific method
    - psi_initializer_random::proj_ao_onkG
      - calculate wavefunction initial guess (before diagonalization) by randomly generating numbers
    - psi_initializer_atomic::proj_ao_onkG
      - calculate wavefunction initial guess (before diagonalization) with atomic pseudo wavefunctions
      - nspin = 4 case
      - nspin = 4 with has_so case
    - psi_initializer_atomic_random::proj_ao_onkG
      - calculate wavefunction initial guess (before diagonalization) with atomic pseudo wavefunctions and random numbers
    - psi_initializer_nao::proj_ao_onkG
      - calculate wavefunction initial guess (before diagonalization) with numerical atomic orbital wavefunctions
    - psi_initializer_nao_random::proj_ao_onkG
      - calculate wavefunction initial guess (before diagonalization) with numerical atomic orbital wavefunctions and random numbers
*/

// here are many empty functions because we have included many headers in file to be tested
// but it does not mean all functions can only be defined for only once in those corresponding files
// so here we define them again to avoid undefined reference error
Atom_pseudo::Atom_pseudo() {}
Atom_pseudo::~Atom_pseudo() {}
#ifdef __MPI
void Atom_pseudo::bcast_atom_pseudo() {}
#endif
pseudo::pseudo() {}
pseudo::~pseudo() {}

pseudopot_cell_vnl::pseudopot_cell_vnl() {}
pseudopot_cell_vnl::~pseudopot_cell_vnl()
{
}
pseudopot_cell_vl::pseudopot_cell_vl() {}
pseudopot_cell_vl::~pseudopot_cell_vl() {}
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}
void output::printM3(std::ofstream &ofs, const std::string &description, const ModuleBase::Matrix3 &m) {}
#ifdef __LCAO
ORB_gaunt_table::ORB_gaunt_table() {}
ORB_gaunt_table::~ORB_gaunt_table() {}
InfoNonlocal::InfoNonlocal() {}
InfoNonlocal::~InfoNonlocal() {}
#endif
Structure_Factor::Structure_Factor() {}
Structure_Factor::~Structure_Factor() {}
void Structure_Factor::setup_structure_factor(const UnitCell* Ucell, const Parallel_Grid&, const ModulePW::PW_Basis* rho_basis) {}
std::complex<double>* Structure_Factor::get_sk(int ik, int it, int ia, ModulePW::PW_Basis_K const*wfc_basis) const
{
    int npw = wfc_basis->npwk[ik];
    std::complex<double> *sk = new std::complex<double>[npw];
    for(int ipw = 0; ipw < npw; ++ipw) { sk[ipw] = std::complex<double>(0.0, 0.0);
}
    return sk;
}

class PsiIntializerUnitTest : public ::testing::Test {
    public:
        Structure_Factor* p_sf = nullptr;
        ModulePW::PW_Basis_K* p_pw_wfc = nullptr;
        UnitCell* p_ucell = nullptr;
        pseudopot_cell_vnl* p_pspot_vnl = nullptr;
        K_Vectors* p_kv = nullptr;
        int random_seed = 1;

        psi_initializer<std::complex<double>>* psi_init;

      private:
      protected:
        void SetUp() override
        {
            // allocate
            this->p_sf = new Structure_Factor();
            this->p_pw_wfc = new ModulePW::PW_Basis_K();
            this->p_ucell = new UnitCell();
            this->p_pspot_vnl = new pseudopot_cell_vnl();
            this->p_kv = new K_Vectors();
            // mock
            PARAM.input.nbands = 1;
            PARAM.input.nspin = 1;
            PARAM.input.orbital_dir = "./support/";
            PARAM.input.pseudo_dir = "./support/";
            PARAM.sys.npol = 1;
            PARAM.input.calculation = "scf";
            PARAM.input.init_wfc = "random";
            PARAM.input.ks_solver = "cg";
            PARAM.sys.domag = false;
            PARAM.sys.domag_z = false;
            // lattice
            this->p_ucell->a1 = {10.0, 0.0, 0.0};
            this->p_ucell->a2 = {0.0, 10.0, 0.0};
            this->p_ucell->a3 = {0.0, 0.0, 10.0};
            this->p_ucell->lat0 = 1.0;
            this->p_ucell->omega = 1000.0;
            this->p_ucell->latvec.e11 = 10.0; this->p_ucell->latvec.e12 = 0.0; this->p_ucell->latvec.e13 = 0.0;
            this->p_ucell->latvec.e21 = 0.0; this->p_ucell->latvec.e22 = 10.0; this->p_ucell->latvec.e23 = 0.0;
            this->p_ucell->latvec.e31 = 0.0; this->p_ucell->latvec.e32 = 0.0; this->p_ucell->latvec.e33 = 10.0;
            this->p_ucell->GT = this->p_ucell->latvec.Inverse();
            this->p_ucell->G = this->p_ucell->GT.Transpose();
            this->p_ucell->GGT = this->p_ucell->G * this->p_ucell->GT;
            this->p_ucell->tpiba = 2.0 * M_PI / this->p_ucell->lat0;
            this->p_ucell->tpiba2 = this->p_ucell->tpiba * this->p_ucell->tpiba;
            // atom
            this->p_ucell->atom_label.shrink_to_fit();
            this->p_ucell->atom_label.resize(1);
            this->p_ucell->atom_label[0] = "Si";
            // atom properties
            this->p_ucell->nat = 1;
            this->p_ucell->ntype = 1;
            this->p_ucell->atoms = new Atom[1];
            this->p_ucell->set_atom_flag = true;
            this->p_ucell->atoms[0].label = "Si";
            this->p_ucell->atoms[0].mass = 28.0855;
            this->p_ucell->atoms[0].na = 1;
            this->p_ucell->atoms[0].angle1.resize(1, 0.0);
            this->p_ucell->atoms[0].angle2.resize(1, 0.0);
            // atom position
            this->p_ucell->atoms[0].tau.resize(1, {0.0, 0.0, 0.0});
            this->p_ucell->atoms[0].taud.resize(1, {0.25, 0.25, 0.25});
            this->p_ucell->atoms[0].mbl.resize(1, {0, 0, 0});
            // atom pseudopotential
            this->p_ucell->pseudo_fn.shrink_to_fit();
            this->p_ucell->pseudo_fn.resize(1);
            this->p_ucell->pseudo_fn[0] = "Si_NCSR_ONCVPSP_v0.5_dojo.upf";
            this->p_ucell->natomwfc = 4;
            this->p_ucell->atoms[0].ncpp.nchi = 2;
            this->p_ucell->atoms[0].ncpp.els = std::vector<std::string>(2, "");
            this->p_ucell->atoms[0].ncpp.mesh = 11;
            this->p_ucell->atoms[0].ncpp.msh = 11;
            this->p_ucell->atoms[0].ncpp.lmax = 2;
            //if(this->p_ucell->atoms[0].ncpp.rab != nullptr) delete[] this->p_ucell->atoms[0].ncpp.rab;
            this->p_ucell->atoms[0].ncpp.rab = std::vector<double>(11, 0.0);
            for(int i = 0; i < 11; ++i) { this->p_ucell->atoms[0].ncpp.rab[i] = 0.01;
}
            //if(this->p_ucell->atoms[0].ncpp.r != nullptr) delete[] this->p_ucell->atoms[0].ncpp.r;
            this->p_ucell->atoms[0].ncpp.r = std::vector<double>(11, 0.0);
            for(int i = 0; i < 11; ++i) { this->p_ucell->atoms[0].ncpp.r[i] = 0.01*i;
}
            this->p_ucell->atoms[0].ncpp.chi.create(2, 11);
            for(int i = 0; i < 2; ++i) { for(int j = 0; j < 11; ++j) { this->p_ucell->atoms[0].ncpp.chi(i, j) = 0.01;
}
}
            //if(this->p_ucell->atoms[0].ncpp.lchi != nullptr) delete[] this->p_ucell->atoms[0].ncpp.lchi;
            this->p_ucell->atoms[0].ncpp.lchi = std::vector<int>(2, 0);
            this->p_ucell->atoms[0].ncpp.lchi[0] = 0;
            this->p_ucell->atoms[0].ncpp.lchi[1] = 1;
            this->p_ucell->lmax_ppwf = 1;
            this->p_ucell->atoms[0].ncpp.oc = std::vector<double>(2, 0.0);
            this->p_ucell->atoms[0].ncpp.oc[0] = 1.0;
            this->p_ucell->atoms[0].ncpp.oc[1] = 1.0;

            this->p_ucell->atoms[0].ncpp.has_so = false;
            this->p_ucell->atoms[0].ncpp.jchi = std::vector<double>(2, 0.0);
            this->p_ucell->atoms[0].ncpp.jchi[0] = 0.5;
            this->p_ucell->atoms[0].ncpp.jchi[1] = 1.5;
            // atom numerical orbital
            this->p_ucell->lmax = 2;
            p_ucell->orbital_fn.shrink_to_fit();
            p_ucell->orbital_fn.resize(1);
            this->p_ucell->orbital_fn[0] = "Si_gga_8au_60Ry_2s2p1d.orb";
            this->p_ucell->atoms[0].nwl = 2;
            this->p_ucell->atoms[0].l_nchi.resize(3);
            this->p_ucell->atoms[0].l_nchi[0] = 2;
            this->p_ucell->atoms[0].l_nchi[1] = 2;
            this->p_ucell->atoms[0].l_nchi[2] = 1;

            
            // can support function PW_Basis::getfftixy2is
            this->p_pw_wfc->nks = 1;
            this->p_pw_wfc->npwk_max = 1;
            if(this->p_pw_wfc->npwk != nullptr) { delete[] this->p_pw_wfc->npwk;
}
            this->p_pw_wfc->npwk = new int[1];
            this->p_pw_wfc->npwk[0] = 1;
            this->p_pw_wfc->fftnxy = 1;
            this->p_pw_wfc->fftnz = 1;
            this->p_pw_wfc->nst = 1;
            this->p_pw_wfc->nz = 1;
            if(this->p_pw_wfc->is2fftixy != nullptr) { delete[] this->p_pw_wfc->is2fftixy;
}
            this->p_pw_wfc->is2fftixy = new int[1];
            this->p_pw_wfc->is2fftixy[0] = 0;
            if(this->p_pw_wfc->fftixy2ip != nullptr) { delete[] this->p_pw_wfc->fftixy2ip;
}
            this->p_pw_wfc->fftixy2ip = new int[1];
            this->p_pw_wfc->fftixy2ip[0] = 0;
            if(this->p_pw_wfc->igl2isz_k != nullptr) { delete[] this->p_pw_wfc->igl2isz_k;
}
            this->p_pw_wfc->igl2isz_k = new int[1];
            this->p_pw_wfc->igl2isz_k[0] = 0;
            if(this->p_pw_wfc->gcar != nullptr) { delete[] this->p_pw_wfc->gcar;
}
            this->p_pw_wfc->gcar = new ModuleBase::Vector3<double>[1];
            this->p_pw_wfc->gcar[0] = {0.0, 0.0, 0.0};
            if(this->p_pw_wfc->igl2isz_k != nullptr) { delete[] this->p_pw_wfc->igl2isz_k;
}
            this->p_pw_wfc->igl2isz_k = new int[1];
            this->p_pw_wfc->igl2isz_k[0] = 0;
            if(this->p_pw_wfc->gk2 != nullptr) { delete[] this->p_pw_wfc->gk2;
}
            this->p_pw_wfc->gk2 = new double[1];
            this->p_pw_wfc->gk2[0] = 0.0;
            this->p_pw_wfc->latvec.e11 = this->p_ucell->latvec.e11; this->p_pw_wfc->latvec.e12 = this->p_ucell->latvec.e12; this->p_pw_wfc->latvec.e13 = this->p_ucell->latvec.e13;
            this->p_pw_wfc->latvec.e21 = this->p_ucell->latvec.e21; this->p_pw_wfc->latvec.e22 = this->p_ucell->latvec.e22; this->p_pw_wfc->latvec.e23 = this->p_ucell->latvec.e23;
            this->p_pw_wfc->latvec.e31 = this->p_ucell->latvec.e31; this->p_pw_wfc->latvec.e32 = this->p_ucell->latvec.e32; this->p_pw_wfc->latvec.e33 = this->p_ucell->latvec.e33;
            this->p_pw_wfc->G = this->p_ucell->G;
            this->p_pw_wfc->GT = this->p_ucell->GT;
            this->p_pw_wfc->GGT = this->p_ucell->GGT;
            this->p_pw_wfc->lat0 = this->p_ucell->lat0;
            this->p_pw_wfc->tpiba = 2.0 * M_PI / this->p_ucell->lat0;
            this->p_pw_wfc->tpiba2 = this->p_pw_wfc->tpiba * this->p_pw_wfc->tpiba;
            if(this->p_pw_wfc->kvec_c != nullptr) { delete[] this->p_pw_wfc->kvec_c;
}
            this->p_pw_wfc->kvec_c = new ModuleBase::Vector3<double>[1];
            this->p_pw_wfc->kvec_c[0] = {0.0, 0.0, 0.0};
            if(this->p_pw_wfc->kvec_d != nullptr) { delete[] this->p_pw_wfc->kvec_d;
}
            this->p_pw_wfc->kvec_d = new ModuleBase::Vector3<double>[1];
            this->p_pw_wfc->kvec_d[0] = {0.0, 0.0, 0.0};

            this->p_pspot_vnl->lmaxkb = 1;

            this->p_kv->ik2iktot.resize(1);
            this->p_kv->ik2iktot[0] = 0;

        }
        void TearDown() override
        {
            delete this->psi_init;
            delete this->p_sf;
            delete this->p_pw_wfc;
            delete this->p_ucell;
            delete this->p_pspot_vnl;
            delete this->p_kv;
         }
};

TEST_F(PsiIntializerUnitTest, ConstructorRandom) {
    this->psi_init = new psi_initializer_random<std::complex<double>>();
    EXPECT_EQ("random", this->psi_init->method());
}

TEST_F(PsiIntializerUnitTest, ConstructorAtomic) {
    this->psi_init = new psi_initializer_atomic<std::complex<double>>();
    EXPECT_EQ("atomic", this->psi_init->method());
}

TEST_F(PsiIntializerUnitTest, ConstructorAtomicRandom) {
    this->psi_init = new psi_initializer_atomic_random<std::complex<double>>();
    EXPECT_EQ("atomic+random", this->psi_init->method());
}

TEST_F(PsiIntializerUnitTest, ConstructorNao) {
    this->psi_init = new psi_initializer_nao<std::complex<double>>();
    EXPECT_EQ("nao", this->psi_init->method());
}

TEST_F(PsiIntializerUnitTest, ConstructorNaoRandom) {
    this->psi_init = new psi_initializer_nao_random<std::complex<double>>();
    EXPECT_EQ("nao+random", this->psi_init->method());
}

TEST_F(PsiIntializerUnitTest, CastToT) {
    this->psi_init = new psi_initializer_random<std::complex<double>>();
    std::complex<double> cd = {1.0, 2.0};
    std::complex<float> cf = {1.0, 2.0};
    double d = 1.0;
    float f = 1.0;
    EXPECT_EQ(this->psi_init->template cast_to_T<std::complex<double>>(cd), cd);
    EXPECT_EQ(this->psi_init->template cast_to_T<std::complex<float>>(cd), cf);
    EXPECT_EQ(this->psi_init->template cast_to_T<double>(cd), d);
    EXPECT_EQ(this->psi_init->template cast_to_T<float>(cd), f);
}

TEST_F(PsiIntializerUnitTest, CalPsigRandom) {
    PARAM.input.init_wfc = "random";
    this->psi_init = new psi_initializer_random<std::complex<double>>();
    this->psi_init->initialize(this->p_sf, 
                               this->p_pw_wfc, 
                               this->p_ucell, 
                               this->p_kv, 
                               this->random_seed,
                               this->p_pspot_vnl,
                               GlobalV::MY_RANK);
    this->psi_init->tabulate(); // always: new, initialize, tabulate, allocate, proj_ao_onkG
    const int nbands_start = this->psi_init->nbands_start();
    const int nbasis = this->p_pw_wfc->npwk_max * PARAM.globalv.npol;
    psi::Psi<std::complex<double>>* psi = new psi::Psi<std::complex<double>>(1, nbands_start, nbasis, nbasis, true);
    this->psi_init->init_psig(psi->get_pointer(), 0);
    EXPECT_NEAR(-0.66187696761064307, psi->operator()(0,0,0).real(), 1e-4);
    delete psi;
}

TEST_F(PsiIntializerUnitTest, CalPsigAtomic) {
    PARAM.input.init_wfc = "atomic";
    this->psi_init = new psi_initializer_atomic<std::complex<double>>();
    this->psi_init->initialize(this->p_sf, 
                               this->p_pw_wfc, 
                               this->p_ucell, 
                               this->p_kv, 
                               this->random_seed,
                               this->p_pspot_vnl,
                               GlobalV::MY_RANK);
    this->psi_init->tabulate(); // always: new, initialize, tabulate, allocate, proj_ao_onkG
    const int nbands_start = this->psi_init->nbands_start();
    const int nbasis = this->p_pw_wfc->npwk_max * PARAM.globalv.npol;
    psi::Psi<std::complex<double>>* psi = new psi::Psi<std::complex<double>>(1, nbands_start, nbasis, nbasis, true);
    this->psi_init->init_psig(psi->get_pointer(), 0);
    EXPECT_NEAR(0, psi->operator()(0,0,0).real(), 1e-12);
    delete psi;
}

TEST_F(PsiIntializerUnitTest, CalPsigAtomicSoc) {
    PARAM.input.init_wfc = "atomic";
    PARAM.input.nspin = 4;
    PARAM.sys.npol = 2;
    this->p_ucell->atoms[0].ncpp.has_so = false;
    this->p_ucell->natomwfc *= 2;
    this->psi_init = new psi_initializer_atomic<std::complex<double>>();
    this->psi_init->initialize(this->p_sf, 
                               this->p_pw_wfc, 
                               this->p_ucell, 
                               this->p_kv, 
                               this->random_seed,
                               this->p_pspot_vnl,
                               GlobalV::MY_RANK);
    this->psi_init->tabulate(); // always: new, initialize, tabulate, allocate, proj_ao_onkG
    const int nbands_start = this->psi_init->nbands_start();
    const int nbasis = this->p_pw_wfc->npwk_max * PARAM.globalv.npol;
    psi::Psi<std::complex<double>>* psi = new psi::Psi<std::complex<double>>(1, nbands_start, nbasis, nbasis, true);
    this->psi_init->init_psig(psi->get_pointer(), 0);
    EXPECT_NEAR(0, psi->operator()(0,0,0).real(), 1e-12);
    PARAM.input.nspin = 1;
    PARAM.sys.npol = 1;
    this->p_ucell->atoms[0].ncpp.has_so = false;
    this->p_ucell->natomwfc /= 2;
    delete psi;
}

TEST_F(PsiIntializerUnitTest, CalPsigAtomicSocHasSo) {
    PARAM.input.init_wfc = "atomic";
    PARAM.input.nspin = 4;
    PARAM.sys.npol = 2;
    this->p_ucell->atoms[0].ncpp.has_so = true;
    this->p_ucell->natomwfc *= 2;
    this->psi_init = new psi_initializer_atomic<std::complex<double>>();
    this->psi_init->initialize(this->p_sf, 
                               this->p_pw_wfc, 
                               this->p_ucell, 
                               this->p_kv, 
                               this->random_seed,
                               this->p_pspot_vnl,
                               GlobalV::MY_RANK);
    this->psi_init->tabulate(); // always: new, initialize, tabulate, allocate, proj_ao_onkG
    const int nbands_start = this->psi_init->nbands_start();
    const int nbasis = this->p_pw_wfc->npwk_max * PARAM.globalv.npol;
    psi::Psi<std::complex<double>>* psi = new psi::Psi<std::complex<double>>(1, nbands_start, nbasis, nbasis, true);
    this->psi_init->init_psig(psi->get_pointer(), 0);
    EXPECT_NEAR(0, psi->operator()(0,0,0).real(), 1e-12);
    PARAM.input.nspin = 1;
    PARAM.sys.npol = 1;
    this->p_ucell->atoms[0].ncpp.has_so = false;
    this->p_ucell->natomwfc /= 2;
    delete psi;
}

TEST_F(PsiIntializerUnitTest, CalPsigAtomicRandom) {
    PARAM.input.init_wfc = "atomic+random";
    this->psi_init = new psi_initializer_atomic_random<std::complex<double>>();
    this->psi_init->initialize(this->p_sf, 
                               this->p_pw_wfc, 
                               this->p_ucell, 
                               this->p_kv, 
                               this->random_seed,
                               this->p_pspot_vnl,
                               GlobalV::MY_RANK);
    this->psi_init->tabulate(); // always: new, initialize, tabulate, allocate, proj_ao_onkG
    const int nbands_start = this->psi_init->nbands_start();
    const int nbasis = this->p_pw_wfc->npwk_max * PARAM.globalv.npol;
    psi::Psi<std::complex<double>>* psi = new psi::Psi<std::complex<double>>(1, nbands_start, nbasis, nbasis, true);
    this->psi_init->init_psig(psi->get_pointer(), 0);
    EXPECT_NEAR(0, psi->operator()(0,0,0).real(), 1e-12);
    delete psi;
}

TEST_F(PsiIntializerUnitTest, CalPsigNao) {
    PARAM.input.init_wfc = "nao";
    this->psi_init = new psi_initializer_nao<std::complex<double>>();
    this->psi_init->initialize(this->p_sf, 
                               this->p_pw_wfc, 
                               this->p_ucell, 
                               this->p_kv, 
                               this->random_seed,
                               this->p_pspot_vnl,
                               GlobalV::MY_RANK);
    this->psi_init->tabulate(); // always: new, initialize, tabulate, allocate, proj_ao_onkG
    const int nbands_start = this->psi_init->nbands_start();
    const int nbasis = this->p_pw_wfc->npwk_max * PARAM.globalv.npol;
    psi::Psi<std::complex<double>>* psi = new psi::Psi<std::complex<double>>(1, nbands_start, nbasis, nbasis, true);
    this->psi_init->init_psig(psi->get_pointer(), 0);
    EXPECT_NEAR(0, psi->operator()(0,0,0).real(), 1e-12);
    delete psi;
}

TEST_F(PsiIntializerUnitTest, CalPsigNaoRandom) {
    PARAM.input.init_wfc = "nao+random";
    this->psi_init = new psi_initializer_nao_random<std::complex<double>>();
    this->psi_init->initialize(this->p_sf, 
                               this->p_pw_wfc, 
                               this->p_ucell, 
                               this->p_kv, 
                               this->random_seed,
                               this->p_pspot_vnl,
                               GlobalV::MY_RANK);
    this->psi_init->tabulate(); // always: new, initialize, tabulate, allocate, proj_ao_onkG
    const int nbands_start = this->psi_init->nbands_start();
    const int nbasis = this->p_pw_wfc->npwk_max * PARAM.globalv.npol;
    psi::Psi<std::complex<double>>* psi = new psi::Psi<std::complex<double>>(1, nbands_start, nbasis, nbasis, true);
    this->psi_init->init_psig(psi->get_pointer(), 0);
    EXPECT_NEAR(0, psi->operator()(0,0,0).real(), 1e-12);
    delete psi;
}

TEST_F(PsiIntializerUnitTest, CalPsigNaoSoc) {
    PARAM.input.init_wfc = "nao";
    PARAM.input.nspin = 4;
    PARAM.sys.npol = 2;
    this->p_ucell->atoms[0].ncpp.has_so = false;
    PARAM.sys.domag = false;
    PARAM.sys.domag_z = false;
    this->psi_init = new psi_initializer_nao<std::complex<double>>();
    this->psi_init->initialize(this->p_sf, 
                               this->p_pw_wfc, 
                               this->p_ucell, 
                               this->p_kv, 
                               this->random_seed,
                               this->p_pspot_vnl,
                               GlobalV::MY_RANK);
    this->psi_init->tabulate(); // always: new, initialize, tabulate, allocate, proj_ao_onkG
    const int nbands_start = this->psi_init->nbands_start();
    const int nbasis = this->p_pw_wfc->npwk_max * PARAM.globalv.npol;
    psi::Psi<std::complex<double>>* psi = new psi::Psi<std::complex<double>>(1, nbands_start, nbasis, nbasis, true);
    this->psi_init->init_psig(psi->get_pointer(), 0);
    EXPECT_NEAR(0, psi->operator()(0,0,0).real(), 1e-12);
    delete psi;
}

TEST_F(PsiIntializerUnitTest, CalPsigNaoSocHasSo) {
    PARAM.input.init_wfc = "nao";
    PARAM.input.nspin = 4;
    PARAM.sys.npol = 2;
    this->p_ucell->atoms[0].ncpp.has_so = true;
    PARAM.sys.domag = false;
    PARAM.sys.domag_z = false;
    this->psi_init = new psi_initializer_nao<std::complex<double>>();
    this->psi_init->initialize(this->p_sf, 
                               this->p_pw_wfc, 
                               this->p_ucell, 
                               this->p_kv, 
                               this->random_seed,
                               this->p_pspot_vnl,
                               GlobalV::MY_RANK);
    this->psi_init->tabulate(); // always: new, initialize, tabulate, allocate, proj_ao_onkG
    const int nbands_start = this->psi_init->nbands_start();
    const int nbasis = this->p_pw_wfc->npwk_max * PARAM.globalv.npol;
    psi::Psi<std::complex<double>>* psi = new psi::Psi<std::complex<double>>(1, nbands_start, nbasis, nbasis, true);
    this->psi_init->init_psig(psi->get_pointer(), 0);
    EXPECT_NEAR(0, psi->operator()(0,0,0).real(), 1e-12);
    delete psi;
}

TEST_F(PsiIntializerUnitTest, CalPsigNaoSocHasSoDOMAG) {
    PARAM.input.init_wfc = "nao";
    PARAM.input.nspin = 4;
    PARAM.sys.npol = 2;
    this->p_ucell->atoms[0].ncpp.has_so = true;
    PARAM.sys.domag = true;
    PARAM.sys.domag_z = false;
    this->psi_init = new psi_initializer_nao<std::complex<double>>();
    this->psi_init->initialize(this->p_sf, 
                               this->p_pw_wfc, 
                               this->p_ucell, 
                               this->p_kv, 
                               this->random_seed,
                               this->p_pspot_vnl,
                               GlobalV::MY_RANK);
    this->psi_init->tabulate(); // always: new, initialize, tabulate, allocate, proj_ao_onkG
    const int nbands_start = this->psi_init->nbands_start();
    const int nbasis = this->p_pw_wfc->npwk_max * PARAM.globalv.npol;
    psi::Psi<std::complex<double>>* psi = new psi::Psi<std::complex<double>>(1, nbands_start, nbasis, nbasis, true);
    this->psi_init->init_psig(psi->get_pointer(), 0);
    EXPECT_NEAR(0, psi->operator()(0,0,0).real(), 1e-12);
    delete psi;
}

int main(int argc, char** argv)
{

#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}