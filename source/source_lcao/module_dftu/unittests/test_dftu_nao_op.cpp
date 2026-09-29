#include "gtest/gtest.h"
#include <chrono>

// mock of DFTU
#include "../dftu_nao_op.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "source_estate/module_dm/density_matrix.h"

Plus_U_Base dftu;

//---------------------------------------
// Unit test of DFTU_onsite operator (dftu_nao_op.cpp).
// Tests constructHR for d2d and d2cd variants.
//---------------------------------------

int test_size = 5;
int test_nw = 5;

class DFTUTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
#ifdef __MPI
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

        ucell.ntype = 1;
        ucell.nat = test_size;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.iat2it = new int[ucell.nat];
        ucell.iat2ia = new int[ucell.nat];
        ucell.atoms[0].tau.resize(ucell.nat);
        ucell.lat0 = 1.0;
        ucell.itia2iat.create(ucell.ntype, ucell.nat);
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2it[iat] = 0;
            ucell.iat2ia[iat] = iat;
            ucell.atoms[0].tau[iat] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
            ucell.itia2iat(0, iat) = iat;
        }
        ucell.atoms[0].na = test_size;
        ucell.atoms[0].nw = test_nw;
        ucell.atoms[0].iw2l.resize(test_nw);
        ucell.atoms[0].iw2m.resize(test_nw);
        ucell.atoms[0].iw2n.resize(test_nw);
        for (int iw = 0; iw < test_nw; ++iw)
        {
            ucell.atoms[0].iw2l[iw] = 2;
            ucell.atoms[0].iw2m[iw] = 0;
            ucell.atoms[0].iw2n[iw] = 0;
        }
        ucell.set_iat2iwt(1);
        init_parav();
        HR = new hamilt::HContainer<double>(ucell, paraV);

        dftu.occmat().data().resize(test_size);
        for (int iat = 0; iat < test_size; iat++)
        {
            dftu.occmat().data()[iat].resize(3);
            for (int l = 0; l < 3; l++)
            {
                dftu.occmat().data()[iat][l].resize(2);
                dftu.occmat().data()[iat][l][0].create(2 * l + 1, 2 * l + 1);
                dftu.occmat().data()[iat][l][1].create(2 * l + 1, 2 * l + 1);
            }
        }
        dftu.u_current = {U_test};
        dftu.l_channel = {orbital_c_test};
    }

    void TearDown() override
    {
        delete HR;
        delete paraV;
        delete[] ucell.atoms;
    }

    double occ_mat_c(int iat, int spin, int icc) const
    {
        return dftu.occmat().data()[iat][2][spin].c[icc];
    }

#ifdef __MPI
    void init_parav()
    {
        int nb = 10;
        int global_row = test_size * test_nw;
        int global_col = test_size * test_nw;
        std::ofstream ofs_running;
        paraV = new Parallel_Orbitals();
        paraV->init(global_row, global_col, nb, MPI_COMM_WORLD);
        paraV->set_atomic_trace(ucell.get_iat2iwt(), test_size, global_row);
    }
#else
    void init_parav()
    {
    }
#endif

    UnitCell ucell;
    hamilt::HContainer<double>* HR;
    Parallel_Orbitals* paraV;
    TwoCenterIntegrator intor_;

    int dsize;
    int my_rank = 0;
    double U_test = 1.0;
    int orbital_c_test = 2;
    double onsite_radius_test = 1.0;
};

TEST_F(DFTUTest, constructHRd2d)
{
    const int nspin = 1;
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<double> hsk(paraV, true);
    hsk.set_zero_hk();
    Grid_Driver gd(0, 0);
    const double factor = 1.0 / test_nw / test_nw / test_size / test_size;
    module_dm::DensityMatrix<double, double> dm(paraV, 1);
    dm.init_dmr(*HR);
    for (int i = 0; i < paraV->nrow; i++)
    {
        for (int j = 0; j < paraV->ncol; j++)
        {
            dm.set_dmk(1, 0, i, j, factor);
        }
    }
    dm.cal_dmr(-1);
    // reset HR
    for (int i = 0; i < HR->get_nnr(); i++)
    {
        HR->get_wrapper()[i] = 0.0;
    }
    hamilt::DFTU_onsite<hamilt::OperatorLCAO<double, double>>
        op(&hsk, kvec_d_in, HR, ucell, &gd, &intor_, {1.0}, &dftu, nspin, onsite_radius_test, &dm);
    op.contributeHR();
    for (int iat = 0; iat < test_size; iat++)
    {
        for (int icc = 0; icc < 25; icc++)
        {
            EXPECT_NEAR(occ_mat_c(iat, 0, icc), 0.5, 1e-10);
        }
    }
    for (int iap = 0; iap < HR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = HR->get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        std::vector<int> indexes1 = paraV->get_indexes_row(iat1);
        std::vector<int> indexes2 = paraV->get_indexes_col(iat2);
        int nwt = indexes1.size() * indexes2.size();
        for (int i = 0; i < nwt; ++i)
        {
            EXPECT_NEAR(tmp.get_pointer(0)[i], -10.0 * test_size, 1e-10);
        }
    }
    op.contributeHk(0);
    double* hk = hsk.get_hk();
    for (int i = 0; i < paraV->get_row_size() * paraV->get_col_size(); ++i)
    {
        EXPECT_NEAR(hk[i], -10.0 * test_size, 1e-10);
    }
}

TEST_F(DFTUTest, constructHRd2cd)
{
    const int nspin = 2;
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(2, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<std::complex<double>> hsk(paraV, true);
    hsk.set_zero_hk();
    Grid_Driver gd(0, 0);
    const double factor = 0.5 / test_nw / test_nw / test_size / test_size;
    std::vector<ModuleBase::Vector3<double>> kvec_d_dm(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    module_dm::DensityMatrix<std::complex<double>, double> dm(paraV, 2, kvec_d_dm, 1);
    dm.init_dmr(*HR);
    for (int is = 1; is <= 2; ++is)
    {
        for (int i = 0; i < paraV->nrow; i++)
        {
            for (int j = 0; j < paraV->ncol; j++)
            {
                dm.set_dmk(is, 0, i, j, std::complex<double>(factor, 0.0));
            }
        }
    }
    dm.cal_dmr(-1);
    // reset HR
    for (int i = 0; i < HR->get_nnr(); i++)
    {
        HR->get_wrapper()[i] = 0.0;
    }
    hamilt::DFTU_onsite<hamilt::OperatorLCAO<std::complex<double>, double>>
        op(&hsk, kvec_d_in, HR, ucell, &gd, &intor_, {1.0}, &dftu, nspin, onsite_radius_test, &dm);
    op.contributeHR();
    for (int iat = 0; iat < test_size; iat++)
    {
        for (int icc = 0; icc < 25; icc++)
        {
            EXPECT_NEAR(occ_mat_c(iat, 0, icc), 0.5, 1e-10);
        }
    }
    for (int iap = 0; iap < HR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = HR->get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        std::vector<int> indexes1 = paraV->get_indexes_row(iat1);
        std::vector<int> indexes2 = paraV->get_indexes_col(iat2);
        int nwt = indexes1.size() * indexes2.size();
        for (int i = 0; i < nwt; ++i)
        {
            EXPECT_NEAR(tmp.get_pointer(0)[i], -10.0 * test_size, 1e-10);
        }
    }
    op.contributeHk(0);
    std::complex<double>* hk = hsk.get_hk();
    for (int i = 0; i < paraV->get_row_size() * paraV->get_col_size(); ++i)
    {
        EXPECT_NEAR(hk[i].real(), -10.0 * test_size, 1e-10);
        EXPECT_NEAR(hk[i].imag(), 0.0, 1e-10);
    }
    op.contributeHR();
    for (int iat = 0; iat < test_size; iat++)
    {
        for (int icc = 0; icc < 25; icc++)
        {
            EXPECT_NEAR(occ_mat_c(iat, 1, icc), 0.5, 1e-10);
        }
    }
}

int main(int argc, char** argv)
{
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
